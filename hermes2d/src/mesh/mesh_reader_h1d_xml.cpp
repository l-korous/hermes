// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include "mesh.h"
#include "mesh_reader_h1d_xml.h"

using namespace std;

namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderH1DXML::MeshReaderH1DXML()
    {
    }

    MeshReaderH1DXML::~MeshReaderH1DXML()
    {
    }

    void MeshReaderH1DXML::load(const char *filename, MeshSharedPtr mesh)
    {
      mesh->free();

      try
      {
        ::xml_schema::flags parsing_flags = 0;
        if (!this->validate)
          parsing_flags = xml_schema::flags::dont_validate;

        std::auto_ptr<XMLMesh1D::mesh> parsed_xml_mesh(XMLMesh1D::mesh_(filename, parsing_flags));

        // Variables //
        unsigned int variables_count = parsed_xml_mesh->variables().present() ? parsed_xml_mesh->variables()->var().size() : 0;
        std::map<std::string, double> variables;
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
#ifdef _MSC_VER
          variables.insert(std::pair<std::string, double>((std::string)parsed_xml_mesh->variables()->var().at(variables_i).name(), (double&&)parsed_xml_mesh->variables()->var().at(variables_i).value()));
#else
          variables.insert(std::pair<std::string, double>((std::string)parsed_xml_mesh->variables()->var().at(variables_i).name(), parsed_xml_mesh->variables()->var().at(variables_i).value()));
#endif

        // Vertices //
        int vertices_count = parsed_xml_mesh->v().size();

        // Initialize mesh.
        int size = HashTable::H2D_DEFAULT_HASH_SIZE;
        while (size < 8 * vertices_count)
          size *= 2;
        mesh->init(size);

        double a = std::numeric_limits<double>::infinity();
        double b = -std::numeric_limits<double>::infinity();

        // Create top-level vertex nodes.
        for (int vertices_i = 0; vertices_i < 2 * vertices_count; vertices_i++)
        {
          Node* node = mesh->nodes.add();
          assert(node->id == vertices_i);
          node->ref = TOP_LEVEL_REF;
          node->type = HERMES_TYPE_VERTEX;
          node->bnd = 0;
          node->p1 = node->p2 = -1;
          node->next_hash = nullptr;

          // variables matching.
          std::string x = parsed_xml_mesh->v().at(vertices_i % vertices_count).x();
          double x_value;

          // variables lookup.
          bool x_found = false;
          if (variables.find(x) != variables.end())
          {
            x_value = variables.find(x)->second;
            x_found = true;
          }

          // test of value if no variable found.
          if (!x_found)
          if (std::strtod(x.c_str(), nullptr) != 0.0)
            x_value = std::strtod(x.c_str(), nullptr);
          else
          {
            // This is a hard part, to find out if it is really zero.
            int dot_position = strchr(x.c_str(), '.') == nullptr ? -1 : strchr(x.c_str(), '.') - x.c_str();
            for (int i = 0; i < dot_position; i++)
            if (strncmp(x.c_str() + i, "0", 1) != 0)
              this->warn("Probably wrong syntax in the x coordinate of vertex no. %i.", vertices_i % vertices_count + 1);
            for (int i = dot_position + 1; i < x.length(); i++)
            if (strncmp(x.c_str() + i, "0", 1) != 0)
              this->warn("Probably wrong syntax in the x coordinate of vertex no. %i.", vertices_i % vertices_count + 1);
            x_value = std::strtod(x.c_str(), nullptr);
          }

          // assignment.
          node->x = x_value;
          if (x_value > b)
            b = x_value;
          if (x_value < a)
            a = x_value;

          if (vertices_i < vertices_count)
            node->y = 0;
          else
            node->y = 1;
        }
        mesh->ntopvert = 2 * vertices_count;

        Node* node;
        for_all_nodes(node, mesh)
        if (node->y == 0)
          node->y = 0;
        else
          node->y = (b - a) / 100;

        // Elements //
        mesh->nbase = mesh->nactive = mesh->ninitial = vertices_count - 1;

        Element* e;
        for (int element_i = 0; element_i < vertices_count - 1; element_i++)
        {
          mesh->element_markers_conversion.insert_marker("H1DMarker");

          int element_marker;
          if (parsed_xml_mesh->v().at(element_i % vertices_count).m().present())
          {
            mesh->element_markers_conversion.insert_marker(parsed_xml_mesh->v().at(element_i % vertices_count).m().get());
            element_marker = mesh->element_markers_conversion.get_internal_marker(parsed_xml_mesh->v().at(element_i % vertices_count).m().get()).marker;
          }
          else
            element_marker = mesh->element_markers_conversion.get_internal_marker("H1DMarker").marker;

          e = mesh->create_quad(element_marker,
            &mesh->nodes[element_i],
            &mesh->nodes[element_i + 1],
            &mesh->nodes[element_i + vertices_count + 1],
            &mesh->nodes[element_i + vertices_count],
            nullptr);

          mesh->boundary_markers_conversion.insert_marker("Unused");

          node = mesh->peek_edge_node(element_i, element_i + 1);
          node->bnd = 1;
          node->marker = mesh->boundary_markers_conversion.get_internal_marker("Unused").marker;
          mesh->nodes[element_i].bnd = 1;
          mesh->nodes[element_i + 1].bnd = 1;

          node = mesh->peek_edge_node(vertices_count + element_i, vertices_count + element_i + 1);
          node->bnd = 1;
          node->marker = mesh->boundary_markers_conversion.get_internal_marker("Unused").marker;
          mesh->nodes[vertices_count + element_i].bnd = 1;
          mesh->nodes[vertices_count + element_i + 1].bnd = 1;
        }

        // Boundaries //
        Node* en;
        int v1_1 = 0;
        int v2_1 = vertices_count;
        int v1_2 = vertices_count - 1;
        int v2_2 = 2 * vertices_count - 1;

        en = mesh->peek_edge_node(v1_1, v2_1);
        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker("Left");
        int marker = mesh->boundary_markers_conversion.get_internal_marker("Left").marker;
        en->marker = marker;
        en->bnd = 1;

        en = mesh->peek_edge_node(v1_2, v2_2);
        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker("Right");
        marker = mesh->boundary_markers_conversion.get_internal_marker("Right").marker;
        en->marker = marker;
        en->bnd = 1;

        mesh->nodes[v1_1].bnd = 1;
        mesh->nodes[v2_1].bnd = 1;

        mesh->nodes[v1_2].bnd = 1;
        mesh->nodes[v2_2].bnd = 1;
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    void MeshReaderH1DXML::save(const char *filename, MeshSharedPtr mesh)
    {
      /// \todo Is this necessary? It is a valid H2D mesh afterall.
    }
  }
}
