// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "element.h"
#include "refmap.h"
#include "global.h"
#include "forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    bool Node::is_constrained_vertex() const
    {
      assert(type == HERMES_TYPE_VERTEX);
      return ref <= 3 && !bnd;
    }

    void Node::ref_element(Element* e)
    {
      if (type == HERMES_TYPE_EDGE)
      {
        // store the element pointer in a free slot of 'elem'
        if (elem[0] == nullptr)
          elem[0] = e;
        else
        {
          if (elem[1] == nullptr)
            elem[1] = e;
          else
            throw Hermes::Exceptions::Exception(false, "No free slot 'elem'");
        }
      }
      ref++;
    }

    void Node::unref_element(HashTable* ht, Element* e)
    {
      if (type == HERMES_TYPE_VERTEX)
      {
        if (!--ref) ht->remove_vertex_node(id);
      }
      else
      {
        // remove the element from the array 'elem'
        if (elem[0] == e) elem[0] = nullptr;
        else if (elem[1] == e) elem[1] = nullptr;

        if (!--ref) ht->remove_edge_node(id);
      }
    }

    void Element::ref_all_nodes()
    {
      for (unsigned int i = 0; i < nvert; i++)
      {
        vn[i]->ref_element();
        en[i]->ref_element(this);
      }
      this->calc_area();
      this->calc_diameter();
    }

    int Element::get_edge_orientation(int ie) const
    {
      return (this->vn[ie]->id < this->vn[this->next_vert(ie)]->id) ? 0 : 1;
    }

    void Element::unref_all_nodes(HashTable* ht)
    {
      for (unsigned int i = 0; i < nvert; i++)
      {
        vn[i]->unref_element(ht);
        en[i]->unref_element(ht, this);
      }
      this->calc_area();
      this->calc_diameter();
    }

    Element::Element() : visited(false), area(0.0), diameter(0.0), center_set(false)
    {
    };
    int Element::next_vert(int i) const
    {
      return (i < (int)nvert - 1) ? i + 1 : 0;
    }

    int Element::prev_vert(int i) const
    {
      return (i > 0) ? i - 1 : nvert - 1;
    }

    bool Element::hsplit() const
    {
      if (active)
        return false;
      return sons[0] != nullptr;
    }

    bool Element::vsplit() const
    {
      if (active)
        return false;
      return sons[2] != nullptr;
    }

    bool Element::bsplit() const
    {
      if (active)
        return false;
      return sons[0] != nullptr && sons[2] != nullptr;
    }

    Element* Element::get_neighbor(int ie) const
    {
      Element** elem = en[ie]->elem;
      if (elem[0] == this)
        return elem[1];
      if (elem[1] == this)
        return elem[0];
      assert(0);
      return nullptr;
    }

    void Element::calc_area(bool precise_for_curvature)
    {
      // First some basic arithmetics.
      double ax, ay, bx, by;
      ax = vn[1]->x - vn[0]->x;
      ay = vn[1]->y - vn[0]->y;
      bx = vn[2]->x - vn[0]->x;
      by = vn[2]->y - vn[0]->y;

      this->area = 0.5*(ax*by - ay*bx);
      if (is_quad())
      {
        ax = bx; ay = by;
        bx = vn[3]->x - vn[0]->x;
        by = vn[3]->y - vn[0]->y;

        this->area = area + 0.5*(ax*by - ay*bx);
      }

      // Either the basic approximation is fine.
      if (!this->is_curved() || !precise_for_curvature)
        return;
      // Or we want to capture the curvature precisely.
      else
      {
        // Utility data.
        RefMap refmap_curv;
        RefMap refmap_straight;
        double3* tan;

        double x_center, y_center;
        this->get_center(x_center, y_center);

        for (int isurf = 0; isurf < this->nvert; isurf++)
        {
          // 0 - prepare data structures.
          int eo = g_quad_2d_std.get_edge_points(isurf, this->get_mode() == HERMES_MODE_TRIANGLE ? g_max_tri : g_max_quad, this->get_mode());
          int np = g_quad_2d_std.get_num_points(eo, this->get_mode());
          double* x_curv = new double[np];
          double* y_curv = new double[np];
          double* x_straight = new double[np];
          double* y_straight = new double[np];

          // 1 - get the x,y coordinates for the curved element.
          refmap_curv.set_active_element(this);
          Geom<double>* geometry = init_geom_surf(&refmap_curv, isurf, this->en[isurf]->marker, eo, tan);
          memcpy(x_curv, geometry->x, np*sizeof(double));
          memcpy(y_curv, geometry->y, np*sizeof(double));
          geometry->free();
          delete geometry;

          // 2. - act if there was no curvature
          CurvMap* cm_temp = this->cm;
          this->cm = nullptr;
          refmap_straight.set_active_element(this);
          geometry = init_geom_surf(&refmap_straight, isurf, this->en[isurf]->marker, eo, tan);
          memcpy(x_straight, geometry->x, np*sizeof(double));
          memcpy(y_straight, geometry->y, np*sizeof(double));
          geometry->free();
          delete geometry;

          // 3. - compare the two, get the updated area.
          double previous_distance;
          for (int i = 0; i < np; i++)
          {
            // Distance between the curved and straight edges.
            double distance_i = std::sqrt(std::pow(x_straight[i] - x_curv[i], 2.0) + std::pow(y_straight[i] - y_curv[i], 2.0));

            // Add to- or Subtract from- the area (depends on the curvature and we shall decide based on distance from the element center).
            double distance_from_center_curved = std::pow(x_center - x_curv[i], 2.0) + std::pow(y_center - y_curv[i], 2.0);
            double distance_from_center_straight = std::pow(x_center - x_straight[i], 2.0) + std::pow(y_center - y_straight[i], 2.0);
            bool add = distance_from_center_curved > distance_from_center_straight;

            // Calculate now the area delta.
            // It depends on the integration point number etc.
            double area_delta;
            if (i == 0)
            {
              double distance_along_edge = std::sqrt(std::pow(x_straight[i] - this->vn[isurf]->x, 2.0) + std::pow(y_straight[i] - this->vn[isurf]->y, 2.0));
              area_delta = distance_i * distance_along_edge * 0.5;
            }
            if (i > 0 && i < np - 1)
            {
              double distance_along_edge = std::sqrt(std::pow(x_straight[i] - x_straight[i - 1], 2.0) + std::pow(y_straight[i] - y_straight[i - 1], 2.0));
              area_delta = 0.5*(distance_i + previous_distance) * distance_along_edge;
            }
            if (i == np - 1)
            {
              double distance_along_edge = std::sqrt(std::pow(x_straight[i] - this->vn[(isurf + 1) % this->nvert]->x, 2.0) + std::pow(y_straight[i] - this->vn[(isurf + 1) % this->nvert]->y, 2.0));
              area_delta = distance_i * distance_along_edge * 0.5;
            }

            if (add)
              area += area_delta;
            else
              area -= area_delta;

            previous_distance = distance_i;
          }

          // 4. - re-add the curvature.
          this->cm = cm_temp;

          // clean up
          delete[] x_curv;
          delete[] y_curv;
          delete[] x_straight;
          delete[] y_straight;
        }
      }
    }

    void Element::get_center(double& x, double& y)
    {
      if (!this->center_set)
      {
#pragma omp critical
        {
          if (!this->center_set)
          {
            // x - coordinate
            this->x_center = this->vn[0]->x + this->vn[1]->x + this->vn[2]->x;
            this->y_center = this->vn[0]->y + this->vn[1]->y + this->vn[2]->y;
            if (this->is_quad())
            {
              this->x_center += this->vn[3]->x;
              this->x_center = this->x_center / 4.0;
              this->y_center += this->vn[3]->y;
              this->y_center = this->y_center / 4.0;
            }
            else
            {
              this->x_center = this->x_center / 3.0;
              this->y_center = this->y_center / 3.0;
            }
            x = this->x_center;
            y = this->y_center;
            this->center_set = true;
          }
        }
      }
      else
      {
        x = this->x_center;
        y = this->y_center;
        return;
      }
    }

    void Element::calc_diameter()
    {
      double max, l;
      if (is_triangle())
      {
        max = 0.0;
        for (int i = 0; i < 3; i++)
        {
          int j = next_vert(i);
          l = sqr(vn[i]->x - vn[j]->x) + sqr(vn[i]->y - vn[j]->y);
          if (l > max)
            max = l;
        }
      }
      else
      {
        max = sqr(vn[0]->x - vn[2]->x) + sqr(vn[0]->y - vn[2]->y);
        l = sqr(vn[1]->x - vn[3]->x) + sqr(vn[1]->y - vn[3]->y);
        if (l > max)
          max = l;
      }
      this->diameter = sqrt(max);
    }

    Node* get_vertex_node(Node* v1, Node* v2)
    {
      // initialize the new_ Node
      Node* newnode = new Node();
      newnode->type = HERMES_TYPE_VERTEX;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = -9999;
      newnode->p2 = -9999;
      newnode->x = (v1->x + v2->x) * 0.5;
      newnode->y = (v1->y + v2->y) * 0.5;

      return newnode;
    }

    Node* get_edge_node()
    {
      // initialize the new_ Node
      Node* newnode = new Node();
      newnode->type = HERMES_TYPE_EDGE;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = -9999;
      newnode->p2 = -9999;
      newnode->marker = 0;
      newnode->elem[0] = newnode->elem[1] = nullptr;

      return newnode;
    }
  }
}
