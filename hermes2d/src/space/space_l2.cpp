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

#include "global.h"
#include "space_l2.h"
#include "matrix.h"
#include "quad_all.h"
#include "shapeset/shapeset_l2_all.h"
#include "boundary_conditions/essential_boundary_conditions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    L2Space<Scalar>::L2Space() : Space<Scalar>()
    {
    }

    template<typename Scalar>
    void L2Space<Scalar>::init(Shapeset* shapeset, int p_init, bool assign_dofs_init)
    {
      if (shapeset == nullptr)
      {
        this->shapeset = new L2Shapeset;
        this->own_shapeset = true;
      }

      // enumerate basis functions
      if (assign_dofs_init)
      {
        // set uniform poly order in elements
        if (p_init < 0)
          throw Hermes::Exceptions::Exception("P_INIT must be >= 0 in an L2 space.");
        else
          this->set_uniform_order_internal(p_init, HERMES_ANY_INT);

        this->assign_dofs();
      }
    }

    template<typename Scalar>
    L2Space<Scalar>::L2Space(MeshSharedPtr mesh, int p_init, Shapeset* shapeset)
      : Space<Scalar>(mesh, shapeset, nullptr)
    {
      init(shapeset, p_init);
    }

    template<typename Scalar>
    L2Space<Scalar>::~L2Space()
    {

    }

    template<typename Scalar>
    void L2Space<Scalar>::copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh)
    {
      Space<Scalar>::copy(space, new_mesh);
    }

    template<typename Scalar>
    void L2Space<Scalar>::set_shapeset(Shapeset *shapeset)
    {
      if (shapeset->get_id() < 40 && shapeset->get_id() > 29)
      {
        this->shapeset = shapeset;
        this->own_shapeset = false;
      }
      else
        throw Hermes::Exceptions::Exception("Wrong shapeset type in L2Space<Scalar>::set_shapeset()");
    }

    template<typename Scalar>
    void L2Space<Scalar>::assign_bubble_dofs()
    {
      Element* e;
      this->bubble_functions_count = 0;
      for_all_active_elements(e, this->mesh)
      {
        typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
        ed->bdof = this->next_dof;
        ed->n = this->shapeset->get_num_bubbles(ed->order, e->get_mode()); //FIXME: this function might return invalid value because retrieved bubble functions for non-uniform orders might be invalid for the given order.
        this->next_dof += ed->n;
        this->bubble_functions_count += ed->n;
      }
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const
    {}

    template<typename Scalar>
    void L2Space<Scalar>::get_element_assembly_list(Element* e, AsmList<Scalar>* al) const
    {
      // add bubble functions to the assembly list
      al->cnt = 0;
      get_bubble_assembly_list(e, al);
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const
    {
      typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
      if (!ed->n) return;

      int* indices = this->shapeset->get_bubble_indices(ed->order, e->get_mode());
      for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof++)
      {
        //printf("triplet: %d, %d, %f\n", *indices, dof, 1.0);
        al->add_triplet(*indices++, dof, 1.0);
      }
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const
    {
      this->get_bubble_assembly_list(e, al);
    }

    template<typename Scalar>
    Scalar* L2Space<Scalar>::get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc)
    {
      throw Hermes::Exceptions::Exception("Method get_bc_projection() called from an L2Space.");
      return nullptr;
    }

    template<typename Scalar>
    L2MarkerWiseConstSpace<Scalar>::L2MarkerWiseConstSpace(MeshSharedPtr mesh) : L2Space<Scalar>(mesh, 0)
    {
    }

    template<typename Scalar>
    void L2MarkerWiseConstSpace<Scalar>::assign_bubble_dofs()
    {
      Element* e;
      this->bubble_functions_count = 0;
      int max_marker = 0;
      for_all_active_elements(e, this->mesh)
      {
        typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
        ed->bdof = this->next_dof + e->marker - 1;
        ed->n = 1;
        max_marker = std::max(max_marker, e->marker);
      }
      this->next_dof += max_marker;
      this->bubble_functions_count = max_marker;
    }

    template HERMES_API class L2Space<double>;
    template HERMES_API class L2Space<std::complex<double> >;

    template HERMES_API class L2MarkerWiseConstSpace<double>;
    template HERMES_API class L2MarkerWiseConstSpace<std::complex<double> >;
  }
}
