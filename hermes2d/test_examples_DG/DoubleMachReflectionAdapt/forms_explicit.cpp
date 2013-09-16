#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "../euler_util.h"

class EulerEquationsWeakFormExplicitDoubleReflection : public WeakForm<double>
{
public:
  double kappa;
  bool fvm_only;
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> prescribed_markers;

  MeshFunctionSharedPtr<double> prev_density;
  MeshFunctionSharedPtr<double> prev_density_vel_x;
  MeshFunctionSharedPtr<double> prev_density_vel_y;
  MeshFunctionSharedPtr<double> prev_energy;

  MeshFunctionSharedPtr<double> exact_density;
  MeshFunctionSharedPtr<double> exact_density_vel_x;
  MeshFunctionSharedPtr<double> exact_density_vel_y;
  MeshFunctionSharedPtr<double> exact_energy;

  // External state.
  Hermes::vector<double> rho_ext;
  Hermes::vector<double> v1_ext;
  Hermes::vector<double> v2_ext;
  Hermes::vector<double> pressure_ext;
  Hermes::vector<double> energy_ext;

  // Fluxes for calculation.
  EulerFluxes* euler_fluxes;

  // Constructor for one inflow with different external states.
  EulerEquationsWeakFormExplicitDoubleReflection(double kappa, 
    Hermes::vector<std::string> solid_wall_markers,
    Hermes::vector<std::string> prescribed_markers,
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y,  MeshFunctionSharedPtr<double> prev_energy, 
    MeshFunctionSharedPtr<double> exact_density, MeshFunctionSharedPtr<double> exact_density_vel_x, MeshFunctionSharedPtr<double> exact_density_vel_y,  MeshFunctionSharedPtr<double> exact_energy, 
    bool fvm_only = false, int num_of_equations = 4) :
  WeakForm<double>(num_of_equations), 
    kappa(kappa), solid_wall_markers(solid_wall_markers), prescribed_markers(prescribed_markers),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), 
    exact_density(exact_density), exact_density_vel_x(exact_density_vel_x), exact_density_vel_y(exact_density_vel_y), exact_energy(exact_energy), 
    fvm_only(fvm_only), 
    euler_fluxes(new EulerFluxes(kappa))
  {
    for(int form_i = 0; form_i < 4; form_i++)
    {
      add_matrix_form(new EulerEquationsBilinearFormTime(form_i));

      add_vector_form(new EulerEquationsLinearFormTime(form_i));

      EulerEquationsVectorFormFlux* formDG = new EulerEquationsVectorFormFlux(form_i, kappa, euler_fluxes);
      add_vector_form_DG(formDG);

      EulerEquationsVectorFormBdyFlux* form_flux = new EulerEquationsVectorFormBdyFlux(form_i, exact_density, exact_density_vel_x, exact_density_vel_y, exact_energy, kappa);
      form_flux->set_areas(prescribed_markers);
      add_vector_form_surf(form_flux);

      add_vector_form_surf(new EulerEquationsVectorFormSolidWall(form_i, solid_wall_markers, kappa));

      for(int form_j = 0; form_j < 4; form_j++)
      {
        if(!fvm_only)
          add_vector_form(new EulerEquationsBilinearForm(form_i, form_j, euler_fluxes));
      }
    }

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  virtual ~EulerEquationsWeakFormExplicitDoubleReflection()
  {
    delete this->euler_fluxes;
  }

  WeakForm<double>* clone() const
  {
    EulerEquationsWeakFormExplicitDoubleReflection* wf;
    wf = new EulerEquationsWeakFormExplicitDoubleReflection(this->kappa, this->solid_wall_markers, this->prescribed_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->exact_density, this->exact_density_vel_x, this->exact_density_vel_y, this->exact_energy, this->fvm_only, this->neq);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->ext[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        wf->ext.push_back(new Solution<double>());
        wf->ext.back()->copy(this->ext[i]);
      }
      else
        wf->ext.push_back(this->ext[i]->clone());
    }

    wf->set_current_time_step(this->get_current_time_step());

    return wf;
  }

  void cloneMembers(const WeakForm<double>* otherWf)
  {
  }

  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>* *ext) const 
    {
      return int_u_v<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return int_u_v<Ord, Ord>(n, wt, u, v);
    }

    MatrixFormVol<double>* clone() const { return new EulerEquationsBilinearFormTime(this->i); }
  };

  class EulerEquationsBilinearForm : public VectorFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int i, int j, EulerFluxes* fluxes)
      : VectorFormVol<double>(i), fluxes(fluxes), j(j) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
      Func<double>* *ext) const
    {
      Func<double>* u = ext[this->j];
      double result = 0.;
      for (int point_i = 0; point_i < n;point_i++) 
      {
        double rho = ext[0]->val[point_i];
        double rho_v_x = ext[1]->val[point_i];
        double rho_v_y = ext[2]->val[point_i];
        double rho_e = ext[3]->val[point_i];

        switch(i)
        {
        case 0:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 1:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 2:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 3:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          }
          break;
        case 1:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 1:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 2:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 3:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          }
          break;
        case 2:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 1:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 2:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          case 3:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
            break;
          }
          break;
        case 3:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_0(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_0(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
            break;
          case 1:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_1(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_1(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
            break;
          case 2:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_2(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_2(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
            break;
          case 3:
            result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_3(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_3(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
            break;
          }
        }
      }

      return result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    VectorFormVol<double>* clone() const
    {
      return new EulerEquationsBilinearForm(this->i, this->j, this->fluxes);
    }

    EulerFluxes* fluxes;
    int j;
  };

  class EulerEquationsVectorFormFlux : public VectorFormDG<double>
  {
  public:
    EulerEquationsVectorFormFlux(int i, double kappa, EulerFluxes* fluxes) 
      : VectorFormDG<double>(i), num_flux(new HLLNumericalFlux(kappa)), fluxes(fluxes) 
    {
    }

    ~EulerEquationsVectorFormFlux()
    {
      delete num_flux;
    }

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double result = 0.;
      double w_L[4], w_R[4];
      for (int point_i = 0; point_i < n; point_i++)
      {
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];

        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, e->nx[point_i], e->ny[point_i]) * v->val[point_i];
      }

      return -result * wf->get_current_time_step();
    }

    VectorFormDG<double>* clone()  const
    { 
      EulerEquationsVectorFormFlux* form = new EulerEquationsVectorFormFlux(this->i, this->num_flux->kappa, this->fluxes);
      form->wf = this->wf;
      return form;
    }

    HLLNumericalFlux* num_flux;
    EulerFluxes* fluxes;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const 
    {
      return int_u_v<double, double>(n, wt, ext[this->i], v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return int_u_v<Ord, Ord>(n, wt, ext[this->i], v);
    }

    VectorFormVol<double>* clone() const { return new EulerEquationsLinearFormTime(this->i); }
  };

  class EulerEquationsVectorFormBdyFlux : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormBdyFlux(int i, MeshFunctionSharedPtr<double> exact_density, MeshFunctionSharedPtr<double> exact_density_vel_x, MeshFunctionSharedPtr<double> exact_density_vel_y,  MeshFunctionSharedPtr<double> exact_energy, double kappa)
      : VectorFormSurf<double>(i),
      num_flux(new HLLNumericalFlux(kappa)),
      exact_density(exact_density), exact_density_vel_x(exact_density_vel_x), exact_density_vel_y(exact_density_vel_y), exact_energy(exact_energy),
      kappa(kappa)
    {

    }

    ~EulerEquationsVectorFormBdyFlux()
    {
      delete this->num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;
      double w_L[4], w_R[4];

      ExactSolutionScalar<double>* sln[4];
      sln[0] = dynamic_cast<ExactSolutionScalar<double>*>(this->exact_density.get());
      sln[1] = dynamic_cast<ExactSolutionScalar<double>*>(this->exact_density_vel_x.get());
      sln[2] = dynamic_cast<ExactSolutionScalar<double>*>(this->exact_density_vel_y.get());
      sln[3] = dynamic_cast<ExactSolutionScalar<double>*>(this->exact_energy.get());

      for (int point_i = 0; point_i < n; point_i++)
      {
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        for(int k = 0; k < 4; k++)
          w_R[k] = sln[k]->value(e->x[point_i], e->y[point_i]);

        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, e->nx[point_i], e->ny[point_i]) * v->val[point_i];
      }

      return -result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* clone()  const
    {
      EulerEquationsVectorFormBdyFlux* form = new EulerEquationsVectorFormBdyFlux(*this);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double kappa;
    HLLNumericalFlux* num_flux;
    MeshFunctionSharedPtr<double> exact_density;
    MeshFunctionSharedPtr<double> exact_density_vel_x;
    MeshFunctionSharedPtr<double> exact_density_vel_y;
    MeshFunctionSharedPtr<double> exact_energy;
  };

  class EulerEquationsVectorFormSolidWall : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSolidWall(int i, Hermes::vector<std::string> markers, double kappa)
      : VectorFormSurf<double>(i), num_flux(new HLLNumericalFlux(kappa)), kappa(kappa) {set_areas(markers);}

    ~EulerEquationsVectorFormSolidWall()
    {
      delete this->num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;;
      double w_L[4], w_R[4];

      for (int point_i = 0; point_i < n; point_i++)
      {
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val[point_i];
        w_R[1] = ext[1]->val[point_i] - 2 * e->nx[point_i] * ((ext[1]->val[point_i] * e->nx[i]) + (ext[2]->val[point_i] * e->ny[i]));
        w_R[2] = ext[2]->val[point_i] - 2 * e->ny[point_i] * ((ext[1]->val[point_i] * e->nx[i]) + (ext[2]->val[point_i] * e->ny[i]));
        w_R[3] = ext[3]->val[point_i];

        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, e->nx[point_i], e->ny[point_i]) * v->val[point_i];
      }

      return -result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* clone()  const
    {
      EulerEquationsVectorFormSolidWall* form = new EulerEquationsVectorFormSolidWall(this->i, this->areas, this->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double kappa;
    HLLNumericalFlux* num_flux;
  };
};
