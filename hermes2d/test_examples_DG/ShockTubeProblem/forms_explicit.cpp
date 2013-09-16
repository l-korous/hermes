#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "../euler_util.h"

class EulerEquationsWeakFormSemiImplicit : public WeakForm<double>
{
public:
  double kappa;
  bool fvm_only;
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> inlet_markers;
  Hermes::vector<std::string> outlet_markers;

  MeshFunctionSharedPtr<double> prev_density;
  MeshFunctionSharedPtr<double> prev_density_vel_x;
  MeshFunctionSharedPtr<double> prev_density_vel_y;
  MeshFunctionSharedPtr<double> prev_energy;

  // External state.
  Hermes::vector<double> rho_ext;
  Hermes::vector<double> v1_ext;
  Hermes::vector<double> v2_ext;
  Hermes::vector<double> pressure_ext;
  Hermes::vector<double> energy_ext;

  // Fluxes for calculation.
  EulerFluxes* euler_fluxes;

  // Discrete indicator in the case of Feistauer limiting.
  bool* discreteIndicator;
  int discreteIndicatorSize;

  // For cache handling.
  class EulerEquationsMatrixFormSurfSemiImplicit;
  class EulerEquationsMatrixFormSemiImplicitInletOutlet;
  bool cacheReadyDG;
  bool cacheReadySurf;
  double** P_plus_cache_DG;
  double** P_minus_cache_DG;
  double** P_plus_cache_surf;
  double** P_minus_cache_surf;

  // Constructor for one inflow with different external states.
  EulerEquationsWeakFormSemiImplicit(double kappa, 
    double rho_ext_inflow, double v1_ext_inflow, double v2_ext_inflow, double pressure_ext_inflow,
    double rho_ext_outflow, double v1_ext_outflow, double v2_ext_outflow, double pressure_ext_outflow,
    Hermes::vector<std::string> solid_wall_markers, Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers, 
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y,  MeshFunctionSharedPtr<double> prev_energy, 
    bool fvm_only = false, int num_of_equations = 4) :
  WeakForm<double>(num_of_equations), 
    kappa(kappa), 
    solid_wall_markers(solid_wall_markers), inlet_markers(inlet_markers), outlet_markers(outlet_markers), 
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), 
    fvm_only(fvm_only), 
    euler_fluxes(new EulerFluxes(kappa)), discreteIndicator(NULL)
  {
    this->rho_ext.push_back(rho_ext_inflow);
    this->v1_ext.push_back(v1_ext_inflow);
    this->v2_ext.push_back(v2_ext_inflow);
    this->pressure_ext.push_back(pressure_ext_inflow);
    this->rho_ext.push_back(rho_ext_outflow);
    this->v1_ext.push_back(v1_ext_outflow);
    this->v2_ext.push_back(v2_ext_outflow);
    this->pressure_ext.push_back(pressure_ext_outflow);
    double energy_ext_inflow = QuantityCalculator::calc_energy(rho_ext_inflow, rho_ext_inflow * v1_ext_inflow, rho_ext_inflow * v2_ext_inflow, pressure_ext_inflow, kappa);
    double energy_ext_outflow = QuantityCalculator::calc_energy(rho_ext_outflow, rho_ext_outflow * v1_ext_outflow, rho_ext_outflow * v2_ext_outflow, pressure_ext_outflow, kappa);

    P_plus_cache_DG = new double*[13];
    P_minus_cache_DG = new double*[13];
    P_plus_cache_surf = new double*[13];
    P_minus_cache_surf = new double*[13];

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
    {
      P_plus_cache_DG[coordinate_i] = new double[16];
      P_minus_cache_DG[coordinate_i] = new double[16];
      P_plus_cache_surf[coordinate_i] = new double[16];
      P_minus_cache_surf[coordinate_i] = new double[16];
    }

    for(int form_i = 0; form_i < 4; form_i++)
    {
      add_matrix_form(new EulerEquationsBilinearFormTime(form_i));

      add_vector_form(new EulerEquationsLinearFormTime(form_i));

      EulerEquationsVectorFormLinearizableSurfSemiImplicit* formDG = new EulerEquationsVectorFormLinearizableSurfSemiImplicit(form_i, kappa, euler_fluxes, &this->cacheReadyDG, this->P_plus_cache_DG, this->P_minus_cache_DG);
      add_vector_form_DG(formDG);
      
      add_vector_form_surf(new EulerEquationsVectorFormSolidWall(form_i, solid_wall_markers, kappa));
        
      for(int form_j = 0; form_j < 4; form_j++)
      {
        if(!fvm_only)
          add_vector_form(new EulerEquationsBilinearForm(form_i, form_j, euler_fluxes));
      }
    }

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  virtual ~EulerEquationsWeakFormSemiImplicit()
  {
    delete this->euler_fluxes;

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
    {
      delete [] P_plus_cache_DG[coordinate_i];
      delete [] P_minus_cache_DG[coordinate_i];
      delete [] P_plus_cache_surf[coordinate_i];
      delete [] P_minus_cache_surf[coordinate_i];
    }

    delete [] P_plus_cache_DG;
    delete [] P_minus_cache_DG;
    delete [] P_plus_cache_surf;
    delete [] P_minus_cache_surf;
  }

  void set_active_edge_state(Element** e, int isurf)
  {
    this->cacheReadySurf = false;
  }

  void set_active_DG_state(Element** e, int isurf)
  {
    this->cacheReadyDG = false;
  }

  WeakForm<double>* clone() const
  {
    EulerEquationsWeakFormSemiImplicit* wf;
    wf = new EulerEquationsWeakFormSemiImplicit(this->kappa, this->rho_ext[0], this->v1_ext[0], this->v2_ext[0], this->pressure_ext[0], this->rho_ext[1], this->v1_ext[1], this->v2_ext[1], this->pressure_ext[1], 
      this->solid_wall_markers, this->inlet_markers, this->outlet_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->fvm_only, this->neq);

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

  class EulerEquationsVectorFormLinearizableSurfSemiImplicit : public VectorFormDG<double>
  {
  public:
    EulerEquationsVectorFormLinearizableSurfSemiImplicit(int i, double kappa, EulerFluxes* fluxes, bool* cacheReady, double** P_plus_cache, double** P_minus_cache) 
      : VectorFormDG<double>(i), num_flux(new LaxFriedrichsNumericalFlux(kappa)), cacheReady(cacheReady), P_plus_cache(P_plus_cache), P_minus_cache(P_minus_cache), fluxes(fluxes) 
    {
    }

    ~EulerEquationsVectorFormLinearizableSurfSemiImplicit() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double w_L[4], w_R[4];
      double result = 0.;

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
      EulerEquationsVectorFormLinearizableSurfSemiImplicit* form = new EulerEquationsVectorFormLinearizableSurfSemiImplicit(this->i, this->num_flux->kappa, this->fluxes, this->cacheReady, this->P_plus_cache, this->P_minus_cache);
      form->wf = this->wf;
      return form;
    }

    bool* cacheReady;
    double** P_plus_cache;
    double** P_minus_cache;
    LaxFriedrichsNumericalFlux* num_flux;
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

  class EulerEquationsVectorFormSolidWall : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSolidWall(int i, Hermes::vector<std::string> markers, double kappa)
      : VectorFormSurf<double>(i), num_flux(new LaxFriedrichsNumericalFlux(kappa)), kappa(kappa) {set_areas(markers);}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double w_L[4], w_R[4];
      double result = 0.;

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
    LaxFriedrichsNumericalFlux* num_flux;
  };
};

