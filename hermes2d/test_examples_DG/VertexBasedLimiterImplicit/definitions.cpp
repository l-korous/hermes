#include "definitions.h"


static double advection_term_cube(double x, double y, double vx, double vy)
{
  return vx + vy;
}

static double advection_term_solid_body_rotation(double x, double y, double vx, double vy)
{
  return vx * (0.5 - y) + vy * (x - 0.5);
}


CustomWeakForm::CustomWeakForm(SolvedExample solvedExample, TimeSteppingType timeSteppingType, int explicitSchemeStep) : WeakForm<double>(1) 
{
  scalar_product_with_advection_direction advection_term;
  if(solvedExample == AdvectedCube)
    advection_term = advection_term_cube;
  else if(solvedExample == SolidBodyRotation)
    advection_term = advection_term_solid_body_rotation;

  if(timeSteppingType == Implicit)
  {
    CustomMatrixFormVolConvection* mf_vol = new CustomMatrixFormVolConvection(0, 0);
    mf_vol->advection_term = advection_term;
    add_matrix_form(mf_vol);

    CustomMatrixFormInterface* mf_interface = new CustomMatrixFormInterface(0, 0);
    mf_interface->advection_term = advection_term;
    add_matrix_form_DG(mf_interface);

    add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
    add_vector_form(new CustomVectorFormVol(0, 0, 1.));
  }
  else
  {
    CustomVectorFormVolConvection* vf_vol = new CustomVectorFormVolConvection(0);
    vf_vol->advection_term = advection_term;
    add_vector_form(vf_vol);

    CustomVectorFormInterface* vf_interface = new CustomVectorFormInterface(0);
    vf_interface->advection_term = advection_term;
    add_vector_form_DG(vf_interface);

    if(explicitSchemeStep == 1)
    {
      add_matrix_form(new CustomMatrixFormVol(0, 0, 1.));
      add_vector_form(new CustomVectorFormVol(0, 0, 1.));
    }
    else if (explicitSchemeStep == 2)
    {
      add_matrix_form(new CustomMatrixFormVol(0, 0, 4.));
      add_vector_form(new CustomVectorFormVol(0, 0, 1.));
      add_vector_form(new CustomVectorFormVol(0, 1, 3.));
    }
    else if (explicitSchemeStep == 3)
    {
      add_matrix_form(new CustomMatrixFormVol(0, 0, 3./2.));
      add_vector_form(new CustomVectorFormVol(0, 0, 1.));
      add_vector_form(new CustomVectorFormVol(0, 1, 1./2.));
    }
  }
}

CustomWeakForm::CustomMatrixFormVolConvection::CustomMatrixFormVolConvection(int i, int j) : MatrixFormVol<double>(i, j)
{
}

double CustomWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                            Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * this->advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return -result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                       Geom<Ord> *e, Func<Ord> **ext) const 
{
  return u->val[0] * v->dx[0];
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}

CustomWeakForm::CustomVectorFormVolConvection::CustomVectorFormVolConvection(int i) : VectorFormVol<double>(i)
{
}

double CustomWeakForm::CustomVectorFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * ext[0]->val[i] * this->advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomVectorFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
{
  return ext[0]->val[0] * v->dx[0];
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVolConvection::clone() const
{
  return new CustomVectorFormVolConvection(*this);
}

CustomWeakForm::CustomMatrixFormVol::CustomMatrixFormVol(int i, int j, double factor) : MatrixFormVol<double>(i, j), factor(factor)
{
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * u->val[i];
  return result * this->factor;
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * u->val[0];
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVol::clone() const
{
  return new CustomWeakForm::CustomMatrixFormVol(*this);
}

CustomWeakForm::CustomVectorFormVol::CustomVectorFormVol(int i, int prev, double factor) : VectorFormVol<double>(i), prev(prev), factor(factor)
{
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * ext[this->prev]->val[i];
  return result * this->factor;
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * ext[this->prev]->val[0];
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVol::clone() const
{
  return new CustomWeakForm::CustomVectorFormVol(*this);
}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
                                                        Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = this->advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * upwind_flux(0., u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * upwind_flux(u->val[i], 0., a_dot_n) * jump_v;
  }
  return result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
                                                   Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
}

MatrixFormDG<double>* CustomWeakForm::CustomMatrixFormInterface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormInterface(*this);
}

double CustomWeakForm::CustomMatrixFormInterface::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::CustomMatrixFormInterface::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}


double CustomWeakForm::CustomVectorFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
                                                        Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = this->advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    double jump_v = -v->val[i];
    result += wt[i] * upwind_flux(ext[1]->val[i], ext[1]->val_neighbor[i], a_dot_n) * jump_v;
  }
  return result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomVectorFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
                                                   Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

VectorFormDG<double>* CustomWeakForm::CustomVectorFormInterface::clone() const
{
  return new CustomWeakForm::CustomVectorFormInterface(*this);
}

double CustomWeakForm::CustomVectorFormInterface::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::CustomVectorFormInterface::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}


InitialConditionAdvectedCube::InitialConditionAdvectedCube(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
{
}

void InitialConditionAdvectedCube::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;
  dy = 0.;
}

double InitialConditionAdvectedCube::value(double x, double y) const
{
  if(x < 0. && y < 0.0 && x > -1. && y > -1.)
    return 1.0;
  else
    return 0.0;
}

Ord InitialConditionAdvectedCube::ord(double x, double y) const 
{
  return Ord(1);
}

MeshFunction<double>* InitialConditionAdvectedCube::clone() const
{
  return new InitialConditionAdvectedCube(this->mesh);

}


void InitialConditionSolidBodyRotation::derivatives(double x, double y, double& dx, double& dy) const 
{

  double radius = 0.;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  {		
    dx = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
    dy = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
  }
  else
  {			
    //cone
    x_0 = 0.5;
    y_0 = 0.25;
    radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
    if((radius< 1.0)&&(x!=x_0)) 
    { 	
      dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
      dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
    }
    else
    {
      dx=0.; dy=0.;
    }	
  }
};

double InitialConditionSolidBodyRotation::value(double x, double y) const 
{

  double result = 0.0;
  double radius;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  { 
    result = (1.0+ std::cos(M_PI*radius))/4.0;
    return result;	
  }
  //slotted cylinder
  x_0 = 0.5;
  y_0 = 0.75;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius <= 1) 
  { 	
    if(fabs((x-x_0))>= 0.025) return 1.0;
    if(y>=0.85) return 1.0;
  }	
  //cone
  x_0 = 0.5;
  y_0 = 0.25;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius<= 1.0) 
  { 	
    result = 1.0-radius;
  }	
  return result;
};

Ord InitialConditionSolidBodyRotation::ord(double x, double y) const 
{
  return Ord(10);
};
MeshFunction<double>* InitialConditionSolidBodyRotation::clone() const
{
  return new InitialConditionSolidBodyRotation(this->mesh);
}