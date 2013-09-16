#ifndef __VERTEX_BASED_DEFS
#define __VERTEX_BASED_DEFS
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::WeakFormsH1;

#pragma region Utilities

enum Algorithm
{
  Multiscale,
  pMultigrid
};

enum SolvedExample
{
  AdvectedCube,
  SolidBodyRotation,
  CircularConvection,
  MovingPeak,
  Benchmark
};

extern double upwind_flux(double u_cent, double u_neib, double a_dot_n);

extern Ord upwind_flux(Ord u_cent, Ord u_neib, double a_dot_n);

#pragma endregion

#pragma region WeakForm

class SmoothingWeakForm  : public WeakForm<double>     
{
public:
  SmoothingWeakForm(SolvedExample solvedExample, bool local, int explicitSchemeStep = 1, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0., double s = 0., double sigma = 0.);
};

class SmoothingWeakFormResidual  : public WeakForm<double>     
{
public:
  SmoothingWeakFormResidual(SolvedExample solvedExample, int explicitSchemeStep = 1, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0., double s = 0., double sigma = 0.);
};

class ExactWeakForm : public WeakForm<double>
{
public:
  ExactWeakForm(SolvedExample solvedExample, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0., double s = 0., double sigma = 0., bool matrix_only = false);
};

class FullImplicitWeakForm : public WeakForm<double>
{
public:
  FullImplicitWeakForm(SolvedExample solvedExample, int explicitSchemeStep = 1, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0.);
};

class ImplicitWeakForm : public WeakForm<double>
{
public:
  ImplicitWeakForm(SolvedExample solvedExample, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0., double s = 0., double sigma = 0.);
};

class ExplicitWeakForm  : public WeakForm<double>     
{
public:
  ExplicitWeakForm(SolvedExample solvedExample, bool add_inlet = false, std::string inlet = "", std::string outlet = "", double diffusivity = 0., double s = 0., double sigma = 0.);
};

class ErrorWeakForm  : public WeakForm<double>     
{
public:
  ErrorWeakForm();
};


#pragma endregion

#pragma region InitialCondition
class InitialConditionAdvectedCube : public ExactSolutionScalar<double>
{
public:
  InitialConditionAdvectedCube(MeshSharedPtr mesh);

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

class InitialConditionSolidBodyRotation : public ExactSolutionScalar<double>
{
public:
  InitialConditionSolidBodyRotation(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

class ExactSolutionCircularConvection : public ExactSolutionScalar<double>
{
public:
  ExactSolutionCircularConvection(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

class ExactSolutionMovingPeak : public ExactSolutionScalar<double>
{
public:
  ExactSolutionMovingPeak(MeshSharedPtr mesh, double diffusivity, double time) : ExactSolutionScalar<double>(mesh), diffusivity(diffusivity) { this->set_current_time(time); }

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;

  double get_current_time() const;
  void set_current_time(double time);
  double current_time;
  double diffusivity, x_hat, y_hat;
};

class InitialConditionBenchmark : public ExactSolutionScalar<double>
{
public:
  InitialConditionBenchmark(MeshSharedPtr mesh, double diffusivity) : ExactSolutionScalar<double>(mesh), diffusivity(diffusivity){};
  
  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
  double diffusivity;
};

class ExactSolutionBenchmark : public ExactSolutionScalar<double>
{
public:
  ExactSolutionBenchmark(MeshSharedPtr mesh, double diffusivity) : ExactSolutionScalar<double>(mesh), diffusivity(diffusivity) {};
  
  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
  double diffusivity;
};

#pragma endregion

double* merge_slns(double* solution_vector_coarse, SpaceSharedPtr<double> space_coarse, double* solution_vector_fine, SpaceSharedPtr<double> space_fine, SpaceSharedPtr<double> space_full);
Hermes::Algebra::Vector<double>* cut_off_linear_part(Hermes::Algebra::Vector<double>* src_vector, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine);

class MyErrorCalculator : public ErrorCalculator<double>
{
public:
  MyErrorCalculator(CalculatedErrorType errorType, int component_count) : ErrorCalculator<double>(errorType)
  {
    DefaultNormFormSurf<double>* form = new DefaultNormFormSurf<double>(0, 0, HERMES_L2_NORM);
    form->set_area("Outlet");
    this->add_error_form(form);
  }
  virtual ~MyErrorCalculator() {}
};

#endif
