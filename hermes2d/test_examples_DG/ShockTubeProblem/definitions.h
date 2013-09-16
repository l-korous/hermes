#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::WeakFormsH1;

enum TimeSteppingType
{
  Explicit,
  Implicit
};

enum SolvedExample
{
  AdvectedCube,
  SolidBodyRotation
};

class CustomWeakForm  : public WeakForm<double>     
{
public:
  CustomWeakForm(SolvedExample solvedExample, TimeSteppingType timeSteppingType = Explicit, int explicitSchemeStep = 1);

  typedef double (*scalar_product_with_advection_direction)(double x, double y, double vx, double vy);

  class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
  {
  public:
    CustomMatrixFormVolConvection(int i, int j);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    MatrixFormVol<double>* clone() const;

    CustomWeakForm::scalar_product_with_advection_direction advection_term;
  };

  class CustomVectorFormVolConvection : public VectorFormVol<double>   
  {
  public:
    CustomVectorFormVolConvection(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    VectorFormVol<double>* clone() const;
    CustomWeakForm::scalar_product_with_advection_direction advection_term;
  };

  class CustomMatrixFormVol : public MatrixFormVol<double>   
  {
  public:
    CustomMatrixFormVol(int i, int j, double factor);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    MatrixFormVol<double>* clone() const;
    double factor;
  };

  class CustomVectorFormVol : public VectorFormVol<double>   
  {
  public:
    CustomVectorFormVol(int i, int prev, double factor);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    VectorFormVol<double>* clone() const;
    int prev;
    double factor;
  };

  class CustomMatrixFormInterface : public MatrixFormDG<double>
  {
  public:
    CustomMatrixFormInterface(int i, int j) : MatrixFormDG<double>(i, j) 
    {
    };

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;

    double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

    Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;

    CustomWeakForm::scalar_product_with_advection_direction advection_term;
  };

  class CustomVectorFormInterface : public VectorFormDG<double>
  {
  public:
    CustomVectorFormInterface(int i) : VectorFormDG<double>(i) 
    {
    };

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    VectorFormDG<double>* clone() const;

    double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

    Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;

    CustomWeakForm::scalar_product_with_advection_direction advection_term;
  };
};

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

int test();