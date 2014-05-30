#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double slope)
            : ExactSolutionScalar<double>(mesh), slope(slope) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord (double x, double y) const;

  MeshFunction<double>* clone() const { return new CustomExactSolution(mesh, slope); }

  double slope;
};

/* Custom function f */

class CustomFunction: public Hermes::Hermes2DFunction<double>
{
public:
  CustomFunction(double slope)
    : Hermes::Hermes2DFunction<double>(), slope(slope) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double slope;
};

class CustomNormFormVol : public NormFormVol<double>
{
public:
  CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
  {
    this->set_area(HERMES_ANY);
  }

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, GeomVol<double> *e) const
  {
    double result = double(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * u->val[i] * v->val[i];
    return result * e->get_area(n, wt);
  }
};