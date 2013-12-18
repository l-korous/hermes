#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class WeakFormMD : public WeakForm<double>
{
public:
  WeakFormMD(double eps, double omega, double tau, std::string left_bdy, std::string right_bdy);

private:
  class MatrixForm : public MatrixFormVol<double>
  {
  public:
    MatrixForm(int i, int j, double eps, double omega, double tau)
      : MatrixFormVol<double>(i, j), eps(eps), omega(omega), tau(tau) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
      Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    double eps;
    double omega;
    double tau;
  };

  class MatrixFormSurfCustom : public MatrixFormSurf<double>
  {
  public:
    MatrixFormSurfCustom(int i, int j, double eps, double omega, double tau)
      : MatrixFormSurf<double>(i, j), eps(eps), omega(omega), tau(tau) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
      Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormSurf<double>* clone() const;

    double eps;
    double omega;
    double tau;
  };

  class VectorForm : public VectorFormVol<double>
  {
  public:
    VectorForm(int i, double eps, double omega, double tau)
      : VectorFormVol<double>(i), eps(eps), omega(omega), tau(tau) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
      Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const;

    double eps;
    double omega;
    double tau;
  };

  class VectorFormHelmholtzEquation_imag : public VectorFormVol<double>
  {
  public:
    VectorFormHelmholtzEquation_imag(int i, double n, double omega)
      : VectorFormVol<double>(i), n(n), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
      Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const;

    double n;
    double omega;
  };
};


class CustomSolution : public ExactSolutionScalar<double>
{
public:
  CustomSolution(MeshSharedPtr mesh, int eqn, double time, double eps, double omega) : 
    ExactSolutionScalar<double>(mesh), eqn(eqn), time(time), eps(eps), omega(omega) {};
  ~CustomSolution();

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual double value(double x, double y) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const;

  int eqn;
  double time;
  double eps;
  double omega;
};