#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

class WeakFormTutorial : public WeakForm
{
public:
  // Problem parameters.
  double CONST_F;

  WeakFormTutorial() : WeakForm(1)
  {
    CONST_F = 2.0;

    add_matrix_form(new MatrixFormVolTutorial(0, 0));
    add_vector_form(new VectorFormVolTutorial(0));
  }

private:
  class MatrixFormVolTutorial : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolTutorial(int i, int j) : WeakForm::MatrixFormVol(i, j) 
    {
      adapt_eval = false;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolTutorial : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolTutorial(int i) : WeakForm::VectorFormVol(i)
    {
      adapt_eval = false;
    }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return static_cast<WeakFormTutorial *>(wf)->CONST_F * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };
};
