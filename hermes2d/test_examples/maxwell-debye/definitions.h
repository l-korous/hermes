#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class WeakFormHelmholtz : public WeakForm<double>
{
public:
  WeakFormHelmholtz(double eps, double omega, std::string left_bdy, std::string right_bdy);

private:
  class MatrixFormHelmholtzEquation_real_real : public MatrixFormVol<double>
  {
  public:
    MatrixFormHelmholtzEquation_real_real(int i, int j, double eps, double omega)
      : MatrixFormVol<double>(i, j), eps(eps), omega(omega) {};

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
  };

  class MatrixFormHelmholtzEquation_real_imag : public MatrixFormVol<double>
  {
  public:
    MatrixFormHelmholtzEquation_real_imag(int i, int j, double eps, double omega)
      : MatrixFormVol<double>(i, j), eps(eps), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const;

    double eps;
    double omega;
  };

  class MatrixFormHelmholtzEquation_imag_real : public MatrixFormVol<double>
  {
  public:
    MatrixFormHelmholtzEquation_imag_real(int i, int j, double eps, double omega)
      : MatrixFormVol<double>(i, j), eps(eps), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const;

    double eps;
    double omega;
  };

  class MatrixFormHelmholtzEquation_imag_imag : public MatrixFormVol<double>
  {
  public:
    MatrixFormHelmholtzEquation_imag_imag(int i, int j, double eps, double omega)
      : MatrixFormVol<double>(i, j), eps(eps), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    double eps;
    double omega;
  };

  class MatrixFormSurfHelmholtz_real_real : public MatrixFormSurf<double>
  {
  public:
    MatrixFormSurfHelmholtz_real_real(int i, int j, double eps, double omega, std::string area, bool opposite_sign)
      : MatrixFormSurf<double>(i, j), eps(eps), omega(omega), opposite_sign(opposite_sign) { this->set_area(area); };

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;
    MatrixFormSurf<double>* clone() const;

    bool opposite_sign;

    double eps;
    double omega;
  };

  class MatrixFormSurfHelmholtz_real_imag : public MatrixFormSurf<double>
  {
  public:
    MatrixFormSurfHelmholtz_real_imag(int i, int j, double eps, double omega, std::string area, bool opposite_sign)
      : MatrixFormSurf<double>(i, j), eps(eps), omega(omega), opposite_sign(opposite_sign) { this->set_area(area); };

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;
    MatrixFormSurf<double>* clone() const;

    bool opposite_sign;

    double eps;
    double omega;
  };

  class MatrixFormSurfHelmholtz_imag_real : public MatrixFormSurf<double>
  {
  public:
    MatrixFormSurfHelmholtz_imag_real(int i, int j, double eps, double omega, std::string area, bool opposite_sign)
      : MatrixFormSurf<double>(i, j), eps(eps), omega(omega), opposite_sign(opposite_sign) { this->set_area(area); };

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormSurf<double>* clone() const;

    bool opposite_sign;

    double eps;
    double omega;
  };

  class MatrixFormSurfHelmholtz_imag_imag : public MatrixFormSurf<double>
  {
  public:
    MatrixFormSurfHelmholtz_imag_imag(int i, int j, double eps, double omega, std::string area, bool opposite_sign)
      : MatrixFormSurf<double>(i, j), eps(eps), omega(omega), opposite_sign(opposite_sign) { this->set_area(area); };

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
      Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormSurf<double>* clone() const;

    bool opposite_sign;

    double eps;
    double omega;
  };

  class VectorFormHelmholtzEquation_real : public VectorFormVol<double>
  {
  public:
    VectorFormHelmholtzEquation_real(int i, double eps, double omega)
      : VectorFormVol<double>(i), eps(eps), omega(omega) {};

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
