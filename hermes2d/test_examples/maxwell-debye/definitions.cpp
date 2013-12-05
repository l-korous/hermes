#include "definitions.h"

WeakFormHelmholtz::WeakFormHelmholtz(double eps, double omega, std::string left_bdy, std::string right_bdy) : WeakForm<double>(2)
{
  // Matrix.
  add_matrix_form(new MatrixFormHelmholtzEquation_real_real(0, 0, eps, omega));
  add_matrix_form(new MatrixFormHelmholtzEquation_real_imag(0, 1, eps, omega));
  add_matrix_form(new MatrixFormHelmholtzEquation_imag_real(1, 0, eps, omega));
  add_matrix_form(new MatrixFormHelmholtzEquation_imag_imag(1, 1, eps, omega));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_real(0, 0, eps, omega, left_bdy, false));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_imag(0, 1, eps, omega, left_bdy, false));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_real(1, 0, eps, omega, left_bdy, false));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_imag(1, 1, eps, omega, left_bdy, false));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_real(0, 0, eps, omega, right_bdy, true));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_imag(0, 1, eps, omega, right_bdy, true));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_real(1, 0, eps, omega, right_bdy, true));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_imag(1, 1, eps, omega, right_bdy, true));

  // Vector.
  add_vector_form(new VectorFormHelmholtzEquation_imag(1, eps, omega));
}

static double c_R(double omega, double eps)
{
  double omega_squared = omega * omega;
  return omega_squared * ((eps + omega_squared) / (1. + omega_squared));
}

static double c_I(double omega, double eps)
{
  double omega_squared = omega * omega;
  return omega_squared * ( omega * (eps - 1.) / (1. + omega_squared));
}

static double c_R_sqrt(double omega, double eps)
{
  std::complex<double> complex_unit(0., 1.);
  std::complex<double> fraction = 1. + (eps - 1.) / (1. - (complex_unit * omega));
  return (std::sqrt(fraction) *  omega).real();
}

static double c_I_sqrt(double omega, double eps)
{
  std::complex<double> complex_unit(0., 1.);
  std::complex<double> fraction = 1. + (eps - 1.) / (1. - (complex_unit * omega));
  return (std::sqrt(fraction) * omega).imag();
}

/* Jacobian forms */

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return int_dudx_dvdx<Real, Scalar>(n, wt, u, v) + c_R(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::clone() const 
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return c_I(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::clone() const 
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return -c_I(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::clone() const 
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return int_dudx_dvdx<Real, Scalar>(n, wt, u, v) + c_R(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::clone() const 
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag(*this);
}

// surf - real real
template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_real::matrix_form_surf(int n, double *wt,
  Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const
{
  Scalar result = -c_I_sqrt(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
  if (opposite_sign)
    result *= -1.;
  return result;
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_real::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_real::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_real::clone() const
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_real(*this);
}

// surf - real imag
template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::matrix_form_surf(int n, double *wt,
  Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const
{
  Scalar result = -c_R_sqrt(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
  if (opposite_sign)
    result *= -1.;
  return result;
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::clone() const
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag(*this);
}

// surf - imag real
template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::matrix_form_surf(int n, double *wt,
  Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const
{
  Scalar result = -c_R_sqrt(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
  if (opposite_sign)
    result *= -1.;
  return result;
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::clone() const
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real(*this);
}

// surf - imag imag
template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_imag::matrix_form_surf(int n, double *wt,
  Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const
{
  Scalar result = c_I_sqrt(omega, eps) * int_u_v<Real, Scalar>(n, wt, u, v);
  if (opposite_sign)
    result *= -1.;
  return result;
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_imag::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_imag::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_imag::clone() const
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_imag(*this);
}

/* Residual forms */

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormHelmholtzEquation_real::vector_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return Scalar(99999);
  /// Tohle se nepouziva
}

double WeakFormHelmholtz::VectorFormHelmholtzEquation_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHelmholtz::VectorFormHelmholtzEquation_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* WeakFormHelmholtz::VectorFormHelmholtzEquation_real::clone() const 
{
  return new WeakFormHelmholtz::VectorFormHelmholtzEquation_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::vector_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  Real result = Real(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * std::sin(this->n * M_PI * e->x[i]);
  return result * (1. / std::sqrt(2. * M_PI));
}

double WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * std::sin(this->n * M_PI * e->x[i]);
  return result * (1. / std::sqrt(2. * M_PI));
}

Ord WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return Ord(24);
}
 
VectorFormVol<double>* WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::clone() const 
{
  return new WeakFormHelmholtz::VectorFormHelmholtzEquation_imag(*this);
}