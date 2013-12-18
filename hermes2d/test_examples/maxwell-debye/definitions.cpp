#include "definitions.h"

WeakFormMD::WeakFormMD(double eps, double omega, double tau, std::string left_bdy, std::string right_bdy) : 
WeakForm<double>(6)
{
  // Matrix.
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      switch (i)
      {
      case 0:
        switch (j)
        {
        case 0:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 2:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 4:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      case 1:
        switch (j)
        {
        case 1:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 3:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 5:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      case 2:
        switch (j)
        {
        case 0:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 2:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      case 3:
        switch (j)
        {
        case 1:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 3:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      case 4:
        switch (j)
        {
        case 2:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 4:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      case 5:
        switch (j)
        {
        case 3:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        case 5:
          add_matrix_form(new MatrixForm(i, j, eps, omega, tau));
          break;
        }
        break;
      }
      
      switch (i)
      {
      case 0:
        switch (j)
        {
        case 2:
          add_matrix_form_surf(new MatrixFormSurfCustom(i, j, eps, omega, tau));
          this->mfsurf.back()->set_areas(Hermes::vector<std::string>(left_bdy, right_bdy));
          break;
        }
        break;
      case 1:
        switch (j)
        {
        case 3:
          add_matrix_form_surf(new MatrixFormSurfCustom(i, j, eps, omega, tau));
          this->mfsurf.back()->set_areas(Hermes::vector<std::string>(left_bdy, right_bdy));
          break;
        }
        break;
      }
    }
    add_vector_form(new VectorForm(i, eps, omega, tau));
  }
}

static double real_expr_sqrt(double omega, double eps)
{
  std::complex<double> complex_unit(0., 1.);
  std::complex<double> expr = (1. + ((eps - 1.) / (1. - (complex_unit * omega))));
  return std::sqrt(expr.real());
}

static double imag_expr_sqrt(double omega, double eps)
{
  std::complex<double> complex_unit(0., 1.);
  std::complex<double> expr = (1. + ((eps - 1.) / (1. - (complex_unit * omega))));
  return std::sqrt(expr.imag());
}

template<typename Real, typename Scalar>
Scalar value_bdy(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = (Scalar)0;
  for (int i_p = 0; i_p < n; i_p++)
    result += wt[i_p] * u->dx[i_p] * v->val[i_p] * e->nx[i_p];
  return result;
}

/* Jacobian forms */

template<typename Real, typename Scalar>
Scalar WeakFormMD::MatrixForm::matrix_form(int n, double *wt,
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  switch (this->i)
  {
  case 0:
    switch (this->j)
    {
    case 0:
      return ((1. / tau) + (eps - 1.)) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 2:
      return int_dudx_dvdx<Real, Scalar>(n, wt, u, v);
      break;
    case 4:
      return -(1/tau) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  case 1:
    switch (this->j)
    {
    case 1:
      return ((1. / tau) + (eps - 1.)) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 3:
      return int_dudx_dvdx<Real, Scalar>(n, wt, u, v);
      break;
    case 5:
      return -(1 / tau) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  case 2:
    switch (this->j)
    {
    case 0:
      return - int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 2:
      return (1. / tau) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  case 3:
    switch (this->j)
    {
    case 1:
      return -int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 3:
      return (1. / tau) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  case 4:
    switch (this->j)
    {
    case 2:
      return - (eps - 1.) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 4:
      return ((1. / tau) + 1.) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  case 5:
    switch (this->j)
    {
    case 3:
      return -(eps - 1.) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    case 5:
      return ((1. / tau) + 1.) * int_u_v<Real, Scalar>(n, wt, u, v);
      break;
    }
    break;
  }
  return Scalar(0);
}

double WeakFormMD::MatrixForm::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormMD::MatrixForm::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormMD::MatrixForm::clone() const
{
  return new WeakFormMD::MatrixForm(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormMD::MatrixFormSurfCustom::matrix_form(int n, double *wt,
  Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const
{
  switch (this->i)
  {
  case 0:
    switch (this->j)
    {
    case 2:
      return - value_bdy<Real, Scalar>(n, wt, u, v, e);
      break;
    }
    break;
  case 1:
    switch (this->j)
    {
    case 2:
      return -value_bdy<Real, Scalar>(n, wt, u, v, e);
      break;
    }
    break;
  }
  return Scalar(0);
}

double WeakFormMD::MatrixFormSurfCustom::value(int n, double *wt, Func<double> *u_ext[],
  Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormMD::MatrixFormSurfCustom::ord(int n, double *wt, Func<Ord> *u_ext[],
  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormMD::MatrixFormSurfCustom::clone() const
{
  return new WeakFormMD::MatrixFormSurfCustom(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormMD::VectorForm::vector_form(int n, double *wt,
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const 
{
  return (1. / tau) * int_u_v<Real, Scalar>(n, wt, ext[this->i], v);
}

double WeakFormMD::VectorForm::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormMD::VectorForm::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* WeakFormMD::VectorForm::clone() const
{
  return new WeakFormMD::VectorForm(*this);
}


void CustomSolution::derivatives(double x, double y, double& dx, double& dy) const
{
};

double CustomSolution::value(double x, double y) const
{
  std::complex<double> complex_unit(0., 1.);
  switch(this->eqn)
  {
  case 0:
    return (-complex_unit * omega * std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x)).real();
    break;
  case 1:
    return (-complex_unit * omega * std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x)).imag();
    break;
  case 2:
    return std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x).real();
    break;
  case 3:
    return std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x).imag();
    break;
  case 4:
    return (-( (eps - 1.) / (1. + complex_unit * omega)) * std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x)).real();
    break;
  case 5:
    return (-((eps - 1.) / (1. + complex_unit * omega)) * std::exp(real_expr_sqrt(omega, eps) * complex_unit * omega * x)).imag();
    break;
  }
}

Ord CustomSolution::ord(double x, double y) const
{
  return Hermes::Ord(20);
}

MeshFunction<double>* CustomSolution::clone() const
{
  return new CustomSolution(mesh, this->eqn, this->time, this->eps, this->omega);
}

CustomSolution::~CustomSolution()
{
}
