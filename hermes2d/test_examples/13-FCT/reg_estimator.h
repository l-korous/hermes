#ifndef __REGESTIMATOR_H
#define __REGESTIMATOR_H


#include "hermes2d.h"
#include <list>


using namespace Hermes;
using namespace Hermes::Solvers;
using namespace Hermes::Algebra;
using namespace Hermes::Hermes2D;



//---------Gradient Reconstruction---------------

class GradientReconstructionMatForm_1 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_1(int i, int j) 
    : MatrixFormVol<double>(i, j){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
    Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const;


  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord>  **ext) const;    
  MatrixFormVol<double>* clone() const;

};


class GradientReconstructionVectorForm_1 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_1(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const ;


  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;
  VectorFormVol<double>* clone() const;

};

class GradientReconstruction_1 : public WeakForm<double>
{
public:
  GradientReconstruction_1(MeshFunctionSharedPtr<double> sln);

};


class GradientReconstructionMatForm_2 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_2(int i, int j) 
    : MatrixFormVol<double>(i, j){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
    Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const;


  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord>  **ext) const;    
  MatrixFormVol<double>* clone() const;
};


class GradientReconstructionVectorForm_2 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_2(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const ;


  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;
  VectorFormVol<double>* clone() const;

};
class GradientReconstruction_2 : public WeakForm<double>
{
public:
  GradientReconstruction_2(MeshFunctionSharedPtr<double> sln);

};




//----------------Regularity_Estimator Class


class Regularity_Estimator
{
public:

  Regularity_Estimator(double epsilon);	
  ~Regularity_Estimator();
  void free();

  int* get_smooth_elems(SpaceSharedPtr<double> new_space,double* coeff_vec,CSCMatrix<double> * mass_matrix=nullptr);
  int* get_smooth_dofs(SpaceSharedPtr<double> new_space,double* coeff_vec,CSCMatrix<double> * mass_matrix=nullptr);


protected:
  double linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c, MeshFunctionSharedPtr<double> sln);
  double linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c, MeshFunctionSharedPtr<double> sln);
  double linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c, MeshFunctionSharedPtr<double> sln);

  void smoothness_indicator(CSCMatrix<double> * mass_matrix = nullptr);

  void set_space(SpaceSharedPtr<double> new_space);

  double epsilon;
  MeshFunctionSharedPtr<double> sln;
  MeshFunctionSharedPtr<double> R_h_1;
  MeshFunctionSharedPtr<double> R_h_2;
  GradientReconstruction_1* grad_1;
  GradientReconstruction_2* grad_2;
  SpaceSharedPtr<double> space;
  SimpleVector<double> * rhs_1;
  SimpleVector<double> * rhs_2;

  int* smooth_elem_patch;
  int* smooth_dof;
  AsmList<double>* al;
};




#endif
