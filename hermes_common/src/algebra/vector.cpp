// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file vector.cpp
\brief Basic vector classes and operations.
*/
#include "common.h"
#include "matrix.h"
#include "callstack.h"

#include "solvers/linear_matrix_solver.h"
#include "solvers/interfaces/umfpack_solver.h"
#include "solvers/interfaces/superlu_solver.h"
#include "solvers/interfaces/amesos_solver.h"
#include "solvers/interfaces/petsc_solver.h"
#include "solvers/interfaces/mumps_solver.h"
#include "solvers/interfaces/aztecoo_solver.h"
#include "solvers/interfaces/paralution_solver.h"
#include "qsort.h"
#include "api.h"
#include "util/memory_handling.h"

namespace Hermes
{
  namespace Algebra
  {
    double inline real(double x)
    {
      return x;
    }

    double inline imag(double x)
    {
      return 0;
    }

    double inline real(std::complex<double> x)
    {
      return x.real();
    }

    double inline imag(std::complex<double> x)
    {
      return x.imag();
    }

    template<typename Scalar>
    Vector<Scalar>::Vector() : size(0)
    {
    }

    template<typename Scalar>
    Vector<Scalar>::Vector(unsigned int size) : size(size)
    {
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->set(i, vec->get(i));
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::set_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->set(i, vec[i]);
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::add_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i));
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i]);
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      if (!v)
        throw Exceptions::MethodNotOverridenException("Vector<Scalar>::export_to_file");

      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
      {
                                        FILE* file = fopen(filename, "w");
                                        if (!file)
                                          throw Exceptions::IOException(Exceptions::IOException::Write, filename);
                                        if (Hermes::Helpers::TypeIsReal<Scalar>::value)
                                          fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
                                        else
                                          fprintf(file, "%%%%MatrixMarket matrix coordinate complex general\n");

                                        fprintf(file, "%d 1 %d\n", this->size, this->size);

                                        for (unsigned int j = 0; j < this->size; j++)
                                        {
                                          Hermes::Helpers::fprint_coordinate_num(file, j + 1, 1, v[j], number_format);
                                          fprintf(file, "\n");
                                        }

                                        fclose(file);
      }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
      {
#ifdef WITH_MATIO
                                       size_t dims[2];
                                       dims[0] = this->size;
                                       dims[1] = 1;

                                       mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);
                                       matvar_t *matvar;

                                       // For complex.
                                       double* v_re = nullptr;
                                       double* v_im = nullptr;

                                       void* data;
                                       if (Hermes::Helpers::TypeIsReal<Scalar>::value)
                                       {
                                         data = v;
                                         matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA);
                                       }
                                       else
                                       {
                                         v_re = malloc_with_check<SimpleVector<Scalar>, double>(this->size, this);
                                         v_im = malloc_with_check<SimpleVector<Scalar>, double>(this->size, this);
                                         struct mat_complex_split_t z = { v_re, v_im };

                                         for (int i = 0; i < this->size; i++)
                                         {
                                           v_re[i] = ((std::complex<double>)(this->v[i])).real();
                                           v_im[i] = ((std::complex<double>)(this->v[i])).imag();
                                           data = &z;
                                         }
                                         matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
                                       }

                                       if (matvar)
                                       {
                                         Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
                                         Mat_VarFree(matvar);
                                       }

                                        free_with_check(v_re);
                                        free_with_check(v_im);
                                       Mat_Close(mat);

                                       if (!matvar)
                                         throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#else
                                       throw Exceptions::Exception("MATIO not included.");
#endif
      }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
      {
                                      FILE* file = fopen(filename, "w");
                                      if (!file)
                                        throw Exceptions::IOException(Exceptions::IOException::Write, filename);
                                      for (unsigned int i = 0; i < this->size; i++)
                                      {
                                        Hermes::Helpers::fprint_num(file, v[i], number_format);
                                        fprintf(file, "\n");
                                      }
                                      fclose(file);
      }
        break;

#ifdef WITH_BSON
      case EXPORT_FORMAT_BSON:
      {
                               // Init bson
                               bson bw;
                               bson_init(&bw);

                               // Matrix size.
                               bson_append_int(&bw, "size", this->size);

                               bson_append_start_array(&bw, "v");
                               for (unsigned int i = 0; i < this->size; i++)
                                 bson_append_double(&bw, "v_i", real(this->v[i]));
                               bson_append_finish_array(&bw);

                               if (!Hermes::Helpers::TypeIsReal<Scalar>::value)
                               {
                                 bson_append_start_array(&bw, "v-imag");
                                 for (unsigned int i = 0; i < this->size; i++)
                                   bson_append_double(&bw, "v_i", imag(this->v[i]));
                                 bson_append_finish_array(&bw);
                               }
                               
                               // Done.
                               bson_finish(&bw);

                               // Write to disk.
                               FILE *fpw;
                               fpw = fopen(filename, "wb");
                               const char *dataw = (const char *)bson_data(&bw);
                               fwrite(dataw, bson_size(&bw), 1, fpw);
                               fclose(fpw);

                               bson_destroy(&bw);
      }
#endif
      }
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_PLAIN_ASCII:
      {
                                      std::vector<Scalar> data;
                                      std::ifstream input(filename);
                                      if (input.bad())
                                        throw Exceptions::IOException(Exceptions::IOException::Read, filename);
                                      std::string lineData;

                                      while (getline(input, lineData))
                                      {
                                        Scalar d;
                                        std::stringstream lineStream(lineData);
                                        lineStream >> d;
                                        data.push_back(d);
                                      }

                                      this->alloc(data.size());
                                      memcpy(this->v, &data[0], sizeof(Scalar)*data.size());
      }
        break;
      case EXPORT_FORMAT_MATLAB_MATIO:
#ifdef WITH_MATIO
        mat_t    *matfp;
        matvar_t *matvar;

        matfp = Mat_Open(filename, MAT_ACC_RDONLY);

        if (!matfp)
        {
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
          return;
        }

        matvar = Mat_VarRead(matfp, var_name);
        if (matvar)
        {
          this->alloc(matvar->dims[0]);
          if (Hermes::Helpers::TypeIsReal<Scalar>::value)
            memcpy(this->v, matvar->data, sizeof(Scalar)*this->size);
          else
          {
            std::complex<double>* complex_data = malloc_with_check<SimpleVector<Scalar>, std::complex<double> >(this->size, this);
            double* real_array = (double*)((mat_complex_split_t*)matvar->data)->Re;
            double* imag_array = (double*)((mat_complex_split_t*)matvar->data)->Im;
            for (int i = 0; i < this->size; i++)
              complex_data[i] = std::complex<double>(real_array[i], imag_array[i]);
            memcpy(this->v, complex_data, sizeof(Scalar)*this->size);
            free_with_check(complex_data);
          }
        }

        Mat_Close(matfp);
        if (!matvar)
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
#else
        throw Exceptions::Exception("MATIO not included.");
#endif
        break;
      case EXPORT_FORMAT_MATRIX_MARKET:
        throw Hermes::Exceptions::MethodNotImplementedException("SimpleVector<Scalar>::import_from_file - Matrix Market");
        break;
#ifdef WITH_BSON
      case EXPORT_FORMAT_BSON:
      {
                               FILE *fpr;
                               fpr = fopen(filename, "rb");

                               // file size:
                               fseek(fpr, 0, SEEK_END);
                               int size = ftell(fpr);
                               rewind(fpr);

                               // allocate memory to contain the whole file:
                               char *datar = malloc_with_check<char>(size);
                               fread(datar, size, 1, fpr);
                               fclose(fpr);

                               bson br;
                               bson_init_finished_data(&br, datar, 0);

                               bson_iterator it;
                               bson sub;
                               bson_find(&it, &br, "size");
                               this->size = bson_iterator_int(&it);

                               this->v = malloc_with_check<SimpleVector<Scalar>, Scalar>(this->size, this);

                               bson_iterator it_coeffs;
                               bson_find(&it_coeffs, &br, "v");
                               bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                               bson_iterator_init(&it, &sub);
                               int index_coeff = 0;
                               while (bson_iterator_next(&it))
                                 this->v[index_coeff++] = bson_iterator_double(&it);

                               if (!Hermes::Helpers::TypeIsReal<Scalar>::value)
                               {
                                 bson_find(&it_coeffs, &br, "v-imag");
                                 bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                                 bson_iterator_init(&it, &sub);
                                 index_coeff = 0;
                                 while (bson_iterator_next(&it))
                                   ((std::complex<double>)this->v[index_coeff++]).imag(bson_iterator_double(&it));
                               }

                               bson_destroy(&br);
                               free_with_check(datar);
      }
        break;
#endif
      }

    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector() : Vector<Scalar>(), v(nullptr)
    {
    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector(unsigned int size) : Vector<Scalar>(size), v(nullptr)
    {
      if (this->size == 0)
        throw Exceptions::ValueException("size", this->size, 1);
      this->alloc(this->size);
    }

    template<typename Scalar>
    SimpleVector<Scalar>::~SimpleVector()
    {
      free();
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      SimpleVector<Scalar>* simple_vec = (SimpleVector<Scalar>*)vec;
      if (simple_vec)
        memcpy(this->v, simple_vec->v, sizeof(Scalar)*this->size);
      else
        Vector<Scalar>::set_vector(vec);
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::set_vector(Scalar* vec)
    {
      memcpy(this->v, vec, sizeof(Scalar)*this->size);
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      this->v = malloc_with_check<SimpleVector<Scalar>, Scalar>(n, this);
      zero();
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++)
        v[i] *= -1.;
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::zero()
    {
      memset(this->v, 0, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::free()
    {
      free_with_check(this->v);
      this->size = 0;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      this->v[idx] = y;
    }

    template<>
    void SimpleVector<double>::add(unsigned int idx, double y)
    {
#pragma omp atomic
      this->v[idx] += y;
    }

    template<>
    void SimpleVector<std::complex<double> >::add(unsigned int idx, std::complex<double> y)
    {
#pragma omp critical (SimpleVector_add)
      this->v[idx] += y;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        this->v[idx[i]] += y[i];
    }

    template<typename Scalar>
    Scalar SimpleVector<Scalar>::get(unsigned int idx) const
    {
      return this->v[idx];
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::extract(Scalar *v) const
    {
      memcpy(v, this->v, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i));
      return this;

    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector_multiple(Vector<Scalar>* vec, Scalar mult)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i) * mult);
      return this;

    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i]);
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector_multiple(Scalar* vec, Scalar mult)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i] * mult);
      return this;
    }

    template<>
    HERMES_API Vector<double>* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
      {
                                    return new SimpleVector<double>;
      }
      case Hermes::SOLVER_AMESOS:
      {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
                                  return new EpetraVector<double>;
#else
                                  throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
                                  break;
      }
      case Hermes::SOLVER_AZTECOO:
      {
                                   if (use_direct_solver)
                                     throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
                                   return new EpetraVector<double>;
#else
                                   throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_MUMPS:
      {
#ifdef WITH_MUMPS
                                 return new SimpleVector<double>;
#else
                                 throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
                                 break;
      }
      case Hermes::SOLVER_PETSC:
      {
                                 if (use_direct_solver)
                                   throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
                                 return new PetscVector<double>;
#else
                                 throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
                                 break;
      }
      case Hermes::SOLVER_UMFPACK:
      {
#ifdef WITH_UMFPACK
                                   return new SimpleVector<double>;
#else
                                   throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
      {
                                          if (use_direct_solver)
                                            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
                                          return new ParalutionVector<double>;
#else
                                          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
                                          break;
      }
      case Hermes::SOLVER_SUPERLU:
      {
#ifdef WITH_SUPERLU
                                   return new SimpleVector<double>;
#else
                                   throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
                                   break;
      }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return nullptr;
    }

    template<>
    HERMES_API Vector<std::complex<double> >* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
      {
                                    return new SimpleVector<std::complex<double> >;
      }

      case Hermes::SOLVER_AMESOS:
      {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
                                  return new EpetraVector<std::complex<double> >;
#else
                                  throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
                                  break;
      }
      case Hermes::SOLVER_AZTECOO:
      {
                                   if (use_direct_solver)
                                     throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");

#if defined HAVE_AZTECOO && defined HAVE_EPETRA
                                   return new EpetraVector<std::complex<double> >;
#else
                                   throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_MUMPS:
      {
#ifdef WITH_MUMPS
                                 return new SimpleVector<std::complex<double> >;
#else
                                 throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
                                 break;
      }
      case Hermes::SOLVER_PETSC:
      {
                                 if (use_direct_solver)
                                   throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
                                 return new PetscVector<std::complex<double> >;
#else
                                 throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
                                 break;
      }
      case Hermes::SOLVER_UMFPACK:
      {
#ifdef WITH_UMFPACK
                                   return new SimpleVector<std::complex<double> >;
#else
                                   throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
      {
                                          if (use_direct_solver)
                                            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
                                          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
                                          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
                                          break;
      }
      case Hermes::SOLVER_SUPERLU:
      {
#ifdef WITH_SUPERLU
                                   return new SimpleVector<std::complex<double> >;
#else
                                   throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
                                   break;
      }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return nullptr;
    }

    template class Vector<double>;
    template class Vector<std::complex<double> >;

    template class SimpleVector<double>;
    template class SimpleVector<std::complex<double> >;
  }
}
