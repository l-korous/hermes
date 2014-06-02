#include "algorithms.h"

// Relative tolerance drop (1e-4 == 4 orders of magnitude drop)
static const double tolerance = 1.e-3;

static const int integrationOrder = 4;

// Utilities
static std::string SolvedExampleString[5] = { "1D", "CircularConvection", "MovingPeak", "AdvectedCube", "SolidBodyRotation" };
static double exact_solver_error;
double initial_error;
MeshFunctionSharedPtr<double> es(new Solution<double>());
double* es_v;

// Uncomment to have OpenGL output throughout calculation.
#define SHOW_OUTPUT

// Under relaxation in Multiscale
#define OMEGA 1.0

// Static logging for output in terminal.
static Hermes::Mixins::Loggable static_log(true);

double calc_l2_error_algebraic(SpaceSharedPtr<double> space, double* v1, double* v2,
                               Hermes::Mixins::Loggable* logger = NULL, int iteration = 0, int init_refs = 0, double D = 0.)
{
  double result = 0.;
  for (int i = 0; i < space->get_num_dofs(); i++)
    result += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  result = std::sqrt(result);
  if (logger)
    logger->info("%d,%f,%d,%f", init_refs, D, iteration, result);
  return result;
}

static
double calc_l2_norm_algebraic(SpaceSharedPtr<double> space, double* v1)
{
  double result = 0.;
  for (int i = 0; i < space->get_num_dofs(); i++)
    result += v1[i]*v1[i];
  result = std::sqrt(result);
  return result;
}

bool error_reduction_condition(double error)
{
  return std::abs(error / initial_error) < tolerance;
}

void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, 
                 MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step,
                 int poly_degree, int init_ref_num)
{
  MeshFunctionSharedPtr<double> exact_solver_sln(new Solution<double>());
  ScalarView* exact_solver_view = new ScalarView("Exact solver solution", new WinGeom(0, 360, 600, 350));

  // Standard L2 space.
  int full_ndofs = space->get_num_dofs();

  // Exact solver
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, initial_sln);
  weakform_exact.set_current_time_step(time_step);
  CSCMatrix<double> matrix_A;
  SimpleVector<double> vector_b;
//   LinearSolver<double> solver_exact(&weakform_exact, space);

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A);

  UMFPackLinearMatrixSolver<double> solver_exact(&matrix_A, &vector_b);
  solver_exact.setup_factorization();
  solver_exact.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);


  // Solve
  dp.assemble(&vector_b);
  vector_b.export_to_file("exact", "b", MatrixExportFormat::EXPORT_FORMAT_PLAIN_ASCII);
  solver_exact.solve();
//   solver_exact.solve();
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, exact_solver_sln);

  // Initial error
  initial_error = get_l2_norm(solver_exact.get_sln_vector(), space->get_num_dofs());

  // es and es_v
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, es);
  es_v = new double[space->get_num_dofs()];
  memcpy(es_v, solver_exact.get_sln_vector(), sizeof(double)* space->get_num_dofs());

  // output
  std::stringstream ss_bmp, ss_vtk;
  ss_bmp.precision(2);
  ss_vtk.precision(2);
  ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);

  ss_bmp << "solution_Exact_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".bmp";
  ss_vtk << "solution_Exact_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".dat";
#ifdef SHOW_OUTPUT
  exact_solver_view->show(es);
  exact_solver_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
  exact_solver_view->get_linearizer()->save_solution_tecplot(es, ss_vtk.str().c_str(), "solution");
}

void exact_solver_timedep(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree,
                          int init_ref_num, double diffusivity, double s, double sigma, double time_step_length,
                          int time_step_count, MeshFunctionSharedPtr<double> previous_solution, 
                          MeshFunctionSharedPtr<double> exact_solution, ScalarView* exact_view, double cfl)
{
  if (!is_timedep(solvedExample))
    return;

  // Standard L2 space.
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int full_ndofs = full_space->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakFormTimedep weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  weakform_exact.set_ext(previous_solution);
  CSCMatrix<double> matrix_A;
  SimpleVector<double> vector_b;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A);

  UMFPackLinearMatrixSolver<double> solver(&matrix_A, &vector_b);
  solver.setup_factorization();
  solver.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Reporting.

  double time = 0.;
  for (int iteration = 1; iteration <= time_step_count; ++iteration)
  {
    static_log.info("Time step: %i, time: %f.", iteration, time+time_step_length);
    dp.assemble(&vector_b);
    solver.solve();
    Solution<double>::vector_to_solution(solver.get_sln_vector(), full_space, previous_solution);
    time += time_step_length;

#ifdef SHOW_OUTPUT
    exact_view->show(previous_solution);
#endif
  }
  
  Solution<double>::vector_to_solution(solver.get_sln_vector(), full_space, es);

  std::stringstream ss_bmpe;
  std::stringstream ss_vtke;
  ss_vtke.precision(2);
  ss_vtke.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_vtke << "solution_Exact_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" 
          << diffusivity << "_CFL=" << cfl << ".dat";

  ss_bmpe.precision(2);
  ss_bmpe.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_bmpe << "solution_Exact_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" 
          << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
  exact_view->show(es);
  exact_view->save_screenshot(ss_bmpe.str().c_str(), true);
  exact_view->close();
#endif
  exact_view->get_linearizer()->save_solution_tecplot(es, ss_vtke.str().c_str(), "exactSolution", 1, 2.0);
}


std::string multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, 
                                     MeshFunctionSharedPtr<double> previous_mean_values, 
                                     MeshFunctionSharedPtr<double> previous_derivatives, 
                                     double diffusivity, double s, double sigma, double time_step_length, int time_step_count, 
                                     MeshFunctionSharedPtr<double> exact_solution, ScalarView* solution_view, 
                                     ScalarView* exact_view, Hermes::Mixins::Loggable& logger, 
                                     Hermes::Mixins::Loggable& logger_details, double cfl)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));

  int ndofs = space->get_num_dofs();
  int const_ndofs = const_space->get_num_dofs();
  int full_ndofs = full_space->get_num_dofs();

  // Utility solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>());

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", 
                               diffusivity, s, sigma, exact_solution);
  MultiscaleWeakForm weakform_implicit(solvedExample, add_inlet(solvedExample), "Inlet", 
                                       diffusivity, s, sigma, exact_solution, false);
  MultiscaleWeakForm weakform_explicit(solvedExample, add_inlet(solvedExample), "Inlet", 
                                       diffusivity, s, sigma, exact_solution, true);
  ExplicitWeakFormOffDiag weakform_explicit_offdiag(solvedExample, add_inlet(solvedExample), "Inlet", 
                                                    diffusivity, s, sigma);
  MassWeakForm weakform_mass;
  weakform_exact.set_current_time_step(time_step_length);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_explicit.set_current_time_step(time_step_length);
  weakform_explicit_offdiag.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_full;
  
  CSCMatrix<double> matrix_A_der;
  CSCMatrix<double> matrix_A_der_offdiag;
  CSCMatrix<double> matrix_M_der;
  SimpleVector<double> vector_b_der;

  CSCMatrix<double> matrix_A_means;
  CSCMatrix<double> matrix_M_means;
  SimpleVector<double> vector_b_means;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_full);

  // Level 1.
  dp.set_space(space);
  dp.set_weak_formulation(&weakform_explicit);
  dp.assemble(&matrix_A_der, &vector_b_der);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_der);
  dp.set_weak_formulation(&weakform_explicit_offdiag);
  dp.assemble(&matrix_A_der_offdiag);

  // Level 0.
  dp.set_space(const_space);
  dp.set_weak_formulation(&weakform_implicit);
  dp.assemble(&matrix_A_means, &vector_b_means);
  vector_b_means.export_to_file("hss", "b", MatrixExportFormat::EXPORT_FORMAT_PLAIN_ASCII);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_means);

  SimpleVector<double> vector_A_der(ndofs);
  SimpleVector<double> vector_A_means(const_ndofs);

  UMFPackLinearMatrixSolver<double> solver_means(&matrix_A_means, &vector_A_means);
  solver_means.setup_factorization();
  solver_means.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_der(&matrix_A_der, &vector_A_der);
  solver_der.setup_factorization();
  solver_der.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  SimpleVector<double> sln_means_k(const_ndofs);
  SimpleVector<double> sln_means_k_long(full_ndofs);
  SimpleVector<double> A_times_sln_means_k_long(full_ndofs);
  SimpleVector<double> sln_der_k(ndofs);
  SimpleVector<double> sln_der_k_long(full_ndofs);
  SimpleVector<double> A_times_sln_der_k_long(full_ndofs);
  SimpleVector<double> sln_der_offdiag(ndofs);
  // For residual calculation
  SimpleVector<double> A_times_u(full_ndofs);

  OGProjection<double>::project_global(const_space, previous_mean_values, sln_means_k.v);
  OGProjection<double>::project_global(space, previous_derivatives, sln_der_k.v);

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;
  double time = 0.;

  double* merged_sln = new double[full_ndofs];

  SimpleVector<double> temp_1(const_ndofs);
  SimpleVector<double> temp_2(ndofs);

  // In the following:
  // u0 - means
  // u0_prev - means from previous iteration
  // u* - derivatives
  // u*_prev - derivatives from previous iteration
  // sln_der_k_long - prolongated (added 0s for means DOFs) derivatives from previous iteration
  // A_times_sln_der_k_long - prolongated (added 0s for means DOFs) derivatives from previous iteration
  // sln_means_k_long - prolongated (added 0s for derivatives DOFs) means from previous iteration
  // A_times_sln_means_k_long - prolongated (added 0s for derivatives DOFs) means from previous iteration
  // M0 - mass matrix for means
  // M* - mass matrix for derivatives
  // M - mass matrix for the complete system
  // A - the full (implicit) A matrix for the complete system
  // A0 - the full (implicit) A matrix for the means
  // A* - the full (implicit) A matrix for the derivatives (not really used)
  // A*~ - only the diagonal (explicit) part of the A matrix for the derivatives
  // b0 - RHS for means
  // b* - RHS for derivatives
  do
  {
    iterations++;
    // 1. means
    // M0 * u0_prev -> vector_A_means (utility variable)
    matrix_M_means.multiply_with_vector(sln_means_k.v, vector_A_means.v, true);
    if (polynomialDegree)
    {
      // Prolongate u*_prev -> sln_der_k_long
      add_means(&sln_der_k, &sln_der_k_long, space, full_space);
      // A * sln_der_k_long -> A_times_sln_der_k_long
      matrix_A_full.multiply_with_vector(sln_der_k_long.v, A_times_sln_der_k_long.v, true);
      // Restrict (take out derivatives DOFs) A_times_sln_der_k_long -> temp_1.
      cut_off_ders(A_times_sln_der_k_long.v, const_space, full_space, temp_1.v);
      // Basically after the following, this will hold: vector_A_means = [M0 * u0_prev - A * u*_prev] (but the multiplication with A had to be done with the full A)
      // vector_A_means -= temp_1
      vector_A_means.add_vector(temp_1.change_sign());
    }
    // vector_A_means += b0
    vector_A_means.add_vector(&vector_b_means);
    // Solve A0 * u0 = vector_A_means for u0
    solver_means.solve();
    // Store the solution into the previous iteration:
    // u0 -> u0_prev
    sln_means_k.set_vector(solver_means.get_sln_vector());

    if (polynomialDegree)
    {
      // 2. corrector
      // M* * u*_prev -> vector_A_der (utility variable)
      matrix_M_der.multiply_with_vector(sln_der_k.v, vector_A_der.v, true);
      // Prolongate u0_prev -> sln_means_k_long
      add_ders(&sln_means_k, &sln_means_k_long, const_space, full_space);
      // A * sln_means_k_long -> A_times_sln_means_k_long
      matrix_A_full.multiply_with_vector(sln_means_k_long.v, A_times_sln_means_k_long.v, true);
      // Restrict (take out means DOFs) A_times_sln_means_k_long -> temp_2.
      cut_off_means(A_times_sln_means_k_long.v, space, full_space, temp_2.v);
      // Basically after the following, this will hold: vector_A_der = [M* * u*_prev + A * u0_prev] (but the multiplication with A had to be done with the full A)
      // vector_A_der -= temp_2
      vector_A_der.add_vector(temp_2.change_sign());
      // (A* - A*~) * u*_prev -> sln_der_offdiag (utility variable)
      matrix_A_der_offdiag.multiply_with_vector(sln_der_k.v, sln_der_offdiag.v, true);
      // vector_A_der -= sln_der_offdiag
      vector_A_der.add_vector(sln_der_offdiag.change_sign());
      // vector_A_der += b*
      vector_A_der.add_vector(&vector_b_der);
      // Solve A*~ * u* = vector_A_der for u*
      solver_der.solve();
      // Store the solution into the previous iteration:
      // u* -> u*_prev   (potentially use relaxation)
      if (OMEGA >= 0.99)
        sln_der_k.set_vector(solver_der.get_sln_vector());
      else
      {
        for (int i = 0; i < ndofs; i++)
          sln_der_k.set(i, (OMEGA * solver_der.get_sln_vector()[i]) + ((1. - OMEGA) * sln_der_k.get(i)));
      }

      // u0 + u* -> u
      merge_slns(sln_means_k.v, const_space, sln_der_k.v, space, full_space, false, merged_sln);
    }

#ifdef SHOW_OUTPUT
    if (polynomialDegree)
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    else
      Solution<double>::vector_to_solution(sln_means_k.v, const_space, solution);

    solution_view->set_title("Time: %f.");
    solution_view->show(solution);
#endif
    // Computation of the current residual
    // !!!!! THIS DOES NOT DO ANYTHING, JUST SOME CALCULATION AND THE OUTPUT GOES NOWHERE
    {
      // A * u -> A_times_u
      matrix_A_full.multiply_with_vector(polynomialDegree ? merged_sln : sln_means_k.v, A_times_u.v, true);
      // Change sign of A_times_u
      A_times_u.change_sign();
      
      // A_times_u += b
      if(polynomialDegree)
        add_means(&vector_b_der, &sln_der_k_long, space, full_space);
      if (polynomialDegree)
        A_times_u.add_vector(&sln_der_k_long);
      
      add_ders(&vector_b_means, &sln_means_k_long, const_space, full_space);
      A_times_u.add_vector(&sln_means_k_long);

      // Algebraic norm of (b - A * u)
      double current_residual = calc_l2_norm_algebraic(polynomialDegree ? full_space : const_space, A_times_u.v);
      std::cout << current_residual << std::endl;
    }

    bool done=error_reduction_condition(calc_l2_error_algebraic(polynomialDegree ? full_space : const_space, 
      polynomialDegree ? merged_sln : sln_means_k.v, es_v, &logger_details, iterations, init_ref_num, diffusivity));
    if (done)
      break;
  } while(true);

  std::stringstream outStream;
  outStream << "Iter=" << iterations;

  return outStream.str();
}

std::string multiscale_decomposition_timedep(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, 
                                             int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values, 
                                             MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, 
                                             double s, double sigma, double time_step_length, int time_step_count,
                                             MeshFunctionSharedPtr<double> exact_solution, ScalarView* solution_view, 
                                             ScalarView* exact_view, Hermes::Mixins::Loggable& logger, 
                                             Hermes::Mixins::Loggable& logger_details, double cfl)
{
  if (!is_timedep(solvedExample))
    return "";

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));

  int ndofs = space->get_num_dofs();
  int const_ndofs = const_space->get_num_dofs();
  int full_ndofs = full_space->get_num_dofs();

  // Utility solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>());

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", 
                               diffusivity, s, sigma, exact_solution);
  MultiscaleWeakForm weakform_implicit(solvedExample, add_inlet(solvedExample), "Inlet", 
                                       diffusivity, s, sigma, exact_solution, false);
  MultiscaleWeakForm weakform_explicit(solvedExample, add_inlet(solvedExample), "Inlet", 
                                       diffusivity, s, sigma, exact_solution, true);
  ExplicitWeakFormOffDiag weakform_explicit_offdiag(solvedExample, add_inlet(solvedExample), "Inlet", 
                                                    diffusivity, s, sigma);
  MassWeakForm weakform_mass;
  weakform_exact.set_current_time_step(time_step_length);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_explicit.set_current_time_step(time_step_length);
  weakform_explicit_offdiag.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_full;
  CSCMatrix<double> matrix_A_means_just_A;
  CSCMatrix<double> matrix_A_der;
  SimpleVector<double> vector_b_der;
  CSCMatrix<double> matrix_M_der;
  CSCMatrix<double> matrix_A_der_offdiag;

  CSCMatrix<double> matrix_A_means;
  SimpleVector<double> vector_b_means;
  CSCMatrix<double> matrix_M_means;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_full);
  dp.set_space(const_space);
  dp.assemble(&matrix_A_means_just_A);

  // Level 1.
  dp.set_space(space);
  dp.set_weak_formulation(&weakform_explicit);
  dp.assemble(&matrix_A_der, &vector_b_der);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_der);
  dp.set_weak_formulation(&weakform_explicit_offdiag);
  dp.assemble(&matrix_A_der_offdiag);

  // Level 0.
  dp.set_space(const_space);
  dp.set_weak_formulation(&weakform_implicit);
  dp.assemble(&matrix_A_means, &vector_b_means);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_means);

  SimpleVector<double> vector_A_der(ndofs);
  SimpleVector<double> vector_A_means(const_ndofs);

  UMFPackLinearMatrixSolver<double> solver_means(&matrix_A_means, &vector_A_means);
  solver_means.setup_factorization();
  solver_means.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_der(&matrix_A_der, &vector_A_der);
  solver_der.setup_factorization();
  solver_der.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  SimpleVector<double> sln_means(const_ndofs);
  SimpleVector<double> sln_means_k(const_ndofs);
  SimpleVector<double> sln_means_tmp(const_ndofs);
  SimpleVector<double> sln_means_long(full_ndofs);
  SimpleVector<double> sln_means_long_temp(full_ndofs);
  SimpleVector<double> sln_der(ndofs);
  SimpleVector<double> sln_der_k(ndofs);
  SimpleVector<double> sln_der_tmp(ndofs);
  SimpleVector<double> sln_der_long(full_ndofs);
  SimpleVector<double> A_times_sln_der_long(full_ndofs);
  SimpleVector<double> sln_der_offdiag(ndofs);

  OGProjection<double>::project_global(const_space, previous_mean_values, sln_means.v);
  OGProjection<double>::project_global(space, previous_derivatives, sln_der.v);

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;
  double time = 0.;

  double* merged_sln = new double[full_ndofs];

  SimpleVector<double> temp_1(const_ndofs);
  SimpleVector<double> temp_2(ndofs);
  
  for (int time_step = 1; time_step <= time_step_count; time_step++)
  {
    double initial_residual, current_residual;

    static_log.info("Time step: %i, time: %f.", time_step, time+time_step_length);

    // Computation of the initial residual
    add_means(&sln_der, &sln_der_long, space, full_space);
    add_ders(&sln_means, &sln_means_long, const_space, full_space);
    sln_der_long.add_vector(&sln_means_long);
    matrix_A_full.multiply_with_vector(sln_der_long.v, A_times_sln_der_long.v, true);
    A_times_sln_der_long.change_sign();
    add_means(&vector_b_der, &sln_der_long, space, full_space);
    add_ders(&vector_b_means, &sln_means_long, const_space, full_space);
    A_times_sln_der_long.add_vector(&sln_means_long);
    A_times_sln_der_long.add_vector(&sln_der_long);
    initial_residual=calc_l2_norm_algebraic(full_space, A_times_sln_der_long.v);

    do
    {
      iterations++;
      // 2. means
      // M
      matrix_M_means.multiply_with_vector(sln_means.v, vector_A_means.v, true);
      // -B
      add_means(&sln_der_k, &sln_der_long, space, full_space);
      matrix_A_full.multiply_with_vector(sln_der_long.v, A_times_sln_der_long.v, true);
      cut_off_ders(A_times_sln_der_long.v, const_space, full_space, temp_1.v);
      temp_1.change_sign();
      vector_A_means.add_vector(&temp_1);

      // b
      vector_A_means.add_vector(&vector_b_means);
      // SOLVE
      solver_means.solve();
      sln_means_k.set_vector(solver_means.get_sln_vector());

      // 3. corrector
      // M
      matrix_M_der.multiply_with_vector(sln_der.v, vector_A_der.v, true);
      // -B
      add_ders(&sln_means_k, &sln_means_long, const_space, full_space);
      matrix_A_full.multiply_with_vector(sln_means_long.v, sln_means_long_temp.v, true);
      cut_off_means(sln_means_long_temp.v, space, full_space, temp_2.v);
      temp_2.change_sign();
      vector_A_der.add_vector(&temp_2);
      // (A - A~) - offdiag
      matrix_A_der_offdiag.multiply_with_vector(sln_der_k.v, sln_der_offdiag.v, true);
      vector_A_der.add_vector(sln_der_offdiag.change_sign());
      // b
      vector_A_der.add_vector(&vector_b_der);
      // SOLVE
      solver_der.solve();

      // Computation of the current residual
      add_means(&sln_der_k, &sln_der_long, space, full_space);
      add_ders(&sln_means_k, &sln_means_long, const_space, full_space);
      sln_der_long.add_vector(&sln_means_long);
      matrix_A_full.multiply_with_vector(sln_der_long.v, A_times_sln_der_long.v, true);
      A_times_sln_der_long.change_sign();
      sln_means_tmp.set_vector(&sln_means_k);
      sln_means_tmp.change_sign()->add_vector(&sln_means);
      matrix_M_means.multiply_with_vector(sln_means_tmp.v, vector_A_means.v, true);
      vector_A_means.add_vector(&vector_b_means);
      sln_der_tmp.set_vector(&sln_der_k);
      sln_der_tmp.change_sign()->add_vector(&sln_der);
      matrix_M_der.multiply_with_vector(sln_der_tmp.v, vector_A_der.v, true);
      vector_A_der.add_vector(&vector_b_der);
      add_ders(&vector_A_means, &sln_means_long, const_space, full_space);
      add_means(&vector_A_der, &sln_der_long, space, full_space);
      A_times_sln_der_long.add_vector(&sln_means_long);
      A_times_sln_der_long.add_vector(&sln_der_long);
      current_residual=calc_l2_norm_algebraic(full_space, A_times_sln_der_long.v);

      if (OMEGA >= 0.99)
        sln_der_k.set_vector(solver_der.get_sln_vector());
      else
      {
        for (int i = 0; i < ndofs; i++)
          sln_der_k.set(i, (OMEGA * solver_der.get_sln_vector()[i]) + ((1. - OMEGA) * sln_der_k.get(i)));
      }
    } while (current_residual/initial_residual>tolerance);

    sln_means.set_vector(&sln_means_k);
    sln_der.set_vector(&sln_der_k);
#ifdef SHOW_OUTPUT
    if (polynomialDegree)
    {
      merge_slns(sln_means.v, const_space, sln_der.v, space, full_space, false, merged_sln);
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    }
    else
      Solution<double>::vector_to_solution(sln_means.v, const_space, solution);

    solution_view->set_title("Time: %f.", time + time_step_length);
    solution_view->show(solution);
#endif

    time += time_step_length;
  }

  std::stringstream outStream;
  outStream << "Iter=" << iterations;
  DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

  std::stringstream ss_vtk;
  std::stringstream ss_bmp;
  ss_vtk.precision(2);
  ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_vtk << "solution_HSS_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";

//   ss_bmp.precision(2);
//   ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
//   ss_bmp << "solution_HSS(" << steps_per_time_step << ")_" << SolvedExampleString[solvedExample] 
//          << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
  //solution_view->show(solution);
  //solution_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
  //solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);

  errorCalculator.calculate_errors(solution, es);

  outStream << "|" << "Err=" << std::sqrt(errorCalculator.get_total_error_squared());

  return outStream.str();
}

std::string p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num,
                        MeshFunctionSharedPtr<double> previous_sln, double diffusivity, double time_step_length,
                        int time_step_count, MeshFunctionSharedPtr<double> exact_solution, ScalarView* solution_view, ScalarView* exact_view,
                        double s, double sigma, Hermes::Mixins::Loggable& logger, int smoothing_steps_per_V_cycle, 
                        double cfl)
{
  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Utility solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>());

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_2;
  SimpleVector<double> vector_b_2;
  CSCMatrix<double> matrix_A_1;
  SimpleVector<double> vector_b_1;
  CSCMatrix<double> matrix_A_0;
  SimpleVector<double> vector_b_0;

  // Matrices (M+A_tilde), vectors -A(u_K)
  SmoothingWeakForm weakform_smoother(solvedExample, true, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  SmoothingWeakForm weakform_smoother_coarse(solvedExample, false, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  weakform_smoother_coarse.set_current_time_step(time_step_length);
  weakform_smoother_coarse.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  MassWeakForm weakform_mass;

  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2(ndofs_2);
  CSCMatrix<double> matrix_MA_tilde_1;
  SimpleVector<double> vector_A_1(ndofs_1);
  CSCMatrix<double> matrix_MA_0;
  SimpleVector<double> vector_A_0(ndofs_0);

  CSCMatrix<double> matrix_M_2;
  CSCMatrix<double> matrix_M_1;
  CSCMatrix<double> matrix_M_0;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(space_2);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_2, &vector_b_2);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_2);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_2);
  if (is_timedep(solvedExample))
    matrix_MA_tilde_2.add_sparse_matrix(&matrix_M_2);

  // Level 1.
  dp.set_space(space_1);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_1, &vector_b_1);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_1);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_1);
  if (is_timedep(solvedExample))
    matrix_MA_tilde_1.add_sparse_matrix(&matrix_M_1);

  // Level 0.
  dp.set_space(space_0);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_0, &vector_b_0);
  dp.set_weak_formulation(&weakform_smoother_coarse);
  dp.assemble(&matrix_MA_0);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_0);
  if (is_timedep(solvedExample))
    matrix_MA_0.add_sparse_matrix(&matrix_M_0);

  UMFPackLinearMatrixSolver<double> solver_2(&matrix_MA_tilde_2, &vector_A_2);
  solver_2.setup_factorization();
  solver_2.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_1(&matrix_MA_tilde_1, &vector_A_1);
  solver_1.setup_factorization();
  solver_1.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_0(&matrix_MA_0, &vector_A_0);
  solver_0.setup_factorization();
  solver_0.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  SimpleVector<double> sln_2(ndofs_2);
  SimpleVector<double> sln_1(ndofs_1);
  SimpleVector<double> sln_0(ndofs_0);

  SimpleVector<double> prev_sln_2(ndofs_2);
  SimpleVector<double> prev_sln_1(ndofs_1);
  SimpleVector<double> prev_sln_0(ndofs_0);

  SimpleVector<double> util_2(ndofs_2), util_21(ndofs_2);
  SimpleVector<double> util_1(ndofs_1), util_11(ndofs_2);
  SimpleVector<double> util_0(ndofs_0), util_01(ndofs_2);

  SimpleVector<double> projected_A_P_1(ndofs_1), sln_2_projected(ndofs_1);
  SimpleVector<double> projected_A_P_0(ndofs_0), sln_1_projected(ndofs_0);
  SimpleVector<double> temp(ndofs_0);

  // Reports.
  int num_coarse = 0;
  int num_2 = 0;
  int num_1 = 0;
  int v_cycles = 0;

  OGProjection<double>::project_global(space_2, previous_sln, &sln_2);
  prev_sln_2.set_vector(&sln_2);
  cut_off_quadratic_part(prev_sln_2.v, space_1, space_2, prev_sln_1.v);
  cut_off_linear_part(prev_sln_1.v, space_0, space_1, prev_sln_0.v);

  double time = 0.;
  if (!is_timedep(solvedExample))
    time_step_count = 1;
  for (int time_step = 1; time_step <= time_step_count; time_step++)
  {
    double initial_residual, current_residual;

    if (is_timedep(solvedExample))
      static_log.info("Time step: %i, time: %f.", time_step, time+time_step_length);
    else
      static_log.info("Time step: %i.", time_step);

    // Computation of the initial residual
    if (is_timedep(solvedExample))
    {
      matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
      vector_A_2.change_sign()->add_vector(&vector_b_2);
      initial_residual=calc_l2_norm_algebraic(space_2, vector_A_2.v);
    }

    do
    {
      v_cycles++;

#pragma region 0 - highest level
      // Store the previous solution.
      if (polynomialDegree > 1)
      {
        for (int smoothing_step = 1; smoothing_step <= smoothing_steps_per_V_cycle; smoothing_step++)
        {
          // Solve for increment.
          matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
          vector_A_2.change_sign()->add_vector(&vector_b_2);

          if (is_timedep(solvedExample))
          {
            util_2.set_vector(&prev_sln_2);
            util_2.change_sign()->add_vector(sln_2.v)->change_sign();
            matrix_M_2.multiply_with_vector(util_2.v, util_21.v, true);

            vector_A_2.add_vector(util_21.v);
          }

          solver_2.solve();
          sln_2.add_vector(solver_2.get_sln_vector());
        }
      }

#pragma endregion

#pragma region 1 - intermediate level

      // f_P1
      SimpleVector<double> f_P1(ndofs_1);
      f_P1.zero();
      // Minus A_P1
      SimpleVector<double> R_P1(ndofs_1);
      // Minus(minus) projected_A_P1
      SimpleVector<double> projected_A_2(ndofs_2);
      matrix_A_2.multiply_with_vector(sln_2.v, projected_A_2.v, true);

      cut_off_quadratic_part(projected_A_2.v, space_1, space_2, projected_A_P_1.v);

      cut_off_quadratic_part(sln_2.v, space_1, space_2, sln_2_projected.v);
      matrix_A_1.multiply_with_vector(sln_2_projected.v, R_P1.v, true);
      
      sln_1.set_vector(&sln_2_projected);

      R_P1.change_sign();
      f_P1.add_vector(&R_P1);
      f_P1.add_vector(&projected_A_P_1);
      f_P1.change_sign();

      for (int smoothing_step = 1; smoothing_step <= smoothing_steps_per_V_cycle; smoothing_step++)
      {
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if (polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);

        if (is_timedep(solvedExample))
        {
          util_1.set_vector(&prev_sln_1);
          util_1.change_sign()->add_vector(sln_1.v)->change_sign();
          matrix_M_1.multiply_with_vector(util_1.v, util_11.v, true);

          vector_A_1.add_vector(util_11.v);
        }

        solver_1.solve();
        sln_1.add_vector(solver_1.get_sln_vector());
      }

#pragma endregion

#pragma region  2 - Solve the problem on the coarse level exactly

      // f_P0
      SimpleVector<double> f_P0(ndofs_0);
      f_P0.zero();
      // Minus A_P0
      SimpleVector<double> R_P0(ndofs_0);
      // Minus(minus) projected_A_P0
      SimpleVector<double> projected_A_1(ndofs_1);
      matrix_A_1.multiply_with_vector(sln_1.v, projected_A_1.v, true);

      cut_off_linear_part(projected_A_1.v, space_0, space_1, projected_A_P_0.v);

      cut_off_linear_part(sln_1.v, space_0, space_1, sln_1_projected.v);
      matrix_A_0.multiply_with_vector(sln_1_projected.v, R_P0.v, true);

      SimpleVector<double> projected_f_P1(ndofs_1);
      projected_f_P1.set_vector(&f_P1);

      sln_0.set_vector(&sln_1_projected);

      R_P0.change_sign();
      f_P0.add_vector(&R_P0);
      f_P0.add_vector(&projected_A_P_0);
      if (polynomialDegree > 1)
      {
        cut_off_linear_part(projected_f_P1.v, space_0, space_1, temp.v);
        temp.change_sign();
        f_P0.add_vector(&temp);
      }
      f_P0.change_sign();

      num_coarse++;

      // A(u_K) - done after the first step.
      if (is_timedep(solvedExample))
      {
        matrix_M_0.multiply_with_vector(prev_sln_0.v, vector_A_0.v, true);
        vector_A_0.add_vector(&f_P0)->add_vector(&vector_b_0);
      }
      else
      {
        matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);
        vector_A_0.change_sign()->add_vector(&f_P0)->add_vector(&vector_b_0);
      }
      solver_0.solve();
      if (is_timedep(solvedExample))
        sln_0.set_vector(solver_0.get_sln_vector());
      else
        sln_0.add_vector(solver_0.get_sln_vector());

#pragma endregion

#pragma region 1 - intermediate level
      // Store the previous solution.
      merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, false, sln_1.v);

      for (int smoothing_step = 1; smoothing_step <= smoothing_steps_per_V_cycle; smoothing_step++)
      {
        // Solve for increment.
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if (polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);

        if (is_timedep(solvedExample))
        {
          util_1.set_vector(&prev_sln_1);
          util_1.change_sign()->add_vector(sln_1.v)->change_sign();
          matrix_M_1.multiply_with_vector(util_1.v, util_11.v, true);

          vector_A_1.add_vector(util_11.v);
        }

        solver_1.solve();
        sln_1.add_vector(solver_1.get_sln_vector());
      }

#pragma endregion

#pragma region 0 - highest level

      if (polynomialDegree > 1)
      {
        // Store the previous solution.
        merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false, sln_2.v);

        for (int smoothing_step = 1; smoothing_step <= smoothing_steps_per_V_cycle; smoothing_step++)
        {
          // Solve for increment.
          matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
          vector_A_2.change_sign()->add_vector(&vector_b_2);

          if (is_timedep(solvedExample))
          {
            util_2.set_vector(&prev_sln_2);
            util_2.change_sign()->add_vector(sln_2.v)->change_sign();
            matrix_M_2.multiply_with_vector(util_2.v, util_21.v, true);
            vector_A_2.add_vector(util_21.v);
          }

          solver_2.solve();
          sln_2.add_vector(solver_2.get_sln_vector());
        }
      }
#pragma endregion

#ifdef SHOW_OUTPUT
      // Make solution & display.
      Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
      solution_view->set_title("Time: %f, V_cycle: %i", time, v_cycles);
      solution_view->show(previous_sln);
#endif

      // Computation of the current residual
      // Decision whether or not to end this loop.
      if (is_timedep(solvedExample))
      {
        matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
        vector_A_2.change_sign()->add_vector(&vector_b_2);
        util_2.set_vector(&sln_2);
        util_2.change_sign()->add_vector(prev_sln_2.v);
        matrix_M_2.multiply_with_vector(util_2.v, util_21.v, true);
        vector_A_2.add_vector(util_21.v);
        current_residual=calc_l2_norm_algebraic(space_2, vector_A_2.v);

        if (current_residual / initial_residual < tolerance)
          break;
      }
      else
      {
        if (error_reduction_condition(calc_l2_error_algebraic(space_2, sln_2.v, es_v)))
          break;
      }

    }
    while (true);

    if (is_timedep(solvedExample))
      time += time_step_length;

    prev_sln_2.set_vector(&sln_2);
    prev_sln_1.set_vector(&sln_1);
    prev_sln_0.set_vector(&sln_0);
  }

  std::stringstream outStream;
  outStream << "Iter=" << v_cycles;
  if (is_timedep(solvedExample))
  {
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

    Solution<double>::vector_to_solution(&sln_2, space_2, solution);

    std::stringstream ss_vtk;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_MG(" << smoothing_steps_per_V_cycle << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);

    std::stringstream ss_bmp;
    ss_bmp.precision(2);
    ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmp << "solution_MG(" << smoothing_steps_per_V_cycle << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
    solution_view->show(solution);
    solution_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
    errorCalculator.calculate_errors(solution, es);
    outStream << "|" << "Err=" << std::sqrt(errorCalculator.get_total_error_squared());
  }

  return outStream.str();
}

// Utilities.
bool add_inlet(SolvedExample solvedExample)
{
  switch (solvedExample)
  {
  case SolidBodyRotation:
  case AdvectedCube:
  case MovingPeak:
    return false;
  case CircularConvection:
  case Benchmark:
    return true;
  }
}

bool is_timedep(SolvedExample solvedExample)
{
  switch (solvedExample)
  {
  case CircularConvection:
  case Benchmark:
    return false;
  case MovingPeak:
  case AdvectedCube:
  case SolidBodyRotation:
    return true;
  }
}

// double end_time(SolvedExample solvedExample)
// {
//   switch (solvedExample)
//   {
//   case CircularConvection:
//   case Benchmark:
//     return 9999999999.;
//   case MovingPeak:
//     return M_PI * 2.;
//   case AdvectedCube:
//     return 1.;
//   case SolidBodyRotation:
//     return M_PI * 2.;
//   }
// }
