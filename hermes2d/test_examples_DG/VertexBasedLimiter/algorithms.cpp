#include "algorithms.h"

static double exact_solver_error;
static const double wrt_exact_solver_tolerance = 1e-6;

double calc_l2_error(SolvedExample solvedExample, MeshSharedPtr mesh, MeshFunctionSharedPtr<double> fn_1, MeshFunctionSharedPtr<double> fn_2, Hermes::Mixins::Loggable& logger)
{
  ErrorWeakForm wf(solvedExample);
  SpaceSharedPtr<double> mspace(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  wf.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(fn_1, fn_2));
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, mspace);
  SimpleVector<double> vector;
  dp->assemble(&vector);
  double result = 0.;
  for(int i = 0; i < vector.get_size(); i++)
    result += vector.get(i);
  result = std::sqrt(result);
  logger.info("L2 Error: %g.", result);
  return result;
}

void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step, Hermes::Mixins::Loggable& logger)
{
  MeshFunctionSharedPtr<double> exact_solver_sln(new Solution<double>());
  ScalarView* exact_solver_view = new ScalarView("Exact solver solution", new WinGeom(0, 360, 600, 350));

  // Exact solver
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", diffusivity, s, sigma, initial_sln);
  weakform_exact.set_current_time_step(time_step);
  LinearSolver<double> solver_exact(&weakform_exact, space);
  solver_exact.solve();
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, exact_solver_sln);
  exact_solver_error = calc_l2_error(solvedExample, space->get_mesh(), exact_solver_sln, exact_solution, logger);
  exact_solver_view->show(exact_solver_sln);
}

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_mean_values, 
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, 
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  ImplicitWeakForm weakform_implicit(solvedExample, true, "Inlet", diffusivity, s, sigma);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, exact_solution));
  ExplicitWeakForm weakform_explicit(solvedExample, true, "Inlet", diffusivity, s, sigma);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, exact_solution));
  weakform_explicit.set_current_time_step(time_step_length);
  LinearSolver<double> solver_implicit(&weakform_implicit, const_space);
  solver_implicit.set_jacobian_constant();
  LinearSolver<double> solver_explicit(&weakform_explicit, space);
  solver_explicit.set_jacobian_constant();

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;

  for(int iteration = 1;; iteration++)
  { 
    iterations++;
    logger.info("Iteration %i.", iteration);
    num_coarse++;

    solver_implicit.solve();
    Solution<double>::vector_to_solution(solver_implicit.get_sln_vector(), const_space, previous_mean_values);

    if(polynomialDegree)
    {
      num_fine++;
      solver_explicit.solve();
      double* merged_sln = merge_slns(solver_implicit.get_sln_vector(), const_space,solver_explicit.get_sln_vector(), space, full_space);
      Solution<double>::vector_to_solution(solver_explicit.get_sln_vector(), space, previous_derivatives);
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
      delete [] merged_sln;
    }
    else
    {
      solution->copy(previous_mean_values);
    }

    solution_view->show(solution);

    double err = calc_l2_error(solvedExample, mesh, solution, exact_solution, logger);

    if(err < wrt_exact_solver_tolerance + exact_solver_error)
      break;
  }

  logger.info("Iterations: %i", iterations);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("Fine systems solved: %i", num_fine);
}

static int smoothing_steps_count(int level, bool pre)
{
  if(level == 2)
    return 2;
  else
    return 5;
}

bool show_intermediate = false;

static double residual_drop(int level, bool pre)
{
  return 1e-1;
}

bool residual_condition(CSCMatrix<double>* mat, SimpleVector<double>* vec, double* sln_vector, double* residual, Hermes::Mixins::Loggable& logger, int iteration, bool output)
{
  mat->multiply_with_vector(sln_vector, residual, true);
  for(int i = 0; i < mat->get_size(); i++)
    residual[i] -= vec->get(i);

  if(output)
  {
    double residual_norm = Hermes2D::get_l2_norm(residual, mat->get_size());
    logger.info("\tIteration: %i, residual norm: %f.", iteration, residual_norm);
  }

  return false;
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                 double diffusivity, double time_step_length, 
                 double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                 ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger)
{

  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_2;
  SimpleVector<double> vector_b_2;
  CSCMatrix<double> matrix_A_1;
  SimpleVector<double> vector_b_1;
  CSCMatrix<double> matrix_A_0;
  SimpleVector<double> vector_b_0;

  // Matrices (M+A_tilde), vectors -A(u_K)
  SmoothingWeakForm weakform_smoother(solvedExample, true, 1, true, "Inlet", diffusivity, s, sigma);
  SmoothingWeakForm weakform_smoother_coarse(solvedExample, false, 1, true, "Inlet", diffusivity, s, sigma);
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  weakform_smoother_coarse.set_current_time_step(time_step_length);
  weakform_smoother_coarse.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2(ndofs_2);
  CSCMatrix<double> matrix_MA_tilde_1;
  SimpleVector<double> vector_A_1(ndofs_1);
  CSCMatrix<double> matrix_MA_0;
  SimpleVector<double> vector_A_0(ndofs_0);

  // Assembler.
  DiscreteProblem<double> dp;
  // Level 2.
  dp.set_space(space_2);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_2, &vector_b_2);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_2);

  // Level 1.
  dp.set_space(space_1);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_1, &vector_b_1);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_1);

  // Level 0.
  dp.set_space(space_0);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_0, &vector_b_0);
  dp.set_weak_formulation(&weakform_smoother_coarse);
  dp.assemble(&matrix_MA_0);

  UMFPackLinearMatrixSolver<double> solver_2(&matrix_MA_tilde_2, &vector_A_2);
  solver_2.setup_factorization();
  solver_2.set_reuse_scheme(MatrixStructureReuseScheme::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_1(&matrix_MA_tilde_1, &vector_A_1);
  solver_1.setup_factorization();
  solver_1.set_reuse_scheme(MatrixStructureReuseScheme::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_0(&matrix_MA_0, &vector_A_0);
  solver_0.setup_factorization();
  solver_0.set_reuse_scheme(MatrixStructureReuseScheme::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  double* residual_2 = new double[ndofs_2];
  SimpleVector<double> sln_2(ndofs_2);
  double* residual_1 = new double[ndofs_1];
  SimpleVector<double> sln_1(ndofs_1);
  double* residual_0 = new double[ndofs_0];
  SimpleVector<double> sln_0(ndofs_0);

  double initial_residual_norm;
  // Reports.
  int num_coarse = 0;
  int num_2 = 0;
  int num_1 = 0;
  int v_cycles = 0;

  for(int step = 1;; step++)
  { 
    logger.info("V-cycle %i.", step);
    v_cycles++;

#pragma region 0 - highest level
    // Store the previous solution.
    OGProjection<double>::project_global(space_2, previous_sln, &sln_2);
    if(polynomialDegree > 1)
    {
      for(int iteration = 1; iteration <= smoothing_steps_count(2, true); iteration++)
      {
        // Solve for increment.
        matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
        vector_A_2.change_sign()->add_vector(&vector_b_2);
        solver_2.solve();
        sln_2.add_vector(solver_2.get_sln_vector());

        // Make solution
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->show(previous_sln);

        // Residual check.
        residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger, iteration, true);
      }
    }
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    OGProjection<double>::project_global(space_1, previous_sln, &sln_1);

    // f_P1
    SimpleVector<double> f_P1(ndofs_1);
    f_P1.zero();
    // Minus A_P1
    SimpleVector<double> R_P1(ndofs_1);
    // Minus(minus) projected_A_P1
    SimpleVector<double> projected_A_2(ndofs_2);
    matrix_A_2.multiply_with_vector(sln_2.v, projected_A_2.v, true);

    SimpleVector<double>* projected_A_P_1
      = (SimpleVector<double>*)cut_off_quadratic_part(projected_A_2.v, space_1, space_2);

    SimpleVector<double>* sln_2_projected = cut_off_quadratic_part(sln_2.v, space_1, space_2);
    matrix_A_1.multiply_with_vector(sln_2_projected->v, R_P1.v, true);

    R_P1.change_sign();
    f_P1.add_vector(&R_P1);
    f_P1.add_vector(projected_A_P_1);
    f_P1.change_sign();

    for(int iteration = 1; iteration <= smoothing_steps_count(1, true); iteration++)
    {
      // Solve for increment.
      SimpleVector<double>* rhs;
      if(iteration == 1)
      {
        if(polynomialDegree > 1)
        {
          rhs = cut_off_quadratic_part(residual_2, space_1, space_2);
          memcpy(vector_A_1.v, rhs->v, ndofs_1 * sizeof(double));
        }
        else
        {
          matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);
          vector_A_1.change_sign()->add_vector(&vector_b_1);
        }
      }
      else
      {
        // A(u_K) - done after the first step.
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if(polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);
      }
      solver_1.solve();
      sln_1.add_vector(solver_1.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      if(show_intermediate)
        solution_view->show(previous_sln);

      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, false);
    }
#pragma endregion

#pragma region  2 - Solve the problem on the coarse level exactly
    OGProjection<double>::project_global(space_0, previous_sln, &sln_0);

    // f_P0
    SimpleVector<double> f_P0(ndofs_0);
    f_P0.zero();
    // Minus A_P0
    SimpleVector<double> R_P0(ndofs_0);
    // Minus(minus) projected_A_P0
    SimpleVector<double> projected_A_1(ndofs_1);
    matrix_A_1.multiply_with_vector(sln_1.v, projected_A_1.v, true);

    SimpleVector<double>* projected_A_P_0
      = (SimpleVector<double>*)cut_off_linear_part(projected_A_1.v, space_0, space_1);

    SimpleVector<double>* sln_1_projected = cut_off_linear_part(sln_1.v, space_0, space_1);
    matrix_A_0.multiply_with_vector(sln_1_projected->v, R_P0.v, true);

    SimpleVector<double> projected_f_P1(ndofs_1);
    projected_f_P1.set_vector(&f_P1);

    R_P0.change_sign();
    f_P0.add_vector(&R_P0);
    f_P0.add_vector(projected_A_P_0);
    f_P0.add_vector(cut_off_linear_part(projected_f_P1.v, space_0, space_1)->change_sign());
    f_P0.change_sign();

    num_coarse++;

    for(int iteration = 1; iteration <= smoothing_steps_count(0, true); iteration++)
    {
      SimpleVector<double>* rhs;
      if(iteration == 1)
      {
        if(polynomialDegree > 1)
        {
          rhs = (SimpleVector<double>*)cut_off_linear_part(residual_1, space_0, space_1)->add_vector(cut_off_linear_part(projected_f_P1.v, space_0, space_1)->change_sign());
          memcpy(vector_A_0.v, rhs->v, ndofs_0 * sizeof(double));
        }
        else
        {
          matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);
          vector_A_0.change_sign()->add_vector(&vector_b_0);
        }
      }
      else
      {
        // A(u_K) - done after the first step.
        matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);

        if(polynomialDegree > 1)
          vector_A_0.change_sign()->add_vector(&f_P0)->add_vector(&vector_b_0);
        else
          vector_A_0.change_sign()->add_vector(&vector_b_0);
      }

      solver_0.solve();
      sln_0.add_vector(solver_0.get_sln_vector());

      if(show_intermediate)
      {
        // Make solution
        Solution<double>::vector_to_solution(&sln_0, space_0, previous_sln);
        solution_view->show(previous_sln);
      }

      // Residual check.
      residual_condition(&matrix_A_0, &vector_b_0, sln_0.v, residual_0, logger, iteration, false);
    }
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    sln_1.set_vector(merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, false));
    Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);

    for(int iteration = 1; iteration <= smoothing_steps_count(1, false); iteration++)
    {
      // Solve for increment.
      matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

      if(polynomialDegree > 1)
        vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
      else
        vector_A_1.change_sign()->add_vector(&vector_b_1);

      solver_1.solve();
      sln_1.add_vector(solver_1.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      if(show_intermediate)
        solution_view->show(previous_sln);

      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, false);
    }
#pragma endregion

#pragma region 0 - highest level

    if(polynomialDegree > 1)
    {
      // Store the previous solution.
      sln_2.set_vector(merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false));
      Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);

      for(int iteration = 1; iteration <= smoothing_steps_count(2, false); iteration++)
      {
        // Solve for increment.
        matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
        vector_A_2.change_sign()->add_vector(&vector_b_2);
        solver_2.solve();
        sln_2.add_vector(solver_2.get_sln_vector());

        // Make solution
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->show(previous_sln);

        // Residual check.
        residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger, iteration, true);
      }
    }
#pragma endregion

    // Error & exact solution display.
    solution_view->show(previous_sln);
    if(calc_l2_error(solvedExample, mesh, previous_sln, exact_solution, logger) < wrt_exact_solver_tolerance + exact_solver_error)
      break;
  }
  logger.info("V-cycles: %i", v_cycles);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("1 systems solved: %i", num_1);
  logger.info("2 systems solved: %i", num_2);
}
