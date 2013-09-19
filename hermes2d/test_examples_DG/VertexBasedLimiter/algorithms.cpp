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
  LinearSolver<double> solver_explicit(&weakform_explicit, space);

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

    if(std::abs(exact_solver_error - err) < wrt_exact_solver_tolerance)
      break;
  }

  logger.info("Iterations: %i", iterations);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("Fine systems solved: %i", num_fine);
}

static int smoothing_steps_count(int level, bool pre)
{
  return 5;
}

static double residual_drop(int level, bool pre)
{
  return 1e-1;
}

static bool use_residual_drop_condition = true;
static bool residual_condition(CSCMatrix<double>* mat, SimpleVector<double>* vec, double* sln_vector, double* residual, Hermes::Mixins::Loggable& logger, int iteration, bool presmoothing)
{
  mat->multiply_with_vector(sln_vector, residual, true);
  for(int i = 0; i < mat->get_size(); i++)
    residual[i] -= vec->get(i);
  
  if(use_residual_drop_condition)
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
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2;
  CSCMatrix<double> matrix_MA_tilde_1;
  SimpleVector<double> vector_A_1;
  
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

  // Utils.
  double* residual_2 = new double[ndofs_2];
  SimpleVector<double> sln_2(ndofs_2);
  double* residual_1 = new double[ndofs_1];
  SimpleVector<double> sln_1(ndofs_1);
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
        dp.set_space(space_2);
        dp.set_weak_formulation(&weakform_smoother);
        dp.assemble(&vector_A_2);
        UMFPackLinearMatrixSolver<double> solver(&matrix_MA_tilde_2, (SimpleVector<double>*)vector_A_2.add_vector(&vector_b_2));
        solver.solve();
        sln_2.add_vector(solver.get_sln_vector());

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

    for(int iteration = 1; iteration <= smoothing_steps_count(1, true); iteration++)
    {
      // Solve for increment.
      SimpleVector<double>* rhs;
      if(iteration == 1)
        rhs = cut_off_quadratic_part(residual_2, space_1, space_2);
      else
      {
        dp.set_space(space_1);
        dp.set_weak_formulation(&weakform_smoother);
        dp.assemble(&vector_A_1);
        SimpleVector<double> sln_2_projected(ndofs_1);
        sln_2_projected.set_vector(cut_off_quadratic_part(sln_2.v, space_1, space_2));
        double* projected_residual_2 = new double[ndofs_1];
        matrix_A_1.multiply_with_vector(sln_2_projected.v, projected_residual_2, true);
        rhs = (SimpleVector<double>*)vector_A_1.add_vector(projected_residual_2);
        delete [] projected_residual_2;
      }
      UMFPackLinearMatrixSolver<double> solver(&matrix_MA_tilde_1, rhs);
      solver.solve();
      sln_1.add_vector(solver.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      solution_view->show(previous_sln);
        
      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, true);
    }
#pragma endregion

#pragma region  2 - Solve the problem on the coarse level exactly
    num_coarse++;
    OGProjection<double>::project_global(space_0, previous_sln, &sln_0);

    //// Solve for increment
    // First term.
    SimpleVector<double>* vector = cut_off_linear_part(residual_1, space_0, space_1);
    // Second term (two projections of 2-residual).
    vector->add_vector(cut_off_linear_part(residual_2, space_0, space_2));
    
    // Third - complicated.
    double* projected_residual_2 = new double[ndofs_1];

    // I(u_2)
    SimpleVector<double> sln_2_projected(ndofs_1);
    sln_2_projected.set_vector(cut_off_quadratic_part(sln_2.v, space_1, space_2));
    // A(I(u_2))
    matrix_A_1.multiply_with_vector(sln_2_projected.v, projected_residual_2, true);
    // -A(I(u_2)) + b_1 (= R(I(u_2))
    SimpleVector<double> projected_residual(ndofs_1);
    projected_residual.set_vector(projected_residual_2);
    delete [] projected_residual_2;
    projected_residual.change_sign();
    projected_residual.add_vector(&vector_b_1);
    // Add I(R(I(u_2)))
    vector->add_vector(cut_off_linear_part(projected_residual.v, space_0, space_1));
    vector->change_sign();

    //cut_off_linear_part(residual_2, space_0, space_2)->export_to_file("a", "a", EXPORT_FORMAT_PLAIN_ASCII);
    //cut_off_linear_part(projected_residual.v, space_0, space_1)->export_to_file("b", "a", EXPORT_FORMAT_PLAIN_ASCII);
    UMFPackLinearMatrixSolver<double> solver(&matrix_A_0, vector);
    solver.solve();
    delete [] vector;
    sln_0.set_vector(solver.get_sln_vector());
    // Make solution
    Solution<double>::vector_to_solution(&sln_0, space_0, previous_sln);
    solution_view->show(previous_sln);
    //View::wait_for_keypress();
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    sln_1.set_vector(merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, true));
    Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);

    for(int iteration = 1; iteration <= smoothing_steps_count(1, true); iteration++)
    {
      // Solve for increment.
      SimpleVector<double>* rhs;
      dp.set_space(space_1);
      dp.set_weak_formulation(&weakform_smoother);
      dp.assemble(&vector_A_1);
      SimpleVector<double> sln_2_projected(ndofs_1);
      sln_2_projected.set_vector(cut_off_quadratic_part(sln_2.v, space_1, space_2));
      double* projected_residual_2 = new double[ndofs_1];
      matrix_A_1.multiply_with_vector(sln_2_projected.v, projected_residual_2, true);
      rhs = (SimpleVector<double>*)vector_A_1.add_vector(projected_residual_2);
      delete [] projected_residual_2;
      UMFPackLinearMatrixSolver<double> solver(&matrix_MA_tilde_1, rhs);
      solver.solve();
      sln_1.add_vector(solver.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      solution_view->show(previous_sln);
        
      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, true);
    }
#pragma endregion

    #pragma region 0 - highest level
    // Store the previous solution.
    sln_2.set_vector(merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false));
    Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);

    if(polynomialDegree > 1)
    {
      for(int iteration = 1; iteration <= smoothing_steps_count(2, true); iteration++)
      {
        // Solve for increment.
        dp.set_space(space_2);
        dp.set_weak_formulation(&weakform_smoother);
        dp.assemble(&vector_A_2);
        UMFPackLinearMatrixSolver<double> solver(&matrix_MA_tilde_2, (SimpleVector<double>*)vector_A_2.add_vector(&vector_b_2));
        solver.solve();
        sln_2.add_vector(solver.get_sln_vector());

        // Make solution
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->show(previous_sln);
        
        // Residual check.
        residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger, iteration, true);
      }
    }
#pragma endregion

    // Error & exact solution display.
    if(std::abs(exact_solver_error - calc_l2_error(solvedExample, mesh, previous_sln, exact_solution, logger)) < wrt_exact_solver_tolerance)
      break;
  }
  logger.info("V-cycles: %i", v_cycles);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("1 systems solved: %i", num_1);
  logger.info("2 systems solved: %i", num_2);
}
