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
  OGProjection<double>::project_global(space, exact_solution, exact_solver_sln);
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
  return 3;
}

static double residual_drop(int level, bool pre)
{
  return 1e-1;
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                 double diffusivity, double time_step_length, 
                 double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                 ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger)
{
  bool use_residual_drop_condition = true;

  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // 1 - solver
  SmoothingWeakForm weakform_1(solvedExample, true, 1, true, "Inlet", diffusivity, s, sigma, false);
  SmoothingWeakForm weakform_2(solvedExample, true, 1, true, "Inlet", diffusivity, s, sigma);
  weakform_1.set_current_time_step(time_step_length);
  weakform_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_solution, exact_solution));
  weakform_2.set_current_time_step(time_step_length);
  weakform_2.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  LinearSolver<double> solver_1(&weakform_1, space_1);
  LinearSolver<double> solver_2(&weakform_2, space_2);
  DiscreteProblem<double> solver_dp_1(&weakform_1, space_1);
  CSCMatrix<double> matrix_1;
  solver_1.set_verbose_output(false);
  solver_1.set_jacobian_constant();
  solver_2.set_verbose_output(false);
  solver_2.set_jacobian_constant();

  // 1 - Residual measurement.
  SmoothingWeakFormResidual weakform_residual_1(solvedExample, 1, true, "Inlet", diffusivity, s, sigma);
  SmoothingWeakFormResidual weakform_residual_2(solvedExample, 1, true, "Inlet", diffusivity, s, sigma);
  weakform_residual_1.set_current_time_step(time_step_length);
  weakform_residual_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_solution, exact_solution));
  weakform_residual_2.set_current_time_step(time_step_length);
  weakform_residual_2.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));

  DiscreteProblem<double> dp_1(&weakform_residual_1, space_1);
  Algebra::SimpleVector<double> vec_1;
  DiscreteProblem<double> dp_2(&weakform_residual_2, space_2);
  Algebra::SimpleVector<double> vec_2;

  // 0 - solver
  FullImplicitWeakForm weakform_0(solvedExample, 1, true, "Inlet", diffusivity);
  weakform_0.set_current_time_step(time_step_length);
  weakform_0.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  DiscreteProblem<double> dp_0(&weakform_0, space_0);
  CSCMatrix<double> matrix_0;

  // Utils.
  double* slnv_2 = new double[ndofs_2];
  double* slnv_1 = new double[ndofs_1];
  double* slnv_1_increment = new double[ndofs_1];
  double* slnv_0 = new double[ndofs_0];
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
    OGProjection<double>::project_global(space_2, previous_sln, slnv_2);
    if(polynomialDegree > 1)
    {
      for(int iteration = 1; iteration <= smoothing_steps_count(2, true); iteration++)
      {
        num_2++;
        // Solve for increment.
        solver_2.solve();
        // Add
        for(int k = 0; k < ndofs_2; k++)
          slnv_2[k] += solver_2.get_sln_vector()[k];

        Solution<double>::vector_to_solution(slnv_2, space_2, previous_sln);
        
        // Show
        solution_view->show(previous_sln);
        
        // Residual check.
        dp_2.assemble(&vec_2);
        if(use_residual_drop_condition)
        {
          double residual_norm = Hermes2D::get_l2_norm(&vec_2);
          logger.info("\tIteration - (P = 2-pre): %i, residual norm: %f.", iteration, residual_norm);
          if(iteration == 1)
            initial_residual_norm = residual_norm;
          else if(residual_norm / initial_residual_norm < residual_drop(2, true))
            break;
        }
      }
    }
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    memset(slnv_1_increment, 0, sizeof(double)*space_1->get_num_dofs());
    
    for(int iteration = 1; iteration <= smoothing_steps_count(1, true); iteration++)
    {
      num_1++;

      // Store the previous solution.
      OGProjection<double>::project_global(space_1, previous_sln, slnv_1);
      if(polynomialDegree > 1 && iteration == 1)
      {
        if(step == 1)
          solver_dp_1.assemble(&matrix_1);
        UMFPackLinearMatrixSolver<double> local_solver_1(&matrix_1, (Algebra::SimpleVector<double>*)cut_off_quadratic_part(&vec_2, space_1, space_2));
        local_solver_1.solve();
        for(int k = 0; k < ndofs_1; k++)
        {
          slnv_1_increment[k] += local_solver_1.get_sln_vector()[k];
          slnv_1[k] += local_solver_1.get_sln_vector()[k];
        }
      }
      else
      {
        solver_1.solve();
        for(int k = 0; k < ndofs_1; k++)
        {
          slnv_1_increment[k] += solver_1.get_sln_vector()[k];
          slnv_1[k] += solver_1.get_sln_vector()[k];
        }
      }

      Solution<double>::vector_to_solution(slnv_1, space_1, previous_sln);

      // Show
      solution_view->show(previous_sln);

      // Residual check.
      dp_1.assemble(&vec_1);
      if(use_residual_drop_condition)
      {
        double residual_norm = Hermes2D::get_l2_norm(&vec_1);
        logger.info("\tIteration - (P = 1-pre): %i, residual norm: %f.", iteration, residual_norm);
        if(iteration == 1)
          initial_residual_norm = residual_norm;
        else if(residual_norm / initial_residual_norm < residual_drop(1, true))
          break;
      }
    }
#pragma endregion

#pragma region  2 - Solve the problem on the coarse level exactly
    num_coarse++;
    //// Solve for increment
    if(step == 1)
      dp_0.assemble(&matrix_0);
    UMFPackLinearMatrixSolver<double> solver_0(&matrix_0, (Algebra::SimpleVector<double>*)cut_off_linear_part(&vec_1, space_0, space_1));
    solver_0.solve();
    //// Show
    ////coarse_solution_view.show(previous_sln);
#pragma endregion

#pragma region 1 - intermediate level
    // 3 - Prolongation and replacement
    double* solution_vector = merge_slns(solver_0.get_sln_vector(), space_0, slnv_1, space_1, space_1, true);
    Solution<double>::vector_to_solution(solution_vector, space_1, previous_sln);
    
    // Store the previous solution.
    memcpy(slnv_1, solution_vector, sizeof(double)*space_1->get_num_dofs());
    delete [] solution_vector;
    for(int iteration = 1; iteration <= smoothing_steps_count(1, false); iteration++)
    {
      num_1++;

      solver_1.solve();
      
      for(int k = 0; k < ndofs_1; k++)
      {
        slnv_1_increment[k] += solver_1.get_sln_vector()[k];
        slnv_1[k] += solver_1.get_sln_vector()[k];
      }

      Solution<double>::vector_to_solution(slnv_1, space_1, previous_sln);

      // Show
      solution_view->show(previous_sln);

      // Residual check.
      dp_1.assemble(&vec_1);
      if(use_residual_drop_condition)
      {
        double residual_norm = Hermes2D::get_l2_norm(&vec_1);
        logger.info("\tIteration - (P = 1-post): %i, residual norm: %f.", iteration, residual_norm);
        if(iteration == 1)
          initial_residual_norm = residual_norm;
        else if(residual_norm / initial_residual_norm < residual_drop(1, false))
          break;
      }
    }
#pragma endregion


#pragma region 0 - highest level
    if(polynomialDegree > 1)
    {
      // 4 - Prolongation and replacement
      solution_vector = merge_slns(slnv_1, space_1, slnv_2, space_2, space_2);
      Solution<double>::vector_to_solution(solution_vector, space_2, previous_sln);
      // Store the previous solution.
      memcpy(slnv_2, solution_vector, sizeof(double)*space_2->get_num_dofs());
      delete [] solution_vector;
    
      for(int iteration = 1; iteration <= smoothing_steps_count(2, false); iteration++)
      {
        num_2++;
       
        // Solve for increment.
        solver_2.solve();

        // Add
        for(int k = 0; k < ndofs_2; k++)
          slnv_2[k] += solver_2.get_sln_vector()[k];
        
        Solution<double>::vector_to_solution(slnv_2, space_2, previous_sln);

        // Show
        solution_view->show(previous_sln);
      
        // Residual check.
        dp_2.assemble(&vec_2);
        if(use_residual_drop_condition)
        {
          double residual_norm = Hermes2D::get_l2_norm(&vec_2);
          logger.info("\tIteration - (P = 2-post): %i, residual norm: %f.", iteration, residual_norm);
          if(iteration == 1)
            initial_residual_norm = residual_norm;
          else if(residual_norm / initial_residual_norm < residual_drop(2, false))
            break;
        }
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
