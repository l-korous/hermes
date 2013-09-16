#include "algorithms.h"

static void calc_l2_error(MeshSharedPtr mesh, MeshFunctionSharedPtr<double> fn_1, MeshFunctionSharedPtr<double> fn_2)
{
  ErrorWeakForm wf;
  SpaceSharedPtr<double> mspace(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  wf.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(fn_1, fn_2));
  DiscreteProblem<double> dp(&wf, mspace);
  UMFPackVector<double> vector;
  dp.assemble(&vector);
  double result = 0.;
  for(int i = 0; i < vector.length(); i++)
    result += vector.get(i);
  std::cout << "L2 Error: " << std::sqrt(result) << std::endl;
}

static void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step)
{
  MeshFunctionSharedPtr<double> exact_solver_sln(new Solution<double>());
  ScalarView* exact_solver_view = new ScalarView("Exact solver solution", new WinGeom(0, 720, 600, 350));

  // Exact solver
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", "Outlet", diffusivity, s, sigma);
  weakform_exact.set_current_time_step(time_step);
  weakform_exact.set_ext(initial_sln);
  LinearSolver<double> solver_exact(&weakform_exact, space);
  solver_exact.solve();
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, exact_solver_sln);
  //OGProjection<double>::project_global(space, exact_solution, exact_solver_sln);
  calc_l2_error(space->get_mesh(), exact_solver_sln, exact_solution);
  exact_solver_view->show(exact_solver_sln);
}

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_mean_values, 
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, 
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);
  
  MeshFunctionSharedPtr<double>initial_sln(new InitialConditionBenchmark(mesh, diffusivity));

  ImplicitWeakForm weakform_implicit(solvedExample, true, "Inlet", "Outlet", diffusivity, s, sigma);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, initial_sln));
  ExplicitWeakForm weakform_explicit(solvedExample, true, "Inlet", "Outlet", diffusivity, s, sigma);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, initial_sln));
  weakform_explicit.set_current_time_step(time_step_length);
  LinearSolver<double> solver_implicit(&weakform_implicit, const_space);
  LinearSolver<double> solver_explicit(&weakform_explicit, space);
  
  solve_exact(solvedExample, full_space, diffusivity, s, sigma, exact_solution, initial_sln, time_step_length);
  exact_view->show(exact_solution);

  double current_time = 0.;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
  for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    std::cout << "Iteration " << time_step << std::endl;
    solver_implicit.solve();
    Solution<double>::vector_to_solution(solver_implicit.get_sln_vector(), const_space, previous_mean_values);

    if(polynomialDegree)
    {
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
    
    calc_l2_error(mesh, solution, exact_solution);

    /*
    {
    Linearizer lin;
    char* filename = new char[100];
    sprintf(filename, "Sln-%i.vtk", time_step);
    lin.save_solution_vtk(solution, filename, "sln");
    delete [] filename;
    }
    */
    current_time += time_step_length;
  }
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                 double diffusivity, double time_step_length, 
                 double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                 ScalarView* solution_view, ScalarView* exact_view, double s, double sigma)
{
  // Spaces
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Previous iteration solution
  MeshFunctionSharedPtr<double>initial_sln(new InitialConditionBenchmark(mesh, diffusivity));
  ScalarView coarse_solution_view("Coarse solution", new WinGeom(0, 360, 600, 350));

  // 1 - solver
  SmoothingWeakForm weakform_1(solvedExample, true, 1, true, "Inlet", "Outlet", diffusivity, s, sigma);
  weakform_1.set_current_time_step(time_step_length);
  weakform_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, initial_sln));
  LinearSolver<double> solver_1(&weakform_1, space_1);
  solver_1.set_verbose_output(false);

  // 1 - Residual measurement.
  SmoothingWeakFormResidual weakform_residual_1(solvedExample, 1, true, "Inlet", "Outlet", diffusivity, s, sigma);
  weakform_residual_1.set_current_time_step(time_step_length);
  weakform_residual_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, initial_sln));
  DiscreteProblem<double> dp(&weakform_residual_1, space_1);
  Algebra::UMFPackVector<double> vec;

  // 0 - solver
  FullImplicitWeakForm weakform_0(solvedExample, 1, true, "Inlet", "Outlet", diffusivity);
  weakform_0.set_current_time_step(time_step_length);
  weakform_0.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, initial_sln));
  DiscreteProblem<double> dp_0(&weakform_0, space_0);
  UMFPackMatrix<double> matrix;
  
  // Exact solver.
  solve_exact(solvedExample, space_1, diffusivity, s, sigma, exact_solution, initial_sln, time_step_length);

  // Utils.
  double* slnv_1 = new double[ndofs_1];
  double* slnv_0 = new double[ndofs_0];
  double initial_residual_norm;
    
  double current_time = 0.;
  for(int step = 0;; step++)
  { 
    std::cout << "V-cycle " << step << std::endl;

    // 1 - pre-smoothing on the 1st level.
    for(int iteration_1 = 1; iteration_1 < 5; iteration_1++)
    {
      // Store the previous solution.
      OGProjection<double>::project_global(space_1, previous_sln, slnv_1);
      
      // Solve for increment.
      solver_1.solve();
      
      // Add
      for(int k = 0; k < ndofs_1; k++)
        slnv_1[k] += solver_1.get_sln_vector()[k];
      Solution<double>::vector_to_solution(slnv_1, space_1, previous_sln);
      
      // Show
      solution_view->show(previous_sln);
      
      // Residual check.
      dp.assemble(&vec);
      double residual_norm = Hermes2D::get_l2_norm(&vec);
      std::cout << "\tIteration - (P = 1-pre): " << iteration_1 << ", residual norm: " << residual_norm << std::endl;
      if(iteration_1 == 1)
        initial_residual_norm = residual_norm;
      else if(residual_norm / initial_residual_norm < 1e-1)
        break;
    }

    // 2 - Solve the problem on the coarse level exactly
    //// Store
    OGProjection<double>::project_global(space_0, previous_sln, slnv_0);
    
    //// Solve for increment
    dp_0.assemble(&matrix);
    UMFPackLinearMatrixSolver<double> solver_0(&matrix, (Algebra::UMFPackVector<double>*)cut_off_linear_part(&vec, space_0, space_1));
    solver_0.solve();
    //// Add
    for(int k = 0; k < ndofs_0; k++)
      slnv_0[k] += solver_0.get_sln_vector()[k];
    Solution<double>::vector_to_solution(slnv_0, space_0, previous_sln);
    //// Show
    coarse_solution_view.show(previous_sln);
    
    // 3 - Prolongation and replacement
    double* solution_vector = merge_slns(slnv_0, space_0, slnv_1, space_1, space_1);
    Solution<double>::vector_to_solution(solution_vector, space_1, previous_sln);
    delete [] solution_vector;

    // 4 - post-smoothing steps
    for(int iteration_1 = 1; iteration_1 < 5; iteration_1++)
    {
      // Store the previous solution.
      OGProjection<double>::project_global(space_1, previous_sln, slnv_1);
      
      // Solve for increment.
      solver_1.solve();
      
      // Add
      for(int k = 0; k < ndofs_1; k++)
        slnv_1[k] += solver_1.get_sln_vector()[k];
      Solution<double>::vector_to_solution(slnv_1, space_1, previous_sln);
      
      // Show
      solution_view->show(previous_sln);
      
      // Residual check.
      dp.assemble(&vec);
      double residual_norm = Hermes2D::get_l2_norm(&vec);
      std::cout << "\tIteration - (P = 1-post): " << iteration_1 << ", residual norm: " << residual_norm << std::endl;
      if(iteration_1 == 1)
        initial_residual_norm = residual_norm;
      else if(residual_norm / initial_residual_norm < 1e-1)
        break;
    }

    // Error & exact solution display.
    //((ExactSolutionMovingPeak*)(exact_solution.get()))->set_current_time(current_time + (M_PI / 2.));
    exact_view->show(exact_solution);
    calc_l2_error(mesh, previous_sln, exact_solution);

    /*
    {
    Linearizer lin;
    char* filename = new char[100];
    sprintf(filename, "Sln-%i.vtk", time_step);
    lin.save_solution_vtk(solution, filename, "sln");
    delete [] filename;
    }
    */
    current_time += time_step_length;
  }
}
