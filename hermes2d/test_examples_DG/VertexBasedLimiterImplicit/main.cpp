#include "definitions.h"
#include "../euler_util.h"

/// Init.
const int polynomialDegree = 1;
const int initialRefinementsCount = 5;
const SolvedExample solvedExample = AdvectedCube;

// Limiting.
bool SHOCK_CAPTURING = true;
const EulerLimiterType limiter_type = VertexBased;

// Views.
bool HermesView = true;
bool VTKView = false;

// Logging.
const double logPercentTimeSteps = .01;
const double time_step_length = 0.01;
const double time_interval_length = solvedExample == SolidBodyRotation ? 2 * M_PI : 1.;
int logPeriod = (int)std::max<double>(1., ((logPercentTimeSteps / 100.) * (time_interval_length / time_step_length)));
Hermes::Mixins::Loggable logger(true);

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load((solvedExample == SolidBodyRotation ? "domain_rotation.xml" : "larger_domain.xml"), mesh);

  for(int i = 0; i < initialRefinementsCount; i++)
    mesh->refine_all_elements();

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* initial_condition;
  if(solvedExample == SolidBodyRotation)
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
  else if (solvedExample == AdvectedCube)
    initial_condition = new InitialConditionAdvectedCube(mesh);

  ExactSolutionScalar<double>* previous_initial_condition;
  if(solvedExample == SolidBodyRotation)
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
  else if (solvedExample == AdvectedCube)
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);

  MeshFunctionSharedPtr<double>previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double>previous_not_updated_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double>previous_solution_time_step(initial_condition);
  // Visualization.
  ScalarView solution_view("Initial condition", new WinGeom(520, 10, 500, 500));
  Linearizer lin;
  if(HermesView)
    solution_view.show(previous_solution);

  // Weak form.
  CustomWeakForm weakform_implicit(solvedExample, Implicit, 0);
  CustomWeakForm weakform_explicit(solvedExample, Explicit, 1);
  weakform_implicit.set_ext(previous_solution_time_step);
  weakform_implicit.set_current_time_step(time_step_length);

  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_not_updated_solution, previous_solution_time_step));
  weakform_explicit.set_current_time_step(time_step_length);

  // Solver.
  LinearSolver<double> solver_implicit(&weakform_implicit, const_space);
  LinearSolver<double> solver_explicit(&weakform_explicit, space);

  Hermes::Mixins::Loggable::set_static_logFile_name("logfile.h2d");

  // Solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);

  double current_time = 0.;
  Element* e;
  AsmList<double> al_fine;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
  for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    if((!(time_step % logPeriod)) || (time_step == number_of_steps))
    {
      logger.info("Time step: %i, time: %f.", time_step, current_time);
    }

    // 0th - step
    double* previous_sln_vector = new double[ndofs];
    double* sln_vector;
    double* mean_values;
    OGProjection<double>::project_global(space, previous_solution_time_step, previous_sln_vector);

    // 1 - solve implicit.
    solver_implicit.solve();
    mean_values = new double[ndofs];
    Solution<double>::vector_to_solution(solver_implicit.get_sln_vector(), const_space, previous_solution_time_step);
    OGProjection<double>::project_global(space, previous_solution_time_step, mean_values);

    // 2 - Update the mean values.
    for_all_active_elements(e, mesh)
    {
      space->get_element_assembly_list(e, &al_fine);
      for(unsigned int shape_i = 0; shape_i < al_fine.cnt; shape_i++)
      {
        int order = space->get_shapeset()->get_order(al_fine.idx[shape_i], e->get_mode());
        if(order == 0)
        {
          int dof = al_fine.dof[shape_i];
          previous_sln_vector[dof] = mean_values[dof];
        }
      }
    }
    
    // Clean up.
    delete [] mean_values;
    
    // 3 - Solve explicit.
    Solution<double>::vector_to_solution(previous_sln_vector, space, previous_solution_time_step);
    delete [] previous_sln_vector;
    solver_explicit.solve();
    sln_vector = solver_explicit.get_sln_vector();

#pragma region *. Get the solution with optional shock capturing.
  if(!SHOCK_CAPTURING)
    Solution<double>::vector_to_solution(sln_vector, space, solution);
  else
  {
    PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, space, sln_vector, 1);
    solution->copy(limiter->get_solution());
    delete limiter;
  }
#pragma endregion

  if((!(time_step % logPeriod)) || (time_step == number_of_steps))
  {
    if(HermesView)
    {
      solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
      solution_view.show(solution);
      //View::wait_for_keypress();
    }
    if(VTKView)
    {
      char* filename = new char[100];
      sprintf(filename, "Sln-%i.vtk", time_step);
      lin.save_solution_vtk(solution, filename, "sln");
      delete [] filename;
    }
  }

  previous_solution_time_step->copy(solution);
  previous_not_updated_solution->copy(solution);
  
  current_time += time_step_length;
}

View::wait();
return 0;
}