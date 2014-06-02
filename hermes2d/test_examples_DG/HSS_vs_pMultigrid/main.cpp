#include <climits>
#include "definitions.h"
#include "../euler_util.h"
#include "algorithms.h"

int polynomialDegree = 2;
int initialRefinementsCount = 4;
const Algorithm algorithm = Multiscale;
SolvedExample solvedExample = MovingPeak;
// For the initial shape of the peak.
const double MovingPeakDiffusivity = 1e-3;
double diffusivity = 1e-3;
double s = -1.;
double sigma_star = 10.;
double CFL = 1e-1;

Hermes::vector<int> smoothing_steps_per_V_cycle;

static std::string SolvedExampleString[5] = { "1D", "CircularConvection", "MovingPeak", "AdvectedCube", "SolidBodyRotation" };
std::string solvedExampleString = SolvedExampleString[solvedExample];
Hermes::Mixins::TimeMeasurable cpu_time;

// Arguments (optional)
// 1 - Initial refinements
// 2 - diffusivity
// 3 - example (index in the array SolvedExampleString)
// 4 - CFL number
int main(int argc, char* argv[])
{
  // Use one thread
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Read the input arguments.
  if (argc > 1)
    initialRefinementsCount = atoi(argv[1]);
  if (argc > 2)
    diffusivity = (double)atof(argv[2]);
  if (argc > 3)
  {
    solvedExample = (SolvedExample)atoi(argv[3]);
    solvedExampleString = SolvedExampleString[solvedExample];
  }
  if (argc > 4)
    CFL = atof(argv[4]);

  // Fix CFL to 128 if stationary.
  // Only 1 smoothing step per V-cycle for time-dependent
  if (!is_timedep(solvedExample))
  {
    CFL = 128.;
    smoothing_steps_per_V_cycle.push_back(2);
    smoothing_steps_per_V_cycle.push_back(3);
    smoothing_steps_per_V_cycle.push_back(5);
    smoothing_steps_per_V_cycle.push_back(10);
    smoothing_steps_per_V_cycle.push_back(15);
  }
  else
  {
    smoothing_steps_per_V_cycle.push_back(1);
  }

  // Load the mesh & set mesh size for CFL -> time step length calculation.
  double mesh_size;
  int time_step_count = INT_MAX;  // Infinite number of time steps for stationary problems
  double time_step_length = CFL * std::pow(2., -(double)initialRefinementsCount);
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  switch (solvedExample)
  {
  case SolidBodyRotation:
    mloader.load("domain_rotation.xml", mesh);
    for (int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 1.;
    time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
    time_step_count = 2 * M_PI / time_step_length;
    time_step_length = 2 * M_PI / time_step_count;
    break;
  case AdvectedCube:
    mloader.load("larger_domain.xml", mesh);
    for (int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 1.;
    time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
    break;
  case CircularConvection:
    mloader.load("domain_circular_convection.xml", mesh);
    for (int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 1.;
    time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
    break;
  case MovingPeak:
    mloader.load("domain.xml", mesh);
    for (int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 2.4;
    time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
    time_step_count = 2 * M_PI / time_step_length;
    time_step_length = 2 * M_PI / time_step_count;
    break;
  case Benchmark:
    mloader.load("domain_benchmark.xml", mesh);
    mesh->refine_all_elements(2);
    for (int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = .6;
    time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
    break;
  }

  // Sigma - needs to know initialRefinementsCount
  double sigma = std::pow(2., (double)(initialRefinementsCount)) * sigma_star * mesh_size;

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition_der = NULL;
  ExactSolutionScalar<double>* initial_solution = NULL;
  ExactSolutionScalar<double>* exact_sln = NULL;
  switch (solvedExample)
  {
  case SolidBodyRotation:
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
    initial_condition_der = new InitialConditionSolidBodyRotation(mesh);
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
    break;
  case AdvectedCube:
    initial_condition = new InitialConditionAdvectedCube(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);
    break;
  case CircularConvection:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    exact_sln = new ExactSolutionCircularConvection(mesh);
    initial_solution = new ExactSolutionCircularConvection(mesh);
    break;
  case MovingPeak:
    initial_condition = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    initial_condition_der = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    previous_initial_condition = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    exact_sln = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, (1. / 2.) * M_PI);
    initial_solution = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, (1. / 2.) * M_PI);
    break;
  case Benchmark:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    exact_sln = new ExactSolutionBenchmark2(mesh, diffusivity);
    initial_solution = new ExactSolutionBenchmark2(mesh, diffusivity);
    break;
  }

  // Solutions.
  MeshFunctionSharedPtr<double> previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double> previous_mean_values(initial_condition);
  MeshFunctionSharedPtr<double> previous_derivatives(initial_condition_der);
  MeshFunctionSharedPtr<double> exact_solution(exact_sln);
  MeshFunctionSharedPtr<double> initial_sln(initial_solution);


  // Visualization classes.
  ScalarView solution_view("Solution", new WinGeom(0, 0, 600, 350));
  ScalarView exact_view("Exact solution", new WinGeom(610, 0, 600, 350));

  // Global logger
  Hermes::Mixins::Loggable logger_global(true, NULL, true);
  logger_global.set_logFile_name(solvedExampleString.append(".h2d"));
  logger_global.set_timestamps(false);

  // Exact solver solution
  {
    SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
    cpu_time.tick();
    if (!is_timedep(solvedExample))
      solve_exact(solvedExample, space, diffusivity, s, sigma, exact_solution, initial_sln,
      time_step_length, polynomialDegree, initialRefinementsCount);
    else
      exact_solver_timedep(mesh, solvedExample, polynomialDegree, initialRefinementsCount, diffusivity, s, sigma, time_step_length, time_step_count, initial_sln, exact_solution, &exact_view, CFL);
    cpu_time.tick();
    std::stringstream ss_global;
    ss_global << "Grid=" << initialRefinementsCount << "|Diff=" << diffusivity << "|CFL=" << CFL << "|Exact|Time=" << cpu_time.last();
    logger_global.info(ss_global.str().c_str());
  }

  return 0;
  // HSS
  if (algorithm == Multiscale || algorithm == Both)
  {
    // Set up loggers.
    Hermes::Mixins::Loggable logger_HSS(true, NULL, false);
    logger_HSS.set_timestamps(false);
    logger_HSS.set_erase_on_beginning(true);
    logger_HSS.set_file_output_only(true);
    Hermes::Mixins::Loggable logger_details(true);
    logger_details.set_timestamps(false);
    logger_details.set_erase_on_beginning(true);
    logger_details.set_file_output_only(false);

    std::stringstream ss;
    ss << "HSS_" << solvedExampleString << "_" << initialRefinementsCount << "_" << diffusivity << "_CFL=" << CFL << ".h2d";
    logger_HSS.set_logFile_name(ss.str());
    std::stringstream ssd;
    ssd << "HSS_detail_" << solvedExampleString << "_" << initialRefinementsCount << "_" << diffusivity << ".h2d";
    logger_details.set_logFile_name(ssd.str());


    // Start measuring CPU time.
    cpu_time.tick();

    // Calculate & return what to put in the log.
    std::string outString = is_timedep(solvedExample) ?
      multiscale_decomposition_timedep(mesh, solvedExample, polynomialDegree, initialRefinementsCount,
      previous_mean_values, previous_derivatives, diffusivity, s, sigma,
      time_step_length, time_step_count, exact_solution, &solution_view, &exact_view, logger_HSS, logger_details, CFL)
      :
      multiscale_decomposition(mesh, solvedExample, polynomialDegree, initialRefinementsCount,
      previous_mean_values, previous_derivatives, diffusivity, s, sigma,
      time_step_length, time_step_count, exact_solution, &solution_view, &exact_view, logger_HSS, logger_details, CFL);

    // Stop measuring time.
    cpu_time.tick();


    // Fill the logs.
    logger_HSS.info("%f|%s", cpu_time.last(), outString.c_str());
    std::stringstream ss_global;
    ss_global << "Grid=" << initialRefinementsCount << "|Diff=" << diffusivity << "|CFL=" << CFL << "|HSS|Time=" << cpu_time.last() << "|" << outString;
    logger_global.info(ss_global.str().c_str());
  }

  // p-Multigrid
  if (algorithm == pMultigrid || algorithm == Both)
  {
    // Set up loggers.
    Hermes::Mixins::Loggable logger_pMultigrid(true, NULL, false);
    logger_pMultigrid.set_timestamps(false);
    logger_pMultigrid.set_erase_on_beginning(true);
    logger_pMultigrid.set_file_output_only(true);

    for (int si = 0; si < smoothing_steps_per_V_cycle.size(); si++)
    {
      std::stringstream ss;
      ss << "MG(" << smoothing_steps_per_V_cycle[si] << ")_" << solvedExampleString << "_" << initialRefinementsCount << "_" << diffusivity << "_CFL=" << CFL << ".h2d";
      logger_pMultigrid.set_logFile_name(ss.str());

      // This is here because of the loop over smoothing_steps_per_V_cycle.
      MeshFunctionSharedPtr<double> previous_solution_local(new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.));


      // Start measuring CPU time.
      cpu_time.tick();

      // Calculate & return what to put in the log.
      std::string outString =
        p_multigrid(mesh, solvedExample, polynomialDegree, initialRefinementsCount,
        previous_solution_local, diffusivity, time_step_length, time_step_count, exact_solution, &solution_view, &exact_view, s, sigma, logger_pMultigrid,
        smoothing_steps_per_V_cycle[si], CFL);

      // Stop measuring time.
      cpu_time.tick();


      // Fill the logs.
      logger_pMultigrid.info("%f|%s", cpu_time.last(), outString.c_str());
      std::stringstream ss_global;
      ss_global << "Grid=" << initialRefinementsCount << "|Diff=" << diffusivity << "|CFL=" << CFL
        << "|MG||Time=" << cpu_time.last() << "|" << outString;
      logger_global.info(ss_global.str().c_str());
    }
  }

  return 0;
}
