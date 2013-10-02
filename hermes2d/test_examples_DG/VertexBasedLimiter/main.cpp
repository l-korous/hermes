#include "definitions.h"
#include "../euler_util.h"
#include "algorithms.h"

const int polynomialDegree = 2;
int initialRefinementsCount = 5;
const Algorithm algorithm = Multiscale;
const SolvedExample solvedExample = Benchmark;
const EulerLimiterType limiter_type = VertexBased;

bool HermesView = true;
bool VTKView = false;

double time_step_length;
double time_interval_length;
Hermes::Mixins::Loggable logger(true);

double diffusivity = 1e-15;
double s = -1;

int main(int argc, char* argv[])
{
  if(argc > 1)
    initialRefinementsCount = atoi(argv[1]);
  if(argc > 2)
    diffusivity = (double)atof(argv[2]);
    
  double sigma = std::pow(2., (double)(initialRefinementsCount)) * (s == -1 ? 10.0 : (s == 1 ? 10. : 0.));

  // test();
  Hermes::Mixins::Loggable logger(true);
  Hermes::Mixins::Loggable logger_details(true);
  std::stringstream ss;
  ss << "logfile_" << initialRefinementsCount << "_eps=" << diffusivity << "_s=" << s << ".h2d";
  logger.set_logFile_name(ss.str());
  
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  switch(solvedExample)
  {
  case AdvectedCube:
    time_step_length = 0.01;
    time_interval_length = 1.;
    break;
  case SolidBodyRotation:
    time_step_length = 0.01;
    time_interval_length = 2 * M_PI;
    break;
  case CircularConvection:
    time_step_length = 1e0;
    time_interval_length = 1e4;
    break;
  case MovingPeak:
    time_step_length = 1e-3;
    time_interval_length = (2. * M_PI) + (time_step_length / 10.);
    break;
  case Benchmark:
    time_step_length = 1e0;
    time_interval_length = 1000. + time_step_length / 10.;
    break;
  }

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    mloader.load("domain_rotation.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case AdvectedCube:
    mloader.load("larger_domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case CircularConvection:
    mloader.load("domain_circular_convection.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case MovingPeak:
    mloader.load("domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case Benchmark:
    mloader.load("domain_benchmark.xml", mesh);
    mesh->refine_all_elements(2);
    mesh->refine_all_elements(2);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  }
  
  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition_der = NULL;
  ExactSolutionScalar<double>* initial_solution = NULL;
  ExactSolutionScalar<double>* exact_sln = NULL;
  switch(solvedExample)
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
    initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    initial_condition_der = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    previous_initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    exact_sln = new ExactSolutionMovingPeak(mesh, diffusivity, (1./2.) * M_PI);
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
  MeshFunctionSharedPtr<double> solution(new Solution<double>);
  MeshFunctionSharedPtr<double> previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double> previous_mean_values(initial_condition);
  MeshFunctionSharedPtr<double> previous_derivatives(initial_condition_der);
  MeshFunctionSharedPtr<double> exact_solution(exact_sln);
  MeshFunctionSharedPtr<double> initial_sln(initial_solution);

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(0, 0, 600, 350));
  ScalarView exact_view("Exact solution", new WinGeom(610, 0, 600, 350));
  //exact_view.show(exact_solution);
  
  // Exact solver solution
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  logger_details.info("Exact solver");
  solve_exact(solvedExample, space, diffusivity, s, sigma, exact_solution, initial_sln, time_step_length, logger, logger_details);
  logger_details.info("\n");
  
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();
  if(algorithm == Multiscale)
  {
    logger.info("Multiscale solver");
    
    multiscale_decomposition(mesh, solvedExample, polynomialDegree, previous_mean_values, previous_derivatives, diffusivity, s, sigma, time_step_length,
    time_interval_length, solution, exact_solution, &solution_view, &exact_view, logger, logger_details);
    
    cpu_time.tick();
    logger.info("Multiscale total: %f", cpu_time.last());
    logger.info("\n");
  }
  
  if(algorithm == pMultigrid)
  {
    int steps[4] = { 2, 3, 5, 10};
    for(int si = 0; si < 4; si++)
    {
      cpu_time.tick();
      logger.info("p-Multigrid solver - %i steps", steps[si]);

      MeshFunctionSharedPtr<double> previous_solution_local(new ZeroSolution<double>(mesh));
  
      p_multigrid(mesh, solvedExample, polynomialDegree, previous_solution_local, diffusivity, time_step_length, time_interval_length, 
        solution, exact_solution, &solution_view, &exact_view, s, sigma, logger, logger_details, steps[si]);
        
      cpu_time.tick();
      logger.info("p-Multigrid solver - %i steps total: %f", steps[si], cpu_time.last());
    }
  }
  
  if(false)
  {
    cpu_time.tick();
    logger.info("only Smoothing solver");

    smoothing(mesh, solvedExample, polynomialDegree, previous_solution, diffusivity, time_step_length, time_interval_length, 
      solution, exact_solution, &solution_view, &exact_view, s, sigma, logger, logger_details, 1);
      
    cpu_time.tick();
    logger.info("only Smoothing total: %s", cpu_time.last_str().c_str());
  }
  
  //View::wait();
  return 0;
}
