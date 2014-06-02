#include "hermes2d.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;
using namespace Hermes::Algebra;

// for stationary examples - returns nothing, just stores the exact solver solution for error calculation
void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step, int poly_degree, int init_ref_num);

// returns log record
std::string multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values,
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, int time_step_count, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, double cfl);

// returns log record
std::string p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_sln,
                              double diffusivity, double time_step_length, int time_step_count, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, int steps, double cfl);

// returns nothing - just stores the exact solver solution for the error calculation
void exact_solver_timedep(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, double diffusivity, double s, double sigma, double time_step_length, int time_step_count,
  MeshFunctionSharedPtr<double> previous_solution, MeshFunctionSharedPtr<double> exact_solution, ScalarView* exact_view, double cfl);

std::string multiscale_decomposition_timedep(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values,
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, int time_step_count, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, double cfl);

// Utilities.
bool add_inlet(SolvedExample solvedExample);
/* double end_time(SolvedExample solvedExample); */
bool is_timedep(SolvedExample solvedExample);
