#include "hermes2d.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_mean_values, 
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, 
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view);

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                              double diffusivity, double time_step_length, 
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, double s, double sigma);
