#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;


// Initial polynomial degree of mesh elements.
const int P_INIT = 1;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;

// Problem parameters.
// Slope of the layer.
double slope = 70.;

int criterion(Element* e)
{
  if (e->vn[0]->x < 0.3 && e->vn[0]->x > 0.15 && e->vn[0]->y > 0.2 && e->vn[0]->y < 0.7)
    return 0;
  else if (e->vn[1]->x < 0.3 && e->vn[1]->x > 0.15 && e->vn[1]->y > 0.2 && e->vn[1]->y < 0.7)
    return 0;
  else if (e->vn[2]->x < 0.3 && e->vn[2]->x > 0.15 && e->vn[2]->y > 0.2 && e->vn[2]->y < 0.7)
    return 0;
  else
    return -1;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  std::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);
  mloader.load("gmsh.xml", meshes);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  
  //mesh->refine_by_criterion(&criterion, 6);

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, slope));

  // Define custom function f.
  CustomFunction f(slope);

  // Initialize the weak formulation.
  Hermes::Hermes1DFunction<double> lambda(1.0);
  WeakFormSharedPtr<double> wf(new DefaultWeakFormPoissonLinear<double>(HERMES_ANY, &f));
  
  // Initialize boundary conditions
  //DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  DefaultEssentialBCNonConst<double> bc_essential({ "0", "1", "2", "3" }, exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  
  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(H2D_HP_ANISO);

  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

  AdaptStoppingCriterionSingleElement<double> criterion(0.85);
  
  AdaptSolverCriterionFixed global_criterion(500);

  AdaptSolver<double, LinearSolver<double> >::view_size = 600;
  AdaptSolver<double, LinearSolver<double> >::scalar_views_switch = false;
  AdaptSolver<double, LinearSolver<double> >::scalar_vtk_switch = true;
  AdaptSolver<double, LinearSolver<double> >::order_views_switch = false;
  AdaptSolver<double, LinearSolver<double> >::order_vtk_switch = true;
  AdaptSolver<double, LinearSolver<double> >::adapt_views_switch = false;
  AdaptSolver<double, LinearSolver<double> >::adapt_vtk_switch = true;
  AdaptSolver<double, LinearSolver<double> >::base_views_switch = false;
  AdaptSolver<double, LinearSolver<double> > adaptSolver(space, wf, &errorCalculator, &criterion, &selector, &global_criterion);

  adaptSolver.switch_visualization(true);
  adaptSolver.switch_vtk(false);
  adaptSolver.set_verbose_output(true);
  adaptSolver.solve(hpAdaptivity);
  return 0;
}