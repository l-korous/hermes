#define HERMES_REPORT_ALL
#include "definitions.h"
#include "bddcml_wrapper.hpp"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int P_INIT = 1;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  Hermes::Hermes2D::H1Space<double>* space = NULL;
  Hermes::Hermes2D::Mesh* mesh = new Mesh();

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("a", FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf(new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));
  
  // This is in a block to test that the instances mesh and space can be deleted after being copied with no harm.
  {
    // Set the number of threads used in Hermes.
    Hermes::HermesCommonApi.set_param_value(Hermes::exceptionsPrintCallstack, 0);
    Hermes::Hermes2D::Hermes2DApi.set_param_value(Hermes::Hermes2D::numThreads, 8);

    // Load the mesh.
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    try
    {
      mloader.load("domain.xml", mesh);

      // Create an H1 space with default shapeset.
      space = new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, 1);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
  }

  Mesh* new_mesh = new Mesh();
  H1Space<double>* new_space = new H1Space<double>();
  new_space->copy(space, new_mesh);

  Views::BaseView<double> b;
  b.show(new_space);
  b.wait_for_close();

  delete space;
  delete mesh;

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  // Initialize linear solver.
  Hermes::Hermes2D::DiscreteProblemLinear<double> dp(&wf, new_space);

  /// Jacobian.
  SparseMatrix<double>* jacobian = Algebra::create_matrix<double>();

  /// Residual.
  Vector<double>* residual = Algebra::create_vector<double>();

  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);

  // Solve the linear problem.
  try
  {
    dp.assemble(jacobian, residual);

    FILE* f = fopen("matrix", "w");
    jacobian->dump(f, "A");
    fclose(f);

    FILE* fa = fopen("rhs", "w");
    residual->dump(fa, "b");
    fclose(fa);

    // No a tady se zavola BDDC.
    // magic();
    matrix_solver->solve();
    Hermes2D::Solution<double>::vector_to_solution(matrix_solver->get_sln_vector(), new_space, &sln);

    // VTK output.
    if(VTK_VISUALIZATION)
    {
      // Output solution in VTK format.
      Hermes::Hermes2D::Views::Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D, 1, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

      // Output mesh and element orders in VTK format.
      Hermes::Hermes2D::Views::Orderizer ord;
      ord.save_mesh_vtk(new_space, "mesh.vtk");
      ord.save_orders_vtk(new_space, "ord.vtk");
    }

    // Visualize the solution.
    Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

    if(HERMES_VISUALIZATION)
    {
      viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_LOW);
      viewS.wait_for_close();
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  return 0;
}
