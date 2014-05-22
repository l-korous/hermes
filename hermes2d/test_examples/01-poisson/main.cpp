#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

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

const bool HERMES_VISUALIZATION = true;   // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;               // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20;        // Fixed temperature on the boundary.

class MyVolumetricIntegralCalculator : public PostProcessing::SurfaceIntegralCalculator<double>
{
public:
  MyVolumetricIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : PostProcessing::SurfaceIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MyVolumetricIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : PostProcessing::SurfaceIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
  {
    for (int i = 0; i < n; i++)
      result[0] += wt[i] * fns[0]->val[i];
  };

  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
    result[0] = Hermes::Ord(21);
  }
};

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_AZTECOO);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);
  mloader.load("sit.msh", meshes);

  MeshSharedPtr eggShell = EggShell::get_egg_shell(mesh, "0", 2);

  MeshView m;
  m.show(mesh);
  m.wait_for_keypress();
  m.show(eggShell);
  m.wait_for_close();
  return 0;
}