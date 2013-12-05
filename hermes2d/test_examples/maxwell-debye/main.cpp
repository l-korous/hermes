#include "definitions.h"

// This example shows how to model harmonic steady state in parallel plate waveguide.
// The complex-valued Helmholtz equation is solved by decomposing it into two equations 
// for the real and imaginary part of the E field. Two typical boundary conditions used in 
// high-frequency problems are demonstrated.
//
// PDE: Helmholtz equation for electric field
//
//    -Laplace E  - omega^2*mu*epsilon*E + j*omega*sigma*mu*E = 0
//
// BC:                     Gamma_1 (perfect)
//                   ----------------------------
//  Gamma_3 (left)   |                           |  Gamma_4 (impedance)
//                   ----------------------------
//                         Gamma_2 (perfect)
//
//     1) Perfect conductor boundary condition Ex = 0 on Gamma_1 and Gamma_2.
//     2) Essential (Dirichlet) boundary condition on Gamma_3
//          Ex(y) = E_0 * cos(y*M_PI/h), where h is height of the waveguide.
//     3) Impedance boundary condition on Gamma_4
//          dE/dn = j*beta*E.
//
// The following parameters can be changed:

// Initial polynomial degree of all elements.
const int P_INIT = 4;
// Number of initial mesh refinements.
const int INIT_REF_NUM = 0;

// Problem parameters.
// Relative permittivity.
const double epsr = 1.0;
// Permittivity of vacuum F/m.
const double eps0 = 8.85418782e-12;
const double eps = epsr * eps0;
// Relative permeablity.
// Frequency MHz.
double n = 10;
// Angular velocity.
double omega = M_PI * n;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc1("Bdy_perfect", 0.0);
  EssentialBCs<double> bcs(&bc1);
  EssentialBCs<double> bcs_im(&bc1);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> e_r_space(new H1Space<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> e_i_space(new H1Space<double>(mesh, &bcs_im, P_INIT));
  Hermes::vector<SpaceSharedPtr<double> > spaces(e_r_space, e_i_space);
  int ndof = Space<double>::get_num_dofs(spaces);

  // Initialize the weak formulation.
  WeakFormHelmholtz wf(eps, omega, "Bdy_left", "Bdy_impedance");

  // Initialize the solutions and solvers.
  MeshFunctionSharedPtr<double> e_r_sln(new Solution<double>), e_i_sln(new Solution<double>);
  Hermes::Hermes2D::LinearSolver<double> solver(&wf, spaces);

  // Views.
  ScalarView viewEr("Er [V/m]", new WinGeom(600, 0, 700, 200));
  viewEr.show_mesh(false);
  ScalarView viewEi("Ei [V/m]", new WinGeom(600, 220, 700, 200));
  ScalarView viewMagnitude("Magnitude of E [V/m]", new WinGeom(600, 440, 700, 200));
  viewMagnitude.show_mesh(false);
  viewMagnitude.show_contours(50., 0.);

  try
  {
    solver.solve();
  }
  catch (Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  };

  // Translate the resulting coefficient vector into Solutions.
  Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(e_r_space, e_i_space),
    Hermes::vector<MeshFunctionSharedPtr<double> >(e_r_sln, e_i_sln));

  // Visualize the solution.
  viewEr.show(e_r_sln);
  viewEr.get_linearizer()->save_solution_tecplot(e_r_sln, "asdf.dat", "asdf");
  // viewEr.save_screenshot("real_part.bmp");

  viewEi.show(e_i_sln);
  // viewEi.save_screenshot("imaginary_part.bmp");

  MeshFunctionSharedPtr<double> magnitude(new MagFilter<double>(Hermes::vector<MeshFunctionSharedPtr<double> >(e_r_sln, e_i_sln)));
  viewMagnitude.show(magnitude);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
