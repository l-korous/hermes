#include "definitions.h"

// Initial polynomial degree of all elements.
const int P_INIT = 3;
// Number of initial mesh refinements.
const int INIT_REF_NUM = 3;

// Problem parameters.
// Relative permittivity.
const double epsr = 1.0;
// Permittivity of vacuum F/m.
const double eps0 = 8.85418782e-12;
const double eps = epsr * eps0;
// Relative permeablity.
// Frequency MHz.
double n = 1;
// Angular velocity.
double omega = M_PI * n;
double tau = 1e-3;
double current_time = 0.;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  ScalarView viewEtr("E_t real", new WinGeom(0, 0, 400, 200));
  ScalarView viewEti("E_t imag", new WinGeom(0, 230, 400, 200));
  ScalarView viewEr("E real", new WinGeom(430, 0, 400, 200));
  ScalarView viewEi("E imag", new WinGeom(430, 230, 400, 200));
  ScalarView viewPr("P real", new WinGeom(860, 0, 400, 200));
  ScalarView viewPi("P imag", new WinGeom(860, 230, 400, 200));
  viewEtr.fix_scale_width();
  viewEti.fix_scale_width();
  viewEr.fix_scale_width();
  viewEi.fix_scale_width();
  viewPr.fix_scale_width();
  viewPi.fix_scale_width();

  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initial conditions.
  MeshFunctionSharedPtr<double> et_r_sln(new CustomSolution(mesh, 0, 0., eps, omega));
  MeshFunctionSharedPtr<double> et_i_sln(new CustomSolution(mesh, 1, 0., eps, omega));
  MeshFunctionSharedPtr<double> e_r_sln(new CustomSolution(mesh, 2, 0., eps, omega));
  MeshFunctionSharedPtr<double> e_i_sln(new CustomSolution(mesh, 3, 0., eps, omega));
  MeshFunctionSharedPtr<double> p_r_sln(new CustomSolution(mesh, 4, 0., eps, omega));
  MeshFunctionSharedPtr<double> p_i_sln(new CustomSolution(mesh, 5, 0., eps, omega));
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(et_r_sln, et_i_sln, e_r_sln, e_i_sln, p_r_sln, p_i_sln);

  // Initialize boundary conditions.
  DefaultEssentialBCNonConst<double> bc_et_r(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), et_r_sln);
  DefaultEssentialBCNonConst<double> bc_et_i(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), et_i_sln);
  DefaultEssentialBCNonConst<double> bc_e_r(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), e_r_sln);
  DefaultEssentialBCNonConst<double> bc_e_i(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), e_i_sln);
  DefaultEssentialBCNonConst<double> bc_p_r(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), p_r_sln);
  DefaultEssentialBCNonConst<double> bc_p_i(Hermes::vector<std::string>("Bdy_left", "Bdy_impedance"), p_i_sln);
  
  EssentialBCs<double> bcs_et_r(&bc_et_r);
  EssentialBCs<double> bcs_et_i(&bc_et_i);
  EssentialBCs<double> bcs_e_r(&bc_e_r);
  EssentialBCs<double> bcs_e_i(&bc_e_i);
  EssentialBCs<double> bcs_p_r(&bc_p_r);
  EssentialBCs<double> bcs_p_i(&bc_p_r);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> et_r_space(new H1Space<double>(mesh, &bcs_et_r, P_INIT));
  SpaceSharedPtr<double> et_i_space(new H1Space<double>(mesh, &bcs_et_i, P_INIT));
  
  SpaceSharedPtr<double> e_r_space(new H1Space<double>(mesh, &bcs_e_r, P_INIT));
  SpaceSharedPtr<double> e_i_space(new H1Space<double>(mesh, &bcs_e_i, P_INIT));

  SpaceSharedPtr<double> p_r_space(new H1Space<double>(mesh, &bcs_p_r, P_INIT));
  SpaceSharedPtr<double> p_i_space(new H1Space<double>(mesh, &bcs_p_i, P_INIT));
  Hermes::vector<SpaceSharedPtr<double> > spaces(et_r_space, et_i_space, e_r_space, e_i_space, p_r_space, p_i_space);
  int ndof = Space<double>::get_num_dofs(spaces);

  // Initialize the weak formulation.
  WeakFormMD wf(eps, omega, tau, "Bdy_left", "Bdy_impedance");
  wf.set_ext(slns);

  Hermes::Hermes2D::LinearSolver<double> solver(&wf, spaces);
  solver.output_matrix();
  solver.output_rhs();

  for (int i = 0; i < 100; i++)
  {
    std::cout << "Time step: " << i << ", time; " << current_time << std::endl;
    try
    {
      solver.solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    };

    // Translate the resulting coefficient vector into Solutions.
    Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, slns);

    // Visualize the solution.
    viewEtr.show(slns[0]);
    viewEti.show(slns[1]);
    viewEr.show(slns[2]);
    viewEi.show(slns[3]);
    viewPr.show(slns[4]);
    viewPi.show(slns[5]);

    current_time += tau;
  }

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
