#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Hermes::Hermes2D::RefinementSelectors;

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -(1/mu)Laplace A + ii*omega*gamma*A - J_ext = 0.
//
//  Domain: Rectangle of height 0.003 and width 0.004. Different
//  materials for the wire, air, and iron (see mesh file domain2.mesh).
//
//  BC: Zero Dirichlet on the top and right edges, zero Neumann
//  elsewhere.
//
//  The following parameters can be changed:
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.95;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
const char* iterative_method = "gmres";           // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers).
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;

typedef std::complex<double> complex;

int main(int argc, char* argv[])
{
  Hermes::Mixins::TimeMeasurable m;
  m.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential("Dirichlet", complex(0.0, 0.0));
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<complex> space(new H1Space<complex>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize the weak formulation.
  CustomWeakForm wf("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON,
    "Wire", MU_0, complex(J_EXT, 0.0), OMEGA);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<complex> sln(new Hermes::Hermes2D::Solution<complex>());
  MeshFunctionSharedPtr<complex> ref_sln(new Hermes::Hermes2D::Solution<complex>());

  // Initialize refinement selector.
  H1ProjBasedSelector<complex> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 600, 350));
  Views::ScalarView sview2("Ref. Solution", new Views::WinGeom(0, 0, 600, 350));
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(610, 0, 520, 350));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  DiscreteProblem<complex> dp(&wf, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<complex> newton(&dp);
    
  Views::MeshView m1, m2;
  Views::OrderView o1, o2;
  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = ref_space_creator.create_ref_space();
    
    newton.set_space(ref_space);

    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.

    // Initial coefficient vector for the Newton's method.
    complex* coeff_vec = new complex[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(complex));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    // For iterative solver.
    if(matrix_solver_type == SOLVER_AZTECOO)
    {
      newton.set_iterative_method(iterative_method);
      newton.set_preconditioner(preconditioner);
    }
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }

    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<complex> ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> real_filter(new RealFilter(sln));
    MeshFunctionSharedPtr<double> rreal_filter(new RealFilter(ref_sln));
    sview.show(real_filter);
    sview2.show(rreal_filter);

    oview.show(space);

    // Calculate element errors and total error estimate.
    DefaultErrorCalculator<complex, HERMES_H1_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
    error_calculator.calculate_errors(sln, ref_sln);

    Adapt<complex> adaptivity(space, &error_calculator);
    adaptivity.set_strategy(AdaptStoppingCriterionLevels, 0.3);
    //adaptivity.set_iterative_improvement(1e-1);

    std::cout << (std::string)"Relative error: " << error_calculator.get_total_error_squared() * 100. << '%' << std::endl;

    // If err_est too large, adapt the mesh->
    if(error_calculator.get_total_error_squared()  * 100. < ERR_STOP)
      done = true;
    else
    {
      std::cout << (std::string)"Adapting..." << std::endl << std::endl;
      adaptivity.adapt(&selector);
    }

    // Clean up.
    delete [] coeff_vec;

    // Increase counter.
    as++;
  }
  while (done == false);

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");

  MeshFunctionSharedPtr<double> real_filter(new RealFilter(ref_sln));
  sview.show(real_filter);

  m.tick();
  std::cout << m.accumulated();

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}