#include "../definitions.h"

// This test makes sure that example 06-bc-newton works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 5;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.

// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 0.0;        // Volume heat sources generated by electric current.
const double ALPHA = 5.0;                  // Heat transfer coefficient.
const double T_EXTERIOR = 50.0;            // Exterior temperature.
const double BDY_A_PARAM = 0.0;
const double BDY_B_PARAM = 0.0;
const double BDY_C_PARAM = 20.0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::Hermes2D::MeshReaderH2D mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL),
    "Copper", new Hermes::Hermes1DFunction<double>(LAMBDA_CU),
    new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC),
    "Outer", ALPHA, T_EXTERIOR);

  // Initialize boundary conditions.
  CustomDirichletCondition bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Left"),
    BDY_A_PARAM, BDY_B_PARAM, BDY_C_PARAM);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  Hermes::Hermes2D::H1Space<double> space(&mesh, &bcs, P_INIT);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  bool success = true;
  for (int p_init = 1; p_init <= 10; p_init++)
  {
    info("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);
    int ndof = space.get_num_dofs();
    info("ndof = %d", ndof);

    // Initialize the FE problem.
    Hermes::Hermes2D::DiscreteProblem<double> dp(&wf, &space);

    // Initial coefficient vector for the Newton's method.
    double* coeff_vec = new double[ndof];
    memset(coeff_vec, 0, ndof*sizeof(double));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    Hermes::Hermes2D::Solution<double> sln;
    Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);
    try{
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
    Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

    double sum = 0;
    for (int i = 0; i < ndof; i++) sum += coeff_vec[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 61.8227) > 1e-1) success = false;
    if (p_init == 2 && fabs(sum - 60.8105) > 1e-1) success = false;
    if (p_init == 3 && fabs(sum - 61.5511) > 1e-1) success = false;
    if (p_init == 4 && fabs(sum - 60.8191) > 1e-1) success = false;
    if (p_init == 5 && fabs(sum - 61.5304) > 1e-1) success = false;
    if (p_init == 6 && fabs(sum - 60.8064) > 1e-1) success = false;
    if (p_init == 7 && fabs(sum - 61.5323) > 1e-1) success = false;
    if (p_init == 8 && fabs(sum - 60.7863) > 1e-1) success = false;
    if (p_init == 9 && fabs(sum - 61.5408) > 1e-1) success = false;
    if (p_init == 10 && fabs(sum - 60.7637) > 1e-1) success = false;

    // Clean up.
    delete [] coeff_vec;
  }

  if (success == 1) {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}
