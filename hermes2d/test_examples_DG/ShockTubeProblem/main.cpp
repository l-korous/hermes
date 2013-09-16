#include "hermes2d.h"
#include "../euler_util.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// Discontinuous Galerkin method of higher order with no adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Solid walls, inlet, outlet.
//
// IC: Constant state identical to inlet.
//

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = true;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 50;

bool SHOCK_CAPTURING = true;
EulerLimiterType limiter_type = CoarseningJumpIndicatorDensityToAll;
bool limit_velocities = false;

// Initial polynomial degree.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
int INIT_REF_NUM_ISO = 5;
int INIT_REF_NUM_ANISO = 2;
// CFL value.
double CFL_NUMBER = 0.1;
// Initial time step.
double time_step_length = 1E-4;

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT_IN = 1.;         
// Inlet density (dimensionless).   
const double RHO_EXT_IN = 1.0;       
// Inlet x-velocity (dimensionless).
const double V1_EXT_IN = 0.;       
// Inlet y-velocity (dimensionless).
const double V2_EXT_IN = 0.0;

// Exterior pressure (dimensionless).
const double P_EXT_OUT = .1;         
// Inlet density (dimensionless).   
const double RHO_EXT_OUT = .125;       
// Inlet x-velocity (dimensionless).
const double V1_EXT_OUT = 0.;       
// Inlet y-velocity (dimensionless).
const double V2_EXT_OUT = 0.0;    
// Kappa.
const double KAPPA = 1.4;

double TIME_INTERVAL_LENGTH = .231;

// Boundary markers.
std::string BDY_INLET = "Left";
std::string BDY_OUTLET = "Right";
std::string BDY_SOLID_WALL_BOTTOM = "Bottom";
std::string BDY_SOLID_WALL_TOP = "Top";

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

std::string filename = "domain.xml";

void set_params(int argc, char* argv[])
{
  if(argc < 3)
    return;

  // argv[1]: EulerLimiterType
  // argv[2]: 1 - quad, 2 - tri-simple, 3 - tri
  // argv[3] ? alpha
  limiter_type = (EulerLimiterType)(atoi(argv[1]));

  if(atoi(argv[1]) < 2)
  {
    if(atoi(argv[2]) == 1)
    {
      INIT_REF_NUM_ISO = 0;
      INIT_REF_NUM_ANISO = 7;
    }
  }
  else
  {
    if(atoi(argv[2]) == 1)
    {
      INIT_REF_NUM_ISO = 5;
      INIT_REF_NUM_ANISO = 2;
    }
  }

  if(atoi(argv[2]) == 2)
  {
    INIT_REF_NUM_ISO = 6;
    INIT_REF_NUM_ANISO = 0;
  }
  if(atoi(argv[2]) == 3)
  {
    INIT_REF_NUM_ISO = 2;
    INIT_REF_NUM_ANISO = 0;
  }

  if(atoi(argv[2]) == 1)
    filename = "/home/staff/korous/hermes/hermes2d/test_examples/ShockTubeProblem/domain.xml";
  else if (atoi(argv[2]) == 2)
    filename = "/home/staff/korous/hermes/hermes2d/test_examples/ShockTubeProblem/domain-tri-simple.xml";
  else
  {
    filename = "/home/staff/korous/hermes/hermes2d/test_examples/ShockTubeProblem/domain-tri.xml";
    BDY_INLET = "0";
    BDY_OUTLET = "1";
    BDY_SOLID_WALL_BOTTOM = "2";
    BDY_SOLID_WALL_TOP = "3";
  }
  if(argc > 3)
  {
    FeistauerPCoarseningLimiter::alpha = atof(argv[3]);
  }
  if(argc >4 && atoi(argv[4]) == 1)
    limit_velocities = true;
}

int main(int argc, char* argv[])
{
#ifndef _WINDOWS
  HermesCommonApi.set_integral_param_value(numThreads, 1);
#endif
  set_params(argc, argv);
  Hermes::Mixins::Loggable logger(true);
  logger.set_logFile_name("computation.log");

  logger.info("Limiter: %i", limiter_type);
  logger.info("Aniso refs: %i", INIT_REF_NUM_ANISO);
  logger.info("Iso refs: %i", INIT_REF_NUM_ISO);
  if(limiter_type > 1)
    logger.info("Feist alpha: %f", FeistauerPCoarseningLimiter::alpha);
  logger.info("File: %s", filename.c_str());
  logger.info("Limit velocities: %i", limit_velocities);

#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);

  try
  {
    if(argc > 2 && atoi(argv[2]) == 3)
      mloader.load(filename.c_str(), meshes);
    else
      mloader.load(filename.c_str(), mesh);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    std::cout << e.what();
  }

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM_ANISO; i++) 
    mesh->refine_all_elements(2);
  for (int i = 0; i < INIT_REF_NUM_ISO; i++) 
    mesh->refine_all_elements(0);
  
  /* // Test code.
  mesh->refine_element_id(3);
  mesh->refine_element_id(4);

  double values[13] = {1, 2, 3, 0, 0, -4, 4, -4, 4, -5, 5, 5, -5};
  MeshFunctionSharedPtr<double> exact(new ExactSolutionConstantArray<double, double>(mesh, values));
  MeshFunctionSharedPtr<double> target_sln(new Solution<double>());
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  double* target_vec = new double[space->get_num_dofs()];
  OGProjection<double>::project_global(space, exact, target_vec);
  FeistauerPCoarseningLimiter flimiter(space, target_vec);
  FeistauerPCoarseningLimiter::alpha = 1.0;
  flimiter.set_verbose_output(true);
  flimiter.set_type(CoarseningJumpIndicatorAllToThemselves);
  target_sln = flimiter.get_solution();

  ScalarView s;
  s.show(target_sln);
  s.wait_for_close();
  */

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);

  int ndof = Space<double>::get_num_dofs(spaces);
  logger.info("Ndof: %d", ndof);
#pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new InitialSolutionShockTube(mesh, 0));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new InitialSolutionShockTube (mesh, 1));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new InitialSolutionShockTube (mesh, 2));
  MeshFunctionSharedPtr<double> prev_e(new InitialSolutionShockTube (mesh, 3));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET);
  Hermes::vector<std::string> inlet_markers;
  //inlet_markers.push_back(BDY_INLET);
  Hermes::vector<std::string> outlet_markers;
  //outlet_markers.push_back(BDY_OUTLET);

  EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT_IN, V1_EXT_IN, V2_EXT_IN, P_EXT_IN, 
    RHO_EXT_OUT, V1_EXT_OUT, V2_EXT_OUT, P_EXT_OUT,
    solid_wall_markers, inlet_markers, outlet_markers,
    prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, (P_INIT == 0));

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(prev_slns));

  ScalarView density_view("Density", new WinGeom(0, 0, 600, 300));
  ScalarView pressure_view("Pressure", new WinGeom(0, 330, 600, 300));
  ScalarView velocity_view("Velocity", new WinGeom(650, 0, 600, 300));
  ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
  ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

  LinearSolver<double> solver(&wf, spaces);
  DiscreteProblemDGAssembler<double>::dg_order = 6;
  solver.set_jacobian_constant();

  FeistauerJumpDetector limiter(spaces, NULL);
  limiter.set_type(limiter_type);

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    logger.info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_current_time_step(time_step_length);

    // Solve.
    solver.solve();

#pragma region *. Get the solution with optional shock capturing.
    if(!SHOCK_CAPTURING)
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, prev_slns);
    else
    {
      limiter.set_solution_vector(solver.get_sln_vector());
      limiter.get_solutions(prev_slns);
      if(limit_velocities)
        limitVelocityAndEnergy(spaces, &limiter, prev_slns);
    }
#pragma endregion

    // Calculate time step according to CFL condition.
    //CFL.calculate(prev_slns, mesh, time_step_length);

#pragma region 4.1. Visualization
    if(((iteration - 1) % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        Mach_number->reinit();
        pressure->reinit();
        velocity->reinit();
        density_view.show(prev_slns[0]);
        pressure_view.show(pressure);
        velocity_view.show(velocity);
        order_view.show(space_rho);
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure->reinit();
        velocity->reinit();
        Linearizer lin;
        char filename[40];
        sprintf(filename, "Pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(pressure, filename, "Pressure", false);
        sprintf(filename, "Velocity-%i.vtk", iteration - 1);
        lin.save_solution_vtk(velocity, filename, "Velocity", false);
        sprintf(filename, "Rho-%i.vtk", iteration - 1);
        lin.save_solution_vtk(prev_rho, filename, "Rho", false);
      }
    }
#pragma endregion
  }
#pragma endregion

  // Done.
  View::wait();
  return 0;
}
