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
const unsigned int EVERY_NTH_STEP = 100;

bool SHOCK_CAPTURING = true;
const EulerLimiterType limiter_type = VertexBased;
bool limit_velocities = false;

// Initial polynomial degree.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
int INIT_REF_NUM = 4;
// Initial time step.
double time_step_length = 1e-6;
double TIME_INTERVAL_LENGTH = 20.;
// CFL value.
double CFL_NUMBER = 2.0;
// Kappa.
const double KAPPA = 1.4;

// Set up CFL calculation class.
CFLCalculation CFL(CFL_NUMBER, KAPPA);

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.5;         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.25;       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        

// Mesh filename.
const std::string MESH_FILENAME = "GAMM-channel.mesh";

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(numThreads, 16);
  PostProcessing::VertexBasedLimiter::wider_bounds_on_boundary = true;


  Hermes::Mixins::Loggable logger(true);
  logger.set_logFile_name("computation.log");

#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh.
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> inlet_markers;
  Hermes::vector<std::string> outlet_markers;
  MeshSharedPtr mesh(new Mesh);

  MeshReaderH2D mloader;
  mloader.load(MESH_FILENAME.c_str(), mesh);
  solid_wall_markers.push_back(BDY_SOLID_WALL_BOTTOM);
  solid_wall_markers.push_back(BDY_SOLID_WALL_TOP);
  inlet_markers.push_back(BDY_INLET);
  outlet_markers.push_back(BDY_OUTLET);

  // Perform initial mesh refinements.
  mesh->refine_element_id(1, 2);
  mesh->refine_all_elements(1);
  mesh->refine_towards_boundary(BDY_SOLID_WALL_BOTTOM, 2, false);
  
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);

  SpaceSharedPtr<double> const_space_rho(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_rho_v_x(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_rho_v_y(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_e(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > const_spaces(const_space_rho, const_space_rho_v_x, const_space_rho_v_y, const_space_e);

  int ndof = Space<double>::get_num_dofs(spaces);
  logger.info("Ndof: %d", ndof);
#pragma endregion

#pragma region 2. Prev slns
  // Set initial conditions.
  MeshFunctionSharedPtr<double> sln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> sln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> sln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> sln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);
  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  // Set initial conditions.
  MeshFunctionSharedPtr<double> updated_prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> updated_prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> updated_prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> updated_prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > updated_prev_slns(updated_prev_rho, updated_prev_rho_v_x, updated_prev_rho_v_y, updated_prev_e);
#pragma endregion

  EulerEquationsWeakFormSemiImplicit wf_implicit(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,
    solid_wall_markers, inlet_markers, outlet_markers,
    prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, true);

  EulerEquationsWeakFormExplicit wf_explicit(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,
    solid_wall_markers, inlet_markers, outlet_markers,
    prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, updated_prev_rho, updated_prev_rho_v_x, updated_prev_rho_v_y, updated_prev_e, (P_INIT == 0));

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(slns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(slns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(slns));

  ScalarView density_view("Density", new WinGeom(0, 0, 600, 300));
  ScalarView pressure_view("Pressure", new WinGeom(0, 330, 600, 300));
  ScalarView M_view("Mach number", new WinGeom(630, 0, 600, 300));
#pragma endregion

  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  LinearSolver<double> solver_implicit(&wf_implicit, const_spaces);
  ((IterSolver<double>*)solver_implicit.get_linear_solver())->set_tolerance(1e-6, IterSolver<double>::AbsoluteTolerance);
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_UMFPACK);
  LinearSolver<double> solver_explicit(&wf_explicit, spaces);
  solver_explicit.set_jacobian_constant();
  wf_explicit.set_current_time_step(time_step_length);
  wf_implicit.set_current_time_step(time_step_length);

  DiscreteProblemDGAssembler<double>::dg_order = 10;

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    logger.info("---- Time step %d, time %3.5f.", iteration, t);

    if(iteration)
    {
      // Solve.
      // 0 - store the sln vector.
      double* previous_sln_vector = new double[ndof];
      OGProjection<double>::project_global(spaces, prev_slns, previous_sln_vector);

      // 1 - solve implicit.
      solver_implicit.solve();
      double* mean_values = new double[ndof];
      Solution<double>::vector_to_solutions(solver_implicit.get_sln_vector(), const_spaces, updated_prev_slns);
      OGProjection<double>::project_global(spaces, updated_prev_slns, mean_values);

      // 2 - Update the mean values.
      Space<double>::assign_dofs(spaces);
      for(int component = 0; component < 4; component++)
      {
        Element* e;
        for_all_active_elements(e, mesh)
        {
          AsmList<double> al_fine;
          spaces[component]->get_element_assembly_list(e, &al_fine);
          for(unsigned int shape_i = 0; shape_i < al_fine.cnt; shape_i++)
          {
            int order = spaces[component]->get_shapeset()->get_order(al_fine.idx[shape_i], e->get_mode());
            if(order == 0)
            {
              int dof = al_fine.dof[shape_i];
              previous_sln_vector[dof] = mean_values[dof];
            }
          }
        }
      }

      // Solve explicit.
      Solution<double>::vector_to_solutions(previous_sln_vector, spaces, updated_prev_slns);

      // Clean up.
      delete [] previous_sln_vector;
      delete [] mean_values;
    }

    solver_explicit.solve();

#pragma region *. Get the solution with optional shock capturing.
    if(!SHOCK_CAPTURING)
      Solution<double>::vector_to_solutions(solver_explicit.get_sln_vector(), spaces, slns);
    else
    {
      PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, spaces, solver_explicit.get_sln_vector(), 1);
      limiter->get_solutions(slns);
      if(limit_velocities)
        limitVelocityAndEnergy(spaces, limiter, slns);
      delete limiter;
    }
#pragma endregion

     // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        Mach_number->reinit();
        pressure->reinit();
        density_view.show(slns[0]);
        pressure_view.show(pressure);
        M_view.show(Mach_number);
      }
 
#pragma region 4.1. Visualization
    if(((iteration - 1) % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
    {
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure->reinit();
        Mach_number->reinit();
        Linearizer lin;
        char filename[40];
        sprintf(filename, "Pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(pressure, filename, "Pressure", false);
        sprintf(filename, "Mach-%i.vtk", iteration - 1);
        lin.save_solution_vtk(Mach_number, filename, "Velocity", false);
        sprintf(filename, "Rho-%i.vtk", iteration - 1);
        lin.save_solution_vtk(prev_rho, filename, "Rho", false);
      }
    }
#pragma endregion
    
    // Calculate time step according to CFL condition.
    CFL.calculate(solver_explicit.get_sln_vector(), spaces, time_step_length);
    wf_explicit.set_current_time_step(time_step_length);
    wf_implicit.set_current_time_step(time_step_length);

    prev_rho->copy(sln_rho);
    prev_rho_v_x->copy(sln_rho_v_x);
    prev_rho_v_y->copy(sln_rho_v_y);
    prev_e->copy(sln_e);

    updated_prev_rho->copy(sln_rho);
    updated_prev_rho_v_x->copy(sln_rho_v_x);
    updated_prev_rho_v_y->copy(sln_rho_v_y);
    updated_prev_e->copy(sln_e);

    iteration++;
  }
#pragma endregion

  // Done.
  View::wait();
  return 0;
}
