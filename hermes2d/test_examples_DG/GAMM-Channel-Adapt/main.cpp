#include "hermes2d.h"
#include "../euler_util.h"

double CFL_NUMBER = 20.0;

#define REFERENCE_SOLUTION_APPROACH

#pragma region Visulization + Utilities
double LAST_UNREF = 0.;
// Initial time step.
double time_step_length = 1e-6;
double TIME_INTERVAL_LENGTH = 20.;

// Kappa.
const double KAPPA = 1.4;

// Set up CFL calculation class.
CFLCalculation CFL(CFL_NUMBER, KAPPA);

// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = true;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 100;

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

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";
#pragma endregion

#pragma region Limiting
// Limiting
bool SHOCK_CAPTURING = true;
const EulerLimiterType limiter_type = VertexBased;
bool limit_velocities = false;
#pragma endregion

// Initial polynomial degree.
const int P_INIT = 1;

// Number of initial uniform mesh refinements.
int INIT_REF_NUM = 1;

// Every UNREF_FREQth time step the mesh is unrefined.
double UNREF_FREQ = 1e-2;

#ifdef REFERENCE_SOLUTION_APPROACH
bool adaptivityErrorStop(int iteration, double time, double error, int ref_ndof)
{
  if(time < 0.1)
    return (error < .005);
  else
    return (error < .0005);
}

// Error calculation.
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(AbsoluteError, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionCumulative<double> stoppingCriterion(.6);
// Initialize refinement selector.
HOnlySelector<double> selector;
#else
#endif

int order_increase = 0;
int max_p = 0;

int main(int argc, char* argv[])
{
  Hermes::Mixins::Loggable logger(true);
  logger.set_logFile_name("computation.log");
  Hermes::Mixins::Loggable::set_static_logFile_name("computation.log");
  PostProcessing::VertexBasedLimiter::wider_bounds_on_boundary = true;

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
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements(0, true);

  mesh->refine_towards_vertex(1, 2);
  mesh->refine_towards_vertex(2, 2);
  mesh->refine_all_elements(0, true);

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
  MeshFunctionSharedPtr<double> sln_rho(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_e(new Solution<double>);
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);

  MeshFunctionSharedPtr<double> rsln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> rsln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

  // For updating of mean values.
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
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(rslns));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 400, 200));
  ScalarView Mach_number_view("Mach number", new WinGeom(400, 0, 400, 200));
  ScalarView eview_element("Error - element", new WinGeom(0, 200, 400, 200));
  ScalarView ref_element("Refined - element", new WinGeom(400, 200, 400, 200));
  ScalarView eview_edge("Error - max-edge", new WinGeom(0, 400, 400, 200));
  ScalarView ref_edge("Refined - max-edge", new WinGeom(400, 400, 400, 200));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

#pragma region Solver setup
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  LinearSolver<double> solver_implicit(&wf_implicit, const_spaces);
  ((IterSolver<double>*)solver_implicit.get_linear_solver())->set_tolerance(1e-6, IterSolver<double>::AbsoluteTolerance);
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_UMFPACK);
  LinearSolver<double> solver_explicit(&wf_explicit, spaces);
  solver_explicit.set_jacobian_constant();
  wf_explicit.set_current_time_step(time_step_length);
  wf_implicit.set_current_time_step(time_step_length);
  DiscreteProblemDGAssembler<double>::dg_order = 20;
#pragma endregion

#ifdef REFERENCE_SOLUTION_APPROACH
#pragma region Adaptivity setup
  Hermes::vector<RefinementSelectors::Selector<double> *> selectors(&selector, &selector, &selector, &selector);
  Adapt<double> adaptivity(space_rho, &errorCalculator, &stoppingCriterion);
#pragma endregion
#endif

#pragma region 3.1 Set up reference mesh and spaces.
  Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
  SpaceSharedPtr<double> ref_space_rho;
  SpaceSharedPtr<double> ref_space_rho_v_x;
  SpaceSharedPtr<double> ref_space_rho_v_y;
  SpaceSharedPtr<double> ref_space_e;
  Hermes::vector<SpaceSharedPtr<double>  > ref_spaces;

  SpaceSharedPtr<double> const_ref_space_rho;
  SpaceSharedPtr<double> const_ref_space_rho_v_x;
  SpaceSharedPtr<double> const_ref_space_rho_v_y;
  SpaceSharedPtr<double> const_ref_space_e;
  Hermes::vector<SpaceSharedPtr<double> > const_ref_spaces;
#pragma endregion

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    logger.info("Time step %d, time %3.5f, time step %3.5f.", iteration, t, time_step_length);

    if(iteration == 100)
      CFL_NUMBER *= 1.5;

    if(iteration == 200)
      CFL_NUMBER *= 1.5;

    if(iteration == 300)
      CFL_NUMBER *= 1.25;

#pragma region 4.1. Periodic global derefinements.
    if (iteration > 1 && t > LAST_UNREF + UNREF_FREQ)
    {
      LAST_UNREF = t;
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      if(order_increase > 0)
      {
        space_rho->unrefine_all_mesh_elements(true);
        space_rho->adjust_element_order(-1, P_INIT, 0, 0);
        space_rho_v_x->adjust_element_order(-1, P_INIT, 0, 0);
        space_rho_v_y->adjust_element_order(-1, P_INIT, 0, 0);
        space_e->adjust_element_order(-1, P_INIT, 0, 0);
      }
      else
      {
        mesh->unrefine_all_elements();
        space_rho->set_uniform_order(P_INIT);
        space_rho_v_x->set_uniform_order(P_INIT);
        space_rho_v_y->set_uniform_order(P_INIT);
        space_e->set_uniform_order(P_INIT);
      }
      Space<double>::assign_dofs(spaces);

      const_space_rho->set_uniform_order(0);
      const_space_rho_v_x->set_uniform_order(0);
      const_space_rho_v_y->set_uniform_order(0);
      const_space_e->set_uniform_order(0);

      Space<double>::assign_dofs(const_spaces);
    }
#pragma endregion

#pragma region 4.2. Adaptivity loop.
    int as = 1;
    do
    {
#pragma region 7.1 Create reference mesh and spaces.
#ifdef REFERENCE_SOLUTION_APPROACH
      MeshSharedPtr ref_mesh = refMeshCreatorFlow.create_ref_mesh();
#else
      MeshSharedPtr ref_mesh(new Mesh);
      ref_mesh->copy(mesh);
#endif
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, order_increase);

      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRho(const_space_rho, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRhoVx(const_space_rho_v_x, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRhoVy(const_space_rho_v_y, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorE(const_space_e, ref_mesh, 0);

      ref_space_rho = refSpaceCreatorRho.create_ref_space();
      ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      ref_space_e = refSpaceCreatorE.create_ref_space();
      ref_spaces = Hermes::vector<SpaceSharedPtr<double> >(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);
      solver_explicit.set_spaces(ref_spaces);

      const_ref_space_rho = const_refSpaceCreatorRho.create_ref_space();
      const_ref_space_rho_v_x = const_refSpaceCreatorRhoVx.create_ref_space();
      const_ref_space_rho_v_y = const_refSpaceCreatorRhoVy.create_ref_space();
      const_ref_space_e = const_refSpaceCreatorE.create_ref_space();
      const_ref_spaces = Hermes::vector<SpaceSharedPtr<double> >(const_ref_space_rho, const_ref_space_rho_v_x, const_ref_space_rho_v_y, const_ref_space_e);
      solver_implicit.set_spaces(const_ref_spaces);
#pragma endregion

      int ref_ndof = Space<double>::get_num_dofs(ref_spaces);
      logger.info("\tAdaptivity step %d, NDOFs: %i.", as, ref_ndof);

      if(iteration)
      {
        // Solve.
        // 0 - store the sln vector.
        double* previous_sln_vector = new double[ref_ndof];
        OGProjection<double>::project_global(ref_spaces, prev_slns, previous_sln_vector);

        // 1 - solve implicit.
        solver_implicit.solve();
        double* mean_values = new double[ref_ndof];
        Solution<double>::vector_to_solutions(solver_implicit.get_sln_vector(), const_ref_spaces, updated_prev_slns);
        OGProjection<double>::project_global(ref_spaces, updated_prev_slns, mean_values);

        // 2 - Update the mean values.
        Space<double>::assign_dofs(ref_spaces);
        for(int component = 0; component < 4; component++)
        {
          Element* e;
          for_all_active_elements(e, ref_mesh)
          {
            AsmList<double> al_fine;
            ref_spaces[component]->get_element_assembly_list(e, &al_fine);
            for(unsigned int shape_i = 0; shape_i < al_fine.cnt; shape_i++)
            {
              int order = ref_spaces[component]->get_shapeset()->get_order(al_fine.idx[shape_i], e->get_mode());
              if(order == 0)
              {
                int dof = al_fine.dof[shape_i];
                previous_sln_vector[dof] = mean_values[dof];
              }
            }
          }
        }

        // Solve explicit.
        Solution<double>::vector_to_solutions(previous_sln_vector, ref_spaces, updated_prev_slns);

        // Clean up.
        delete [] previous_sln_vector;
        delete [] mean_values;
      }

      solver_explicit.solve();

#pragma region *. Get the solution with optional shock capturing.
      if(!SHOCK_CAPTURING)
        Solution<double>::vector_to_solutions(solver_explicit.get_sln_vector(), ref_spaces, rslns);
      else
      {
        PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, ref_spaces, solver_explicit.get_sln_vector(), 1);
        limiter->get_solutions(rslns);
        if(limit_velocities)
          limitVelocityAndEnergy(ref_spaces, limiter, rslns);
        delete limiter;
      }
#pragma endregion

      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        pressure->reinit();
        Mach_number->reinit();
        pressure_view.show(pressure);
        Mach_number_view.show(Mach_number);
      }
        
#pragma region 7.3.1 Visualization
      if((iteration % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
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
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(Mach_number, filename, "Velocity", false);
          Orderizer ord;
          sprintf(filename, "Mesh-%i.vtk", iteration - 1);
          ord.save_mesh_vtk(ref_spaces[0], filename);
        }
      }
#pragma endregion

      CFL.calculate(solver_explicit.get_sln_vector(), ref_spaces, time_step_length);

#ifdef REFERENCE_SOLUTION_APPROACH
#pragma region 7.4 Project to coarse mesh -> error estimation -> space adaptivity
      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("\t\tProjecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(spaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(sln_rho, rsln_rho);
      double density_error = errorCalculator.get_error_squared(0) * 100;
      if(iteration)
      {
        std::cout << errorCalculator.get_element_error_squared(0, 38) << std::endl;
        std::cout << errorCalculator.get_element_error_squared(0, 21) << std::endl;
      }
      if(HERMES_VISUALIZATION)
        eview_element.show(errorCalculator.get_errorMeshFunction(0));

      // Report results.
      logger.info("\t\tDensity error: %g%%.", density_error);

      // If err_est too large, adapt the mesh.
      if (adaptivityErrorStop(iteration, t, density_error, ref_ndof))
        break;
      else
      {
        Hermes::Mixins::Loggable::Static::info("\t\tAdapting coarse mesh.");
        adaptivity.adapt(&selector);
        space_rho_v_x->copy(space_rho, mesh);
        space_rho_v_y->copy(space_rho, mesh);
        space_e->copy(space_rho, mesh);
        Space<double>::assign_dofs(spaces);

        const_space_rho->copy(space_rho, mesh);
        const_space_rho_v_x->copy(space_rho, mesh);
        const_space_rho_v_y->copy(space_rho, mesh);
        const_space_e->copy(space_rho, mesh);

        const_space_rho->set_uniform_order(0);
        const_space_rho_v_x->set_uniform_order(0);
        const_space_rho_v_y->set_uniform_order(0);
        const_space_e->set_uniform_order(0);

        Space<double>::assign_dofs(const_spaces);
      }
#pragma endregion

#else
      DensityErrorCalculator eC;
      double element_error_threshold = 1e-4;
      double element_size_threshold = 1e-4;
      double element_length_threshold = 1e-2;
      double edge_error_power = -1.0;
      double edge_error_threshold = 1e-6;

      eC.process(rsln_rho, ref_spaces[0]);

      std::set<int> element_numbers_to_be_refined_0;
      std::set<int> element_numbers_to_be_refined_1;
      std::set<int> element_numbers_to_be_refined_2;

      bool* refinements_element = (bool*)calloc(ref_mesh->get_max_element_id(), sizeof(bool));
      bool* refinements_edge = (bool*)calloc(ref_mesh->get_max_element_id(), sizeof(bool));

      Element* e;
      for_all_active_elements(e, ref_mesh)
      {
        // Adjust errors.
        double element_size = e->get_area();
        double element_diam = e->get_diameter();
        eC.element_errors[e->id] *= element_size;

        for(int edge = 0; edge < e->nvert; edge++)
          eC.edge_errors[e->id][edge] /= std::pow(element_diam, edge_error_power);

        // Take a look at element sizes.
        bool element_too_small = false;
        if(element_size < element_size_threshold)
          element_too_small = true;
        for(int edge = 0; edge < e->nvert; edge++)
        {
          if(std::abs(e->vn[edge]->x - e->vn[(edge+1)%e->nvert]->x) < element_length_threshold)
            element_too_small = true;
          if(std::abs(e->vn[edge]->y - e->vn[(edge+1)%e->nvert]->y) < element_length_threshold)
            element_too_small = true;
        }
        if(element_too_small)
          continue;

        if(eC.element_errors[e->id] > element_error_threshold)
        {
          element_numbers_to_be_refined_0.insert(e->id);
          refinements_element[e->id] = true;
          continue;
        }

        if(e->is_quad())
        {
          bool jump_vertical = false;
          bool jump_horizontal = false;
          if((eC.edge_errors[e->id][0] > edge_error_threshold) || (eC.edge_errors[e->id][2] > edge_error_threshold))
            jump_horizontal = true;
          if((eC.edge_errors[e->id][1] > edge_error_threshold) || (eC.edge_errors[e->id][3] > edge_error_threshold))
            jump_vertical = true;

          if(!(jump_vertical || jump_horizontal))
            continue;
          if(jump_vertical && jump_horizontal)
            element_numbers_to_be_refined_0.insert(e->id);
          else if(jump_horizontal)
            element_numbers_to_be_refined_1.insert(e->id);
          else if(jump_vertical)
            element_numbers_to_be_refined_2.insert(e->id);
          refinements_edge[e->id] = true;
        }
        else
        {
          for(int edge = 0; edge < 3; edge++)
            if(eC.edge_errors[e->id][edge] > edge_error_threshold)
            {
              element_numbers_to_be_refined_0.insert(e->id);
              refinements_edge[e->id] = true;
            }
        }
      }

      double* max_errors_edge = new double[ref_mesh->get_max_element_id()];
      memset(max_errors_edge, 0, sizeof(double) * ref_mesh->get_max_element_id());
      for_all_active_elements(e, ref_mesh)
      {
        for(int j = 0; j < e->nvert; j++)
          if(eC.edge_errors[e->id][j] > max_errors_edge[e->id])
            max_errors_edge[e->id] = eC.edge_errors[e->id][j];
      }

      MeshFunctionSharedPtr<double> error_fn_element(new ExactSolutionConstantArray<double, double>(ref_mesh, eC.element_errors));
      MeshFunctionSharedPtr<double> error_fn_edge(new ExactSolutionConstantArray<double, double>(ref_mesh, max_errors_edge));

      MeshFunctionSharedPtr<double> ref_fn_element(new ExactSolutionConstantArray<double, bool>(ref_mesh, refinements_element));
      MeshFunctionSharedPtr<double> ref_fn_edge(new ExactSolutionConstantArray<double, bool>(ref_mesh, refinements_edge));

      ref_element.show(ref_fn_element);
      ref_edge.show(ref_fn_edge);

      eview_element.show(error_fn_element);
      eview_edge.show(error_fn_edge);

      delete [] max_errors_edge;
      ::free(refinements_element);
      ::free(refinements_edge);

      if(element_numbers_to_be_refined_0.size() == 0 &&
        element_numbers_to_be_refined_1.size() == 0 &&
        element_numbers_to_be_refined_2.size() == 0)
        break;
      else
      {
        as++;

        for(std::set<int>::iterator it = element_numbers_to_be_refined_0.begin(); it != element_numbers_to_be_refined_0.end(); it++)
          mesh->refine_element_id(*it, 0);
        for(std::set<int>::iterator it = element_numbers_to_be_refined_1.begin(); it != element_numbers_to_be_refined_1.end(); it++)
          mesh->refine_element_id(*it, 1);
        for(std::set<int>::iterator it = element_numbers_to_be_refined_2.begin(); it != element_numbers_to_be_refined_2.end(); it++)
          mesh->refine_element_id(*it, 2);

        space_rho->set_uniform_order(P_INIT);
        space_rho_v_x->set_uniform_order(P_INIT);
        space_rho_v_y->set_uniform_order(P_INIT);
        space_e->set_uniform_order(P_INIT);
        Space<double>::assign_dofs(spaces);

        const_space_rho->set_uniform_order(0);
        const_space_rho_v_x->set_uniform_order(0);
        const_space_rho_v_y->set_uniform_order(0);
        const_space_e->set_uniform_order(0);
        Space<double>::assign_dofs(const_spaces);
      }
#endif
    }
    while (true);

    wf_explicit.set_current_time_step(time_step_length);
    wf_implicit.set_current_time_step(time_step_length);

    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);

    updated_prev_rho->copy(rsln_rho);
    updated_prev_rho_v_x->copy(rsln_rho_v_x);
    updated_prev_rho_v_y->copy(rsln_rho_v_y);
    updated_prev_e->copy(rsln_e);

    iteration++;
  }
#pragma endregion

  // Done.
  View::wait();
  return 0;
}
