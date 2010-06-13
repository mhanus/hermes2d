#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

//  This test makes sure that the example "iron-water" works correctly.
//  It is only executed if WITH_EXODUSII is set to YES

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;       // Number of initial uniform mesh refinements
const double THRESHOLD = 0.6;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 1e-2;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters
double L = 30;                        // edge of square
double L0 = 0.75*0.5*L;               // end of first water layer
double L1 = 0.5*L;                    // end of second water layer
double L2 = 0.75*L;                   // end of iron layer
double Q_EXT = 1.0;                   // neutron source (nonzero in domain 1 only)
double SIGMA_T_WATER = 3.33;          // total cross-section
double SIGMA_T_IRON = 1.33;
double C_WATER = 0.994;               // scattering ratio
double C_IRON = 0.831;
double D_WATER = 1./(3.*SIGMA_T_WATER);  // diffusion coefficient
double D_IRON = 1./(3.*SIGMA_T_IRON);
double SIGMA_A_WATER = SIGMA_T_WATER - C_WATER*SIGMA_T_WATER;  // absorbing cross-section
double SIGMA_A_IRON = SIGMA_T_IRON - C_IRON*SIGMA_T_IRON;

// Boundary condition types
BCType bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Check input parameters.
  // If true, coarse mesh FE problem is solved in every adaptivity step.
  // If false, projection of the fine mesh solution on the coarse mesh is used. 
  bool SOLVE_ON_COARSE_MESH = false;
  if (argc > 1 && strcmp(argv[1], "-coarse_mesh") == 0)
    SOLVE_ON_COARSE_MESH = true;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  ExodusIIReader mloader;
  if (!mloader.load("iron-water.e", &mesh)) error("ExodusII mesh load failed.");

  // Initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset.
  H1Shapeset shapeset;

  // Create an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions.
  int ndof = assign_dofs(&space);

  // Initialize the weak formulation.
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form_water, bilinear_form_ord, H2D_SYM, 1);
  wf.add_biform(0, 0, bilinear_form_water, bilinear_form_ord, H2D_SYM, 2);
  wf.add_biform(0, 0, bilinear_form_iron, bilinear_form_ord, H2D_SYM, 3);
  wf.add_liform(0, linear_form_source, linear_form_ord, 1);

  // Matrix solver.
  UmfpackSolver solver;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &solver, &space);

  // Adaptivity loop:
  int as = 1; bool done = false;
  Solution sln_coarse, sln_fine;
  do
    {
    info("---- Adaptivity step %d:", as);

    // Assemble and solve the fine mesh problem.
    info("Solving on fine mesh.");
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      info("Solving on coarse mesh.");
      ls.assemble();
      ls.solve(&sln_coarse);
    }
    else {
      info("Projecting fine mesh solution on coarse mesh.");
      ls.project_global(&sln_fine, &sln_coarse);
    }

    // Time measurement.
    cpu_time.tick();

    // Skip exact error calculation time. 
    cpu_time.tick(H2D_SKIP);

    // Calculate error estimate wrt. fine mesh solution.
    info("Calculating error (est).");
    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%",
      space.get_num_dofs(), rs.get_num_dofs(), err_est);

    // Add entries to DOF convergence graphs.
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 2550;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);// ndofs was 2491 at the time this test was created
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
