#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This is the IAEA EIR-2 benchmark problem. Note the way of handling different material
//  parameters. This is an alternative to how this is done in tutorial examples 07 and 12
//  and in example "iron-water".
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources.
//
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh).
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary").
//       Zero Neumann for the left and bottom edges ("reflection boundary").
//
//  The following parameters can be changed:

const int P_INIT = 1;                    // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;              // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.6;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
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
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters
double LH = 96;                          // Total horizontal length.
double LH0 = 18;                         // First horizontal length.
double LH1 = 48;                         // Second horizontal length.
double LH2 = 78;                         // Third horizontal length.
double LV = 96;                          // Total vertical length.
double LV0 = 18;                         // First vertical length.
double LV1 = 48;                         // Second vertical length.
double LV2 = 78;                         // Third vertical length.
double SIGMA_T_1 = 0.60;                 // Total cross-sections.
double SIGMA_T_2 = 0.48;
double SIGMA_T_3 = 0.70;
double SIGMA_T_4 = 0.85;
double SIGMA_T_5 = 0.90;
double SIGMA_S_1 = 0.53;                 // Scattering cross sections.
double SIGMA_S_2 = 0.20;
double SIGMA_S_3 = 0.66;
double SIGMA_S_4 = 0.50;
double SIGMA_S_5 = 0.89;
double Q_EXT_1 = 1;                      // Nonzero sources in domains 1 and 3 only,
double Q_EXT_3 = 1;                      // sources in other domains are zero.

// Additional constants
double D_1 = 1/(3.*SIGMA_T_1);           // Diffusion coefficients.
double D_2 = 1/(3.*SIGMA_T_2);
double D_3 = 1/(3.*SIGMA_T_3);
double D_4 = 1/(3.*SIGMA_T_4);
double D_5 = 1/(3.*SIGMA_T_5);
double SIGMA_A_1 = SIGMA_T_1 - SIGMA_S_1;    // Absorption coefficients.
double SIGMA_A_2 = SIGMA_T_2 - SIGMA_S_2;
double SIGMA_A_3 = SIGMA_T_3 - SIGMA_S_3;
double SIGMA_A_4 = SIGMA_T_4 - SIGMA_S_4;
double SIGMA_A_5 = SIGMA_T_5 - SIGMA_S_5;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form_1, bilinear_form_ord, H2D_SYM, 1);
  wf.add_matrix_form(bilinear_form_2, bilinear_form_ord, H2D_SYM, 2);
  wf.add_matrix_form(bilinear_form_3, bilinear_form_ord, H2D_SYM, 3);
  wf.add_matrix_form(bilinear_form_4, bilinear_form_ord, H2D_SYM, 4);
  wf.add_matrix_form(bilinear_form_5, bilinear_form_ord, H2D_SYM, 5);
  wf.add_vector_form(linear_form_1, linear_form_ord, 1);
  wf.add_vector_form(linear_form_3, linear_form_ord, 3);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY);

  // Geometry and position of visualization windows.
  WinGeom* sln_win_geom = new WinGeom{0, 0, 440, 350};
  WinGeom* mesh_win_geom = new WinGeom{450, 0, 400, 350};

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  bool verbose = true;     // Prinf info during adaptivity.
  solve_linear_adapt(&space, &wf, sln, SOLVER_UMFPACK, ref_sln, H2D_H1_NORM, 
                     &selector, &apt, sln_win_geom, mesh_win_geom, verbose);

  // Show the final result.
  ScalarView sview("Final solution", sln_win_geom);
  OrderView  oview("Final mesh", mesh_win_geom);
  sview.show_mesh(false);
  sview.show(sln);
  oview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
