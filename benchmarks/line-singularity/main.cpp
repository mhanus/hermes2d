#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This is a simple elliptic problem with known exact solution where one
//  can compare isotropic and anisotropic refinements.
//
//  PDE: -Laplace u = f.
//  where f is dictated by the exact solution.
//
//  Exact solution: u(x,y) = cos(K*y)    for x < 0,
//                  u(x,y) = cos(K*y) + pow(x, alpha)   for x > 0   where alpha > 0.
//
//  Domain: square, see the file singpert.mesh.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;              // Number of initial mesh refinements (the original mesh is just one element)
const int P_INIT = 2;                    // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values
                                         // are H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
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
const double ERR_STOP = 0.0001;          // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;            // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Equation parameters.
const double K = M_PI/2;
const double ALPHA = 2.01;

// Exact solution.
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Eessential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return fn(x, y);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_vector_form(linear_form, linear_form_ord);

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
  ExactSolution exact(&mesh, fndd);
  bool verbose = true;     // Prinf info during adaptivity.
  solve_linear_adapt(&space, &wf, sln, SOLVER_UMFPACK, ref_sln, H2D_H1_NORM, 
                     &selector, &apt, sln_win_geom, mesh_win_geom, verbose, &exact);

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
