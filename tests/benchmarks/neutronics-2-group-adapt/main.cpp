#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_neutronics_2g Benchmarks/Neutronics 2-group
 *  \{
 *  \brief This test makes sure that the benchmark "neutronics-2-group-adapt" works correctly.
 *
 *  \section s_params Parameters
 *  - P_INIT=1 for both solution components
 *  - INIT_REF_NUM=1 for both solution components
 *  - THERSHOLD=0.3
 *  - STRATEGY=0
 *  - CAND_LIST=HP_ANISO_P
 *  - MESH_REGULARITY=-1
 *  - ERR_STOP=4
 *  - CONV_EXP=1.0
 *  - NDOF_STOP=40000
 *  - ERROR_WEIGHTS=default values
 *
 *  \section s_res Expected results
 *  - DOFs: 331, 2876   (for the two solution components)
 *  - Iterations: 19    (the last iteration at which ERR_STOP is fulfilled)
 *  - Error:  4.33036%  (H1 norm of error with respect to the exact solution)
 *  - Negatives: 0      (number of negative values) 
 */
 
// INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptivity control:

const int P_INIT[2] =
  {1, 1};                                  // Initial polynomial orders for the individual solution components.
const int INIT_REF_NUM[2] =
  {1, 1};                                  // Initial uniform mesh refinement for the individual solution components.

const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_P; // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 4;                 // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // reference and coarse mesh solution in percent).
const int NDOF_STOP = 40000;               // Adaptivity process stops when the number of degrees of freedom grows over
                                           // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 50;              // Adaptivity process stops when the number of adaptation steps grows over
                                           // this limit.
const int ADAPTIVITY_NORM = 2;             // Specifies the norm used by H1Adapt to calculate the error and norm.
                                           // ADAPTIVITY_NORM = 0 ... H1 norm.
                                           // ADAPTIVITY_NORM = 1 ... norm defined by the diagonal parts of the bilinear form.
                                           // ADAPTIVITY_NORM = 2 ... energy norm defined by the full (non-symmetric) bilinear form.

// Variables used for reporting of results
TimePeriod cpu_time;			                 // Time measurements.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem data:

// Two-group material properties for the 4 macro regions.
const double D[4][2] = {	{1.12, 0.6},
                        	{1.2, 0.5},
                        	{1.35, 0.8},
                        	{1.3, 0.9}	};
const double Sr[4][2] = {	{0.011, 0.13},
                        	{0.09, 0.15},
                        	{0.035, 0.25},
                        	{0.04, 0.35}	};
const double nSf[4][2] ={	{0.0025, 0.15},
                        	{0.0, 0.0},
                        	{0.0011, 0.1},
                        	{0.004, 0.25}	};
const double chi[4][2] ={	{1, 0},
                        	{1, 0},
                        	{1, 0},
                        	{1, 0}	};
const double Ss[4][2][2] = {	{	{ 0.0, 0.0 },
                             		{ 0.05, 0.0 }	},
                             	{	{ 0.0, 0.0 },
                             		{ 0.08, 0.0 }	},
                             	{	{ 0.0, 0.0 },
                             		{ 0.025, 0.0 }	},
                             	{	{ 0.0, 0.0 },
                             		{ 0.014, 0.0 }	}	};

double a = 0., b = 1., c = (a+b)/2.;

inline int get_material(double x, double y)
{
  if (x >= a && x <= c && y >= a && y <= c) return 0;
  if (x >= c && x <= b && y >= a && y <= c) return 1;
  if (x >= c && x <= b && y >= c && y <= b) return 2;
  if (x >= a && x <= c && y >= c && y <= b) return 3;
}

double Q1(double x, double y)
{
  int q = get_material(x,y);

  double exfl1 = exp(-4*sqr(x))*(y/2.-sqr(y/2.));
  double exfl2 = exfl1 * (1 + sqr(sin(4*M_PI*x)) * sqr(sin(4*M_PI*y))) / 10.0;

  double L = 0.5*exp(-4*sqr(x))*(1+4*(8*sqr(x)-1)*y*(y-2))*D[q][0];
  return L + Sr[q][0]*exfl1 - chi[q][0]*nSf[q][0]*exfl1 - chi[q][1]*nSf[q][1]*exfl2;
}

double Q2(double x, double y)
{
  int q = get_material(x,y);

  double yym2 = (y-2)*y;
  double pi2 = sqr(M_PI), x2 = sqr(x), pix = M_PI*x, piy = M_PI*y;
  double cy2 = sqr(cos(4*piy)),
	 sy2 = sqr(sin(4*piy)),
	 sx2 = sqr(sin(4*pix)),
	 em4x2 = exp(-4*x2);

  double exfl1 = em4x2*(y/2.-sqr(y/2.));
  double exfl2 = exfl1 * (1 + sx2 * sy2) / 10.0;

  double L = 1./20.*em4x2*D[q][1]*(
	     1+4*(8*x2-1)*yym2+16*pi2*yym2*cy2*sx2 + 0.5*sy2*(1-4*(1+4*pi2-8*x2)*yym2 +
             (4*(1+12*pi2-8*x2)*yym2-1)*cos(8*pix) - 64*pix*yym2*sin(8*pix)) + 8*M_PI*(y-1)*sx2*sin(8*piy) );
  return L + Sr[q][1]*exfl2 - Ss[q][1][0]*exfl1;
}

double g1_D(double x, double y) {
  return 0.0;
}

double g2_D(double x, double y) {
  return 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution:

static double exact_flux1(double x, double y, double& dx, double& dy)
{
  double em4x2 = exp(-4*sqr(x));
  dx =  2.0*em4x2*x*y*(y-2);
  dy = -0.5*em4x2*(y-1);
  return em4x2*(y/2.-sqr(y/2.));
}

static double exact_flux2(double x, double y, double& dx, double& dy)
{
  double em4x2 = exp(-4*sqr(x));
  dx = 0.1*em4x2*y*(y-2)*(2*x+(2*x*sqr(sin(4*M_PI*x))-M_PI*sin(8*M_PI*x))*sqr(sin(4*M_PI*y)));
  dy = 0.05*em4x2*(1-y+sqr(sin(4*M_PI*x))*(-(y-1)*sqr(sin(4*M_PI*y))-2*M_PI*y*(y-2)*sin(8*M_PI*y)));
  return em4x2*(y/2.-sqr(y/2.)) * (1 + sqr(sin(4*M_PI*x)) * sqr(sin(4*M_PI*y))) / 10.0;
}

// APPROXIMATE SOLUTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boundary conditions:

const int bc_flux = 1;
const int bc_gamma = 2;
const int bc_symmetry = 3;

BCType bc_types(int marker)
{
	if (marker == bc_flux)
  	return BC_ESSENTIAL;
  else
  	return BC_NATURAL;
}

scalar essential_bc_values_1(int marker, double x, double y)
{
  return g1_D(x, y);
}
scalar essential_bc_values_2(int marker, double x, double y)
{
  return g2_D(x, y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Weak forms:

#include "forms.cpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for calculating errors:

double error_total(double (*efn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*),
       double (*nfn)(MeshFunction*, RefMap*), Tuple<Solution*>& slns1, Tuple<Solution*>& slns2  )
{
  Tuple<Solution*>::iterator it1, it2;
  double error = 0.0, norm = 0.0;

  for (it1=slns1.begin(), it2=slns2.begin(); it1 < slns1.end(); it1++, it2++) {
    assert(it2 < slns2.end());
    error += sqr(calc_abs_error(efn, *it1, *it2));
    if (nfn) norm += sqr(calc_norm(nfn, *it2));
  }

  return (nfn ? sqrt(error/norm) : sqrt(error));
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
	Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();

  int n = 0;

  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);

		for (int i = 0; i < np; i++)
			if (uval[i] < -1e-12)
				n++;
  }

  return n;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh1, mesh2;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh1);

  // Obtain meshes for the 2nd group by cloning the mesh loaded for the 1st group.
  mesh2.copy(&mesh1);

  // Initial uniform refinements.
  for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();

  // Solution variables.
  Solution sln1, sln2;		      // Coarse mesh solution.
  Solution ref_sln1, ref_sln2;	// Reference solution.

  // Create H1 space with default shapesets.
  H1Space space1(&mesh1, bc_types, essential_bc_values_1, P_INIT[0]);
  H1Space space2(&mesh2, bc_types, essential_bc_values_2, P_INIT[0]);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(biform_0_0), H2D_SYM);
  wf.add_matrix_form(0, 1, callback(biform_0_1));
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(1, 1, callback(biform_1_1), H2D_SYM);
  wf.add_vector_form(0, liform_0, liform_0_ord);
  wf.add_vector_form(1, liform_1, liform_1_ord);
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), bc_gamma);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), bc_gamma);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  //selector.set_option(H2D_PREFER_SYMMETRIC_MESH, false);
  //selector.set_error_weights(2.25, 1, sqrt(2.0));

  //////////////////////////////  Adaptivity loop  /////////////////////////////

  // Start time measurement.
  cpu_time.tick();

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(SOLVER_UMFPACK, get_num_dofs(Tuple<Space *>(&space1, &space2)), mat, rhs, solver);
  
  double cta;
  int order_increase = 1;
  int iadapt;
  double error_h1;
  for (iadapt = 0; iadapt < MAX_ADAPT_NUM; iadapt++) {
    
    int ndof = get_num_dofs(Tuple<Space *>(&space1, &space2));
    if (ndof >= NDOF_STOP) break;

    cpu_time.tick();
    info("!---- Adaptivity step %d ---------------------------------------------", iadapt);
    cpu_time.tick(HERMES_SKIP);

    // Construct globally refined reference meshes.
    Mesh ref_mesh1, ref_mesh2;
    ref_mesh1.copy(&mesh1);
    ref_mesh2.copy(&mesh2);
    ref_mesh1.refine_all_elements();
    ref_mesh2.refine_all_elements();

    // Setup spaces for the reference solution.
    Space *ref_space1 = space1.dup(&ref_mesh1);
    Space *ref_space2 = space2.dup(&ref_mesh2);
    int order_increase = 1;
    ref_space1->copy_orders(&space1, order_increase);
    ref_space2->copy_orders(&space2, order_increase);

    cpu_time.tick();
    int ref_ndof = get_num_dofs(Tuple<Space *>(ref_space1, ref_space2));
    info("------------------ Reference solution; NDOF=%d -------------------", ref_ndof);
    cpu_time.tick(HERMES_SKIP);

    // Solve the reference problem.
    solve_linear(Tuple<Space *>(ref_space1, ref_space2), &wf,
                 Tuple<Solution *>(&ref_sln1, &ref_sln2), SOLVER_UMFPACK);
    
    // Project the reference solution on the new coarse mesh.
    cpu_time.tick();
    info("---- Projecting reference solution on new coarse mesh; NDOF=%d ----", ndof);
    cpu_time.tick(HERMES_SKIP);
    project_global(Tuple<Space *>(&space1, &space2), Tuple<int>(H2D_H1_NORM, H2D_H1_NORM), 
                   Tuple<MeshFunction*>(&ref_sln1, &ref_sln2),
                   Tuple<Solution*>(&sln1, &sln2));

    // Calculate element errors and total error estimate.

    Adapt hp(Tuple<Space*>(&space1, &space2), Tuple<int>(H2D_H1_NORM, H2D_H1_NORM));  
    if (ADAPTIVITY_NORM == 2) {
      hp.set_error_form(0, 0, callback(biform_0_0));
      hp.set_error_form(0, 1, callback(biform_0_1));
      hp.set_error_form(1, 0, callback(biform_1_0));
      hp.set_error_form(1, 1, callback(biform_1_1));
    } else {
      if (ADAPTIVITY_NORM == 1) {
        hp.set_error_form(0, 0, callback(biform_0_0));
        hp.set_error_form(1, 1, callback(biform_1_1));
      }
    }

    Tuple<Solution*> slns(&sln1, &sln2);
    Tuple<Solution*> slns_ref(&ref_sln1, &ref_sln2);

    hp.set_solutions(slns, slns_ref);

    double err_est = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
    double err_est_h1 = error_total(error_fn_h1, norm_fn_h1, slns, slns_ref) * 100;

    // Report results.
    cpu_time.tick();
    cta = cpu_time.accumulated();

    // Error w.r.t. the exact solution.
    ExactSolution ex1(&mesh1, exact_flux1), ex2(&mesh2, exact_flux2);
    DiffFilter err_distrib_1(&ex1, &sln1);
    DiffFilter err_distrib_2(&ex2, &sln2);

    double err_exact_h1_1 = calc_rel_error(&ex1, &sln1, H2D_H1_NORM) * 100;
    double err_exact_h1_2 = calc_rel_error(&ex2, &sln2, H2D_H1_NORM) * 100;;

    Tuple<Solution*> exslns(&ex1, &ex2);
    error_h1 = error_total(error_fn_h1, norm_fn_h1, slns, exslns) * 100;

    info("Per-component error wrt. exact solution (H1 norm): %g%%, %g%%", err_exact_h1_1, err_exact_h1_2);
    info("Total error wrt. exact solution (H1 norm): %g%%", error_h1);
    info("Total error wrt. ref. solution  (H1 norm): %g%%", err_est_h1);
    info("Total error wrt. ref. solution  (E norm):  %g%%", err_est);
    
    // Report the number of negative values.
    info("Num. of negative values: %d, %d", get_num_of_neg(&sln1), get_num_of_neg(&sln2));

    cpu_time.tick(HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) break;
    else hp.adapt(Tuple<RefinementSelectors::Selector*>(&selector,&selector), 
                  THRESHOLD, STRATEGY,  MESH_REGULARITY);
  }

  cpu_time.tick();
  cta = cpu_time.accumulated();
  info("Total running time: %g s", cta);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_1 = space1.get_num_dofs(), 
      n_dof_2 = space2.get_num_dofs();
  int n_dof_1_allowed = 400, 
      n_dof_2_allowed = 2900;
      
  int n_neg = get_num_of_neg(&sln1) + get_num_of_neg(&sln2);
  int n_neg_allowed = 0;
  
  int n_iter = iadapt;
  int n_iter_allowed = 20;

  double error = error_h1;
  double error_allowed = 4.5;
  
  printf("n_dof_actual  = %d,%d\n", n_dof_1, n_dof_2);
  printf("n_dof_allowed = %d,%d\n", n_dof_1_allowed, n_dof_2_allowed);
  printf("n_iter_actual = %d\n", n_iter);
  printf("n_iter_allowed= %d\n", n_iter_allowed);
  printf("n_neg_actual  = %d\n", n_neg);
  printf("n_neg_allowed = %d\n", n_neg_allowed);
  printf("error_actual  = %g\n", error);
  printf("error_allowed = %g\n", error_allowed);  
  
  if (   n_dof_1 <= n_dof_1_allowed && n_dof_2 <= n_dof_2_allowed 
      && n_neg <= n_neg_allowed
      && n_iter <= n_iter_allowed
      && error <= error_allowed )   {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
