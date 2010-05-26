#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This example uses automatic adaptivity to solve a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics in a general coordinate system is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. Identifying r = x, z = y, the gradient in the weak form has the same components as in the 
// x-y system, while all integrands are multiplied by 2\pi x (determinant of the transformation matrix).
//
// BC:
//
// homogeneous neumann on symmetry axis
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//

#include "h1_adapt_norm_maxim.h"


// Adaptivity control

const int P_INIT[4] = 
	{1, 1, 1, 1};						 			 // Initial polynomial orders for the individual solution components
const int INIT_REF_NUM[4] = 
	{1, 1, 1, 1};			 						 // Initial uniform mesh refinement for the individual solution components
	
const double THRESHOLD = 0.3;    // This is a quantitative parameter of the adapt(...) function and
                                 // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          // Adaptive strategy:
                                 // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                 //   error is processed. If more elements have similar errors, refine
                                 //   all to keep the mesh symmetric.
                                 // STRATEGY = 1 ... refine all elements whose error is larger
                                 //   than THRESHOLD times maximum element error.
                                 // STRATEGY = 2 ... refine all elements whose error is larger
                                 //   than THRESHOLD.
                                 // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;			   // Maximum allowed element degree
const double ERR_STOP = 0.01;    // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPTATIONS = 100; // Adaptivity process stops when the number of adaptation steps grows over
                                 // this limit.
const NormType NORM_TYPE = NORM_EUCLEDIAN; // When H1 norm is used to measure the error, this parameter specifies the
																 // norm used for normalizing the error (so that solution components with different
																 // scales get same attention). Ignored when energetic norm is used to estimate the error.
																 // NORM_EUCLEDIAN = 0	...	H1 norm
																 // NORM_MAXIMA = 1	...	L_infty norm
const bool PI_ON_COARSE_MESH = false;
																 
const bool jfnk = false;     // true = jacobian-free method,
                            // false = Newton
const bool precond = true;  // preconditioning by jacobian in case of jfnk,
                            // default ML proconditioner in case of Newton

// Macro for simpler definition of bilinear forms in the energy norm
#define callback_egnorm(a)     a<scalar, scalar>, a<Ord, Ord>


// Area markers
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary indices
const int bc_vacuum = 1;
const int bc_sym = 2;

// Boundary condition types
int bc_types(int marker)
{
  return BC_NATURAL;
}

// Reflector properties (0) core properties (1),
const double D[2][4] = {{0.0164, 0.0085, 0.00832, 0.00821},
                        {0.0235, 0.0121, 0.0119, 0.0116}};
const double Sa[2][4] = {{0.00139, 0.000218, 0.00197, 0.0106},
                         {0.00977, 0.162, 0.156, 0.535}};
const double Sr[2][4] = {{1.77139, 0.533218, 3.31197, 0.0106},
                         {1.23977, 0.529, 2.436, 0.535}};
const double Sf[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.00395, 0.0262, 0.0718, 0.346}};
const double nu[2][4] = {{0.0, 0.0, 0.0, 0.0}, {2.49, 2.43, 2.42, 2.42}};
const double chi[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.9675, 0.03250, 0.0, 0.0}};
const double Ss[2][4][4] = {{{ 0.0,   0.0,  0.0, 0.0},
                             {1.77,   0.0,  0.0, 0.0},
                             { 0.0, 0.533,  0.0, 0.0},
                             { 0.0,   0.0, 3.31, 0.0}},
                            {{ 0.0,   0.0,  0.0, 0.0},
                             {1.23,   0.0,  0.0, 0.0},
                             { 0.0, 0.367,  0.0, 0.0},
                             { 0.0,   0.0, 2.28, 0.0}}};

// Power iteration control
																
double k_eff = 1.0;					// initial eigenvalue approximation
double TOL_PIT_CM = 1e-5;		// tolerance for eigenvalue convergence when solving on coarse mesh
double TOL_PIT_RM = 1e-6;		// tolerance for eigenvalue convergence when solving on reference mesh

TimePeriod cpu_time;

// Weak forms
#include "forms.cpp"

#include "CustomProjector.cpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Source function
void source_fn(int n, scalar* a, scalar* b, scalar* c, scalar* d, scalar* out)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = (nu[1][0] * Sf[1][0] * a[i] +
        nu[1][1] * Sf[1][1] * b[i] +
        nu[1][2] * Sf[1][2] * c[i] +
        nu[1][3] * Sf[1][3] * d[i]);
  }
}

// Integral over the active core
double integrate(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }

  return 2.0 * M_PI * integral;
}

// Calculate number of negative solution values
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
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
			
		for (int i = 0; i < np; i++)
			if (uval[i] < -1e-12)
				n++;
  }
  
  return n;
}

// Power iteration using "sys" as the linear system associated with the generalized 
// eigenvalue problem and "iter_X" as the initial guess for eigenvector; "sys" 
// is assumed to be already assembled and only rhs updates are performed; converged
// eigenvectors are stored in "sln_X" and the eigenvalue in global variable "k_eff"
void power_iteration( Solution *sln1, Solution *sln2, Solution *sln3, Solution *sln4,
											Solution *iter1, Solution *iter2, Solution *iter3, Solution *iter4,
											FeProblem *fep, NoxSolver *solver, double *vec, double tol)
{
	bool eigen_done = false; int it = 0;
	
	Space *space1 = fep->get_space(0);
	Space *space2 = fep->get_space(1);
	Space *space3 = fep->get_space(2);
	Space *space4 = fep->get_space(3);
	
	PrecalcShapeset *pss1 = fep->get_pss(0);
	PrecalcShapeset *pss2 = fep->get_pss(1);
	PrecalcShapeset *pss3 = fep->get_pss(2);
	PrecalcShapeset *pss4 = fep->get_pss(3);
	
	iter1->set_fe_solution(space1, pss1, vec);
  iter2->set_fe_solution(space2, pss2, vec);
  iter3->set_fe_solution(space3, pss3, vec);
  iter4->set_fe_solution(space4, pss4, vec);
    
	do
	{
	  info("------------------------ Power iteration %d ------------------------", it++);
	  cpu_time.tick(H2D_SKIP);
	  
		// solve for new eigenvectors
		solver->set_init_sln(vec);
info("initialized");
	  bool solved = solver->solve();
info("solved");
		if (solved) {
			vec = solver->get_solution();
      sln1->set_fe_solution(space1, pss1, vec);
  		sln2->set_fe_solution(space2, pss2, vec);
		  sln3->set_fe_solution(space3, pss3, vec);
		  sln4->set_fe_solution(space4, pss4, vec);

			cpu_time.tick();
      info("Number of nonlin iters: %d (norm of residual: %g)",
          solver->get_num_iters(), solver->get_residual());
      info("Total number of iters in linsolver: %d (achieved tolerance in the last step: %g)",
          solver->get_num_lin_iters(), solver->get_achieved_tol());
      cpu_time.tick(H2D_SKIP);
      
			// update fission sources
			SimpleFilter source(source_fn, sln1, sln2, sln3, sln4);
			SimpleFilter source_prev(source_fn, iter1, iter2, iter3, iter4);

		  // compute eigenvalue
	  	double k_new = k_eff * (integrate(&source, marker_core) / integrate(&source_prev, marker_core));
	  
			cpu_time.tick();
			info("Dominant eigenvalue (est): %g, rel error: %g", k_new, fabs((k_eff - k_new) / k_new));
			cpu_time.tick(H2D_SKIP);

			// stopping criterion
			if (fabs((k_eff - k_new) / k_new) < tol) eigen_done = true;

			// store eigenvectors for next iteration
			iter1 = sln1; iter2 = sln2; iter3 = sln3; iter4 = sln4;
			// iter1->copy(sln1);    iter2->copy(sln2);
			// iter3->copy(sln3);    iter4->copy(sln4);
	  
			// update eigenvalue
			k_eff = k_new;
		}
		else
      error("Power iteration failed.");
	 
	  cpu_time.tick();
	}
	while (!eigen_done);
}
										

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh1, mesh2, mesh3, mesh4;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh1);

	// obtain meshes for the 2nd to 4th group by cloning the mesh loaded for the 1st group
	mesh2.copy(&mesh1);
	mesh3.copy(&mesh1);
	mesh4.copy(&mesh1);

  // initial uniform refinements
  for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[2]; i++) mesh3.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[3]; i++) mesh4.refine_all_elements();
  	
  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);
  PrecalcShapeset pss3(&shapeset);
  PrecalcShapeset pss4(&shapeset);

  // solution variables
  Solution iter1, iter2, iter3, iter4;													// previous iterates
  Solution sln1, sln2, sln3, sln4;	// current iterate on coarse mesh
	Solution sln1_ref, sln2_ref, sln3_ref, sln4_ref;					// current iterate on fine mesh
	// set initial iteration guess
  iter1.set_const(&mesh1, 1.00);
  iter2.set_const(&mesh2, 1.00);
  iter3.set_const(&mesh3, 1.00);
  iter4.set_const(&mesh4, 1.00);

  // create finite element spaces
  H1Space space1(&mesh1, &shapeset);
  H1Space space2(&mesh2, &shapeset);
  H1Space space3(&mesh3, &shapeset);
  H1Space space4(&mesh4, &shapeset);
  space1.set_bc_types(bc_types);
  space2.set_bc_types(bc_types);
  space3.set_bc_types(bc_types);
  space4.set_bc_types(bc_types);
  space1.set_uniform_order(P_INIT[0]);
  space2.set_uniform_order(P_INIT[1]);
  space3.set_uniform_order(P_INIT[2]);
  space4.set_uniform_order(P_INIT[3]);

  // initialize the weak formulation
  WeakForm wf(4, jfnk);
  wf.add_jacform(0, 0, callback(jacobian_0_0));
  wf.add_jacform(1, 1, callback(jacobian_1_1));
  wf.add_jacform(1, 0, callback(jacobian_1_0));
  wf.add_jacform(2, 2, callback(jacobian_2_2));
  wf.add_jacform(2, 1, callback(jacobian_2_1));
  wf.add_jacform(3, 3, callback(jacobian_3_3));
  wf.add_jacform(3, 2, callback(jacobian_3_2));
  wf.add_resform(0, callback(residual_0), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_resform(1, callback(residual_1), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_resform(2, callback(residual_2), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_resform(3, callback(residual_3), marker_core, 4, &iter1, &iter2, &iter3, &iter4);
  wf.add_jacform_surf(0, 0, callback(jacobian_surf_0_0), bc_vacuum);
  wf.add_jacform_surf(1, 1, callback(jacobian_surf_1_1), bc_vacuum);
  wf.add_jacform_surf(2, 2, callback(jacobian_surf_2_2), bc_vacuum);
  wf.add_jacform_surf(3, 3, callback(jacobian_surf_3_3), bc_vacuum);

  // visualization
  ScalarView view1("Neutron flux 1", 0, 0, 320, 400);
  ScalarView view2("Neutron flux 2", 330, 0, 320, 400);
  ScalarView view3("Neutron flux 3", 660, 0, 320, 400);
  ScalarView view4("Neutron flux 4", 990, 0, 320, 400);
  OrderView oview1("Mesh for group 1", 0, 450, 320, 500);
  OrderView oview2("Mesh for group 2", 330, 450, 320, 500);
  OrderView oview3("Mesh for group 3", 660, 450, 320, 500);
  OrderView oview4("Mesh for group 4", 990, 450, 320, 500);
/* 
  ScalarView view1("Neutron flux 1", 0, 0, 320, 400);
//  ScalarView view2("Neutron flux 2", 330, 0, 320, 400);
  ScalarView view3("Neutron flux 3", 330, 0, 320, 400);
  ScalarView view4("Neutron flux 4", 660, 0, 320, 400);
  OrderView oview1("Mesh for group 1", 0, 450, 320, 500);
  OrderView oview2("Mesh for group 2", 330, 450, 320, 500);
  OrderView oview3("Mesh for group 3", 0, 1000, 320, 500);
  OrderView oview4("Mesh for group 4", 330, 1000, 320, 500);
  /*view1.show_mesh(false);
  view2.show_mesh(false);
  view3.show_mesh(false);
  view4.show_mesh(false);
  */																 

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu, graph_dof_keff, graph_cpu_keff;

  // prepare selector
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, MAX_ORDER, &shapeset);

  // adaptivity loop
	// start time measurement
	
	cpu_time.tick();
	
	// enumerate DoF
	int ndof = assign_dofs(4, &space1, &space2, &space3, &space4);
		
	// initial coarse mesh solution
	FeProblem fep(&wf);
	fep.set_spaces(4, &space1, &space2, &space3, &space4);
	fep.set_pss(4, &pss1, &pss2, &pss3, &pss4);

	// Obtain solution vector for initial guess
	cpu_time.tick();
	info("Projecting initial solution");
	cpu_time.tick(H2D_SKIP);
  
  CustomProjector proj(&iter1, &iter2, &iter3, &iter4,
  									 	 &space1, &space2, &space3, &space4,
  									 	 &pss1, &pss2, &pss3, &pss4	);

  UmfpackSolver umfpack;
  proj.set_solver(&umfpack);
  double* vec = proj.project();


  // Solver + preconditioner
  NoxSolver solver(&fep);
  MlPrecond pc("sa");
  if (precond)
  {
    if (jfnk) solver.set_precond(&pc);
    else solver.set_precond("ML");
  }
	solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest | 
                            NOX::Utils::LinearSolverDetails);
	power_iteration(&sln1, &sln2, &sln3, &sln4, 
									&iter1, &iter2, &iter3, &iter4,
									&fep, &solver, vec, TOL_PIT_CM);
									
	for (int iadapt = 0; iadapt < MAX_ADAPTATIONS; iadapt++) {
		cpu_time.tick();
		
		info("!---- Adaptivity step %d ---------------------------------------------", iadapt);
			
		// view the solution and meshes

    info("flux1_dof=%d, flux2_dof=%d, flux3_dof=%d, flux4_dof=%d", 
    			space1.get_num_dofs(), space2.get_num_dofs(), space3.get_num_dofs(), space4.get_num_dofs());
    view1.show(&sln1);
    view2.show(&sln2);
    view3.show(&sln3);
    view4.show(&sln4);
    oview1.show(&space1);
    oview2.show(&space2);
    oview3.show(&space3);
    oview4.show(&space4);

		info("Num. of negative values: %d, %d, %d, %d", 
					get_num_of_neg(&sln1), get_num_of_neg(&sln2), get_num_of_neg(&sln3), get_num_of_neg(&sln4));		
		
    cpu_time.tick(H2D_SKIP);
    
    // solve the fine mesh problem

    Mesh refmesh1, refmesh2, refmesh3, refmesh4;
    
    refmesh1.copy(&mesh1);
    refmesh2.copy(&mesh2);
    refmesh3.copy(&mesh3);
    refmesh4.copy(&mesh4);
    refmesh1.refine_all_elements();
    refmesh2.refine_all_elements();
    refmesh3.refine_all_elements();
    refmesh4.refine_all_elements();
    
    H1Space refspace1(&refmesh1, &shapeset);
    H1Space refspace2(&refmesh2, &shapeset);
    H1Space refspace3(&refmesh3, &shapeset);
    H1Space refspace4(&refmesh4, &shapeset);
    
    refspace1.set_bc_types(bc_types);
    refspace1.copy_orders(&space1, 1);
    refspace2.set_bc_types(bc_types);
    refspace2.copy_orders(&space2, 1);
    refspace3.set_bc_types(bc_types);
    refspace3.copy_orders(&space3, 1);
    refspace4.set_bc_types(bc_types);
    refspace4.copy_orders(&space4, 1);
        
    // enumerate degrees of freedom
    int ndof_ref = assign_dofs(4, &refspace1, &refspace2, &refspace3, &refspace4);

    FeProblem reffep(&wf);
    reffep.set_spaces(4, &refspace1, &refspace2, &refspace3, &refspace4);
    reffep.set_pss(4, &pss1, &pss2, &pss3, &pss4);
    
    NoxSolver refsolver(&reffep);
    if (precond)
    {
      if (jfnk) refsolver.set_precond(&pc);
      else refsolver.set_precond("ML");
    }
    
    cpu_time.tick();    
		info("---------- Reference mesh power iteration; NDOF=%d ----------------", ndof_ref);	
    cpu_time.tick(H2D_SKIP);
    
    
   	if (iadapt == 0) {
    	CustomProjector proj(	&sln1, &sln2, &sln3, &sln4,
      											&refspace1, &refspace2, &refspace3, &refspace4,
      											&pss1, &pss2, &pss3, &pss4 );
		  
		  UmfpackSolver umfpack;
		  proj.set_solver(&umfpack);
		  
			vec = proj.project();
    }
    else {
    	CustomProjector proj(	&sln1_ref, &sln2_ref, &sln3_ref, &sln4_ref,
      											&refspace1, &refspace2, &refspace3, &refspace4,
      											&pss1, &pss2, &pss3, &pss4 );
		  
		  UmfpackSolver umfpack;
		  proj.set_solver(&umfpack);
		  
			vec = proj.project();
    }

		power_iteration(&sln1_ref, &sln2_ref, &sln3_ref, &sln4_ref,
										&iter1, &iter2, &iter3, &iter4,
										&reffep, &refsolver, vec, TOL_PIT_RM);
		
		// calculate element errors and total error estimate

		// H1 normalized by H1 or L_inf
    H1AdaptNormMaxim hp(NORM_TYPE, Tuple<Space*>(&space1, &space2, &space3, &space4));
    /*hp.set_biform(0, 0, callback_egnorm(biform_0_0));
		hp.set_biform(1, 1, callback_egnorm(biform_1_1));
		hp.set_biform(1, 0, callback_egnorm(biform_1_0));
		hp.set_biform(2, 2, callback_egnorm(biform_2_2));
		hp.set_biform(2, 1, callback_egnorm(biform_2_1));
		hp.set_biform(3, 3, callback_egnorm(biform_3_3));
		hp.set_biform(3, 2, callback_egnorm(biform_3_2));
		*/
		/*for (int i = 0; i < 4; i++)
			hp.set_biform(i,i,callback_egnorm(projection_biform));
		*/	
		// how to include surface elements in the energy norm ?
		
    hp.set_solutions(Tuple<Solution*>(&sln1, &sln2, &sln3, &sln4),
      Tuple<Solution*>(&sln1_ref, &sln2_ref, &sln3_ref, &sln4_ref));
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  
	  // eigenvalue error w.r.t. solution obtained on a 3x uniformly refined mesh
  	// with uniform distribution of polynomial degrees (=3), converged to within
  	// tolerance of 1e-7; in units of percent-milli (pcm)
  	double keff_err = 1e5*fabs(k_eff - 1.140911)/1.140911;
  	
    // report results
    
    cpu_time.tick();        
    
    double cta = cpu_time.accumulated();

    info("Estimate of error: %g%%; exact eigenvalue error: %g [pcm]", err_est, keff_err);
  					
    // add entry to DOF convergence graph
    graph_dof.add_values(ndof, err_est);
    graph_dof.save("conv_dof.dat");
  
    // add entry to CPU convergence graph
    graph_cpu.add_values(cta, err_est);
    graph_cpu.save("conv_cpu.dat");
    
    // add entry to DOF convergence graph w.r.t. dominant eigenvalue
    graph_dof_keff.add_values(ndof, keff_err);
    graph_dof_keff.save("conv_dof_keff.dat");
    
    // add entry to CPU convergence graph w.r.t. dominant eigenvalue
    graph_cpu_keff.add_values(cta, keff_err);
    graph_cpu_keff.save("conv_cpu_keff.dat");
    
    cpu_time.tick(H2D_SKIP);
    
    // alternative error measure 
    //if (keff_err < 0.1)
    //	break;
    	
    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) {
    	break;
    }
    else {
    	hp.adapt(&selector, THRESHOLD, STRATEGY,  MESH_REGULARITY);
    	
    	ndof = assign_dofs(4, &space1, &space2, &space3, &space4);
    	
    	if (ndof >= NDOF_STOP)
    		break;
    		
  		// project the fine mesh solution on the new coarse mesh
			cpu_time.tick();
      info("---- Projecting fine mesh solution on new coarse mesh -----------------");
      cpu_time.tick(H2D_SKIP);
      
      CustomProjector proj(	&sln1_ref, &sln2_ref, &sln3_ref, &sln4_ref,
		    										&space1, &space2, &space3, &space4, 
		    										&pss1, &pss2, &pss3, &pss4 );
		  
		  UmfpackSolver umfpack;
		  proj.set_solver(&umfpack);
		  
			vec = proj.project();
    	
    	if (PI_ON_COARSE_MESH) {
				// run the power iteration on the new coarse mesh
				cpu_time.tick();	
				info("----------- Coarse mesh power iteration; NDOF=%d -----------------", ndof);	  
				cpu_time.tick(H2D_SKIP);
				
				fep.set_spaces(4, &space1, &space2, &space3, &space4);
				fep.set_pss(4, &pss1, &pss2, &pss3, &pss4);
				
				NoxSolver solver(&fep);
				if (precond)
				{
				  if (jfnk) solver.set_precond(&pc);
				  else solver.set_precond("ML");
				}
   				
				power_iteration(&sln1, &sln2, &sln3, &sln4, 
												&iter1, &iter2, &iter3, &iter4,
												&fep, &solver, vec, TOL_PIT_CM);
			}
			else {
				sln1.set_fe_solution(&space1, &pss1, vec);
  			sln2.set_fe_solution(&space2, &pss2, vec);
		  	sln3.set_fe_solution(&space3, &pss3, vec);
		  	sln4.set_fe_solution(&space4, &pss4, vec);
			}		
		}
	}
	
	cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // wait for all views to be closed
  View::wait();
  return 0;
};

