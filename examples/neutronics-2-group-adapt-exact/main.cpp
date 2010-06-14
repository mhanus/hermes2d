#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"
#include "solver_umfpack.h"
#include "h1_adapt_custom_norm.h"
#include "TikZGraph.h"
#include "ScaledScalarView.h"

using namespace RefinementSelectors;

// This example uses automatic adaptivity to solve a 2-group neutron diffusion equation with a known exact solution.
// The solution reflects the typical behavior observed in real cases, where one component is very smooth and the 
// other more oscillating. Typical boundary conditions prescribed in real models have also been chosen.
//
// EQUATION:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g 
//		- \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'}	- \sum_{g'} \nu\Sigma_f^{g'} \phi_{g'} = Q_g
//
// BC:
//
// homogeneous neumann on symmetry axes
// homogeneous dirichlet on zero flux boundary
// -d D_g\phi_g / d n = 8 \phi_g   on albedo boundary (homogeneous Robin)
//

// INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptivity control

const bool SOLVE_ON_COARSE_MESH = false; 					// If true, coarse mesh FE problem is solved in every adaptivity step.
                                         					// If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT[2] = 
	{1, 1};						 			 				// Initial polynomial orders for the individual solution components
const int INIT_REF_NUM[2] = 
	{1, 1};			 						 				// Initial uniform mesh refinement for the individual solution components
	
const double THRESHOLD = 0.3;    	// This is a quantitative parameter of the adapt(...) function and
                                 	// it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          	// Adaptive strategy:
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
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;	         	 // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
const NormType NORM_TYPE = NORM_DEFAULT; // Specifies norm used by the class H1AdaptCustomNorm to calculate relative error 
																				 // NORM_DEFAULT = 0 ... the same norm as used by H1Adapt 
																				 //		(defined either by the supplied bilinear forms or by the standard H1 norm)
																				 // NORM_MAX = 1 ... maximum (L_infty) norm
const int ENERGY_NORM = 2;							 // Specifies the norm used by H1AdaptCustomNorm to calculate the error 
																				 // (and its normalization if NORM_TYPE == NORM_DEFAULT)
																				 // ENERGY_NORM = 0 ... H1 norm
																				 // ENERGY_NORM = 1 ... norm defined by the diagonal parts of the bilinear form																				 
																				 // ENERGY_NORM = 2 ... energy norm defined by the full non-symmetric bilinear form

// Macro for simpler definition of bilinear forms in the energy norm
#define callback_egnorm(a)     a<scalar, scalar>, a<Ord, Ord>
		
// Variables used for reporting of results														 
TimePeriod cpu_time;							// Time measurements
const int ERR_PLOT = 0;						// Row in the convergence graphs for exact errors 
const int ERR_EST_PLOT = 1;				// Row in the convergence graphs for error estimates 
const int GROUP_1 = 0;						// Row in the NDOFs evolution graph for group 1
const int GROUP_2 = 1;						// Row in the NDOFs evolution graph for group 2


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem data

// Two-group material properties for the 4 macroregions
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
const double Ss[4][2][2] = {	{	{	0.0, 0.0 },
                             		{ 0.05, 0.0 }	},
                             	{	{	0.0, 0.0 },
                             		{ 0.08, 0.0 }	},
                             	{	{	0.0, 0.0 },
                             		{ 0.025, 0.0 }	},
                             	{	{	0.0, 0.0 },
                             		{ 0.014, 0.0 }	}	};

double a = 0., b = 1., c = (a+b)/2.;

inline int get_material(double x, double y) {
	if (x >= a && x <= c && y >= a && y <= c) return 0;
	if (x >= c && x <= b && y >= a && y <= c) return 1;
	if (x >= c && x <= b && y >= c && y <= b) return 2;
	if (x >= a && x <= c && y >= c && y <= b) return 3;
}

double Q1(double x, double y) {
	int q = get_material(x,y);
	
	double exfl1 = exp(-4*sqr(x))*(y/2.-sqr(y/2.));
	double exfl2 = exfl1 * (1 + sqr(sin(4*M_PI*x)) * sqr(sin(4*M_PI*y))) / 10.0;
	
	double L = 0.5*exp(-4*sqr(x))*(1+4*(8*sqr(x)-1)*y*(y-2))*D[q][0];
	return L + Sr[q][0]*exfl1 - chi[q][0]*nSf[q][0]*exfl1 - chi[q][1]*nSf[q][1]*exfl2;
}
double Q2(double x, double y) {
	int q = get_material(x,y);

	double yym2 = (y-2)*y;
	double pi2 = sqr(M_PI), x2 = sqr(x), pix = M_PI*x, piy = M_PI*y;
	double 	cy2 = sqr(cos(4*piy)),
					sy2 = sqr(sin(4*piy)),
					sx2 = sqr(sin(4*pix)),
					em4x2 = exp(-4*x2);
					
	double exfl1 = em4x2*(y/2.-sqr(y/2.));
	double exfl2 = exfl1 * (1 + sx2 * sy2) / 10.0;
	
	double L = 1./20.*em4x2*D[q][1]*( 
		1+4*(8*x2-1)*yym2+16*pi2*yym2*cy2*sx2 + 0.5*sy2*(1-4*(1+4*pi2-8*x2)*yym2 + (4*(1+12*pi2-8*x2)*yym2-1)*cos(8*pix) - 64*pix*yym2*sin(8*pix)) + 8*M_PI*(y-1)*sx2*sin(8*piy) );
	return L + Sr[q][1]*exfl2 - Ss[q][1][0]*exfl1;
}

double g1_D(double x, double y) {
  return 0.0;
}
double g2_D(double x, double y) {
  return 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution

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
// Boundary conditions

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

scalar bc_values1(int marker, double x, double y)
{
  return g1_D(x, y);
}
scalar bc_values2(int marker, double x, double y)
{
  return g2_D(x, y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Weak forms

#include "forms.cpp"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Additional functions

inline void err_fn(int n, scalar* uh, scalar* u, scalar* out)
{
  for (int i = 0; i < n; i++)  out[i] = u[i] - uh[i];
}

inline double error_fn_h1_semi(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
	return int_h1_semi_error<scalar>(sln1, sln2, ru, rv);
}

inline double norm_fn_h1_semi(MeshFunction* sln, RefMap* ru)
{
	return int_h1_seminorm<scalar>(sln, ru);
}

inline double h1_semi_error(MeshFunction* sln1, MeshFunction* sln2)
{
  return calc_error(error_fn_h1_semi, sln1, sln2) / calc_norm(norm_fn_h1_semi, sln2);
}

inline double h1_semi_norm(MeshFunction* sln)
{
  return calc_norm(norm_fn_h1_semi, sln);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	///////////////////////////////  initialization  /////////////////////////////

  // load the mesh
  Mesh mesh1, mesh2;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh1);

	// obtain meshes for the 2nd group by cloning the mesh loaded for the 1st group
	mesh2.copy(&mesh1);

  // initial uniform refinements
  for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();
  
  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss1(&shapeset);
  PrecalcShapeset pss2(&shapeset);

  // solution variables
  Solution sln1, sln2;					// coarse solution
	Solution sln1_ref, sln2_ref;	// refined solution

  // matrix solver
  UmfpackSolver umfpack;

  // create finite element spaces
  H1Space space1(&mesh1, &shapeset);
  H1Space space2(&mesh2, &shapeset);
  space1.set_bc_types(bc_types);
  space2.set_bc_types(bc_types);
  space1.set_essential_bc_values(bc_values1);
  space2.set_essential_bc_values(bc_values2);
  space1.set_uniform_order(P_INIT[0]);
  space2.set_uniform_order(P_INIT[1]);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(biform_0_0), H2D_SYM);
  wf.add_biform(0, 1, callback(biform_0_1));
  wf.add_biform(1, 0, callback(biform_1_0));
  wf.add_biform(1, 1, callback(biform_1_1), H2D_SYM);
  wf.add_liform(0, liform_0, liform_0_ord);
  wf.add_liform(1, liform_1, liform_1_ord);
  wf.add_biform_surf(0, 0, callback(biform_surf_0_0), bc_gamma);
  wf.add_biform_surf(1, 1, callback(biform_surf_1_1), bc_gamma);

  // visualization
  ScalarView view1("Neutron flux 1", 0, 0, 500, 460);
  ScaledScalarView view2("Neutron flux 2", 500, 0, 500, 460);
  ScalarView view3("Error in neutron flux 1", 1000, 0, 500, 460);
  ScalarView view4("Error in neutron flux 2", 1500, 0, 500, 460);
  //ScalarView view3("Error in neutron flux 1", 0, 0, 500, 460);
  //ScalarView view4("Error in neutron flux 2", 500, 0, 500, 460);
  OrderView oview1("Mesh for group 1", 0, 500, 360, 300);
  OrderView oview2("Mesh for group 2", 360, 500, 360, 300);
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);
    
  // DOF and CPU convergence graphs
  TikZGraph graph_dof("Error convergence", "Degrees of freedom", "Error [\\%]");
  graph_dof.add_row("exact error", "b", "-", "o");
  graph_dof.add_row("est. error ", "r", "-", "s");
  graph_dof.set_log_y();
  graph_dof.show_legend(); 
  graph_dof.show_grid();
  
  TikZGraph graph_dof_evol("Evolution of NDOFs", "Adaptation step", "Num. DOFs");
  graph_dof_evol.add_row("group 1", "b", "-", "o");
  graph_dof_evol.add_row("group 2 ", "r", "-", "s"); 
  graph_dof_evol.set_log_y();
  graph_dof_evol.show_legend(); 
  graph_dof_evol.show_grid();
  
  TikZGraph graph_cpu("Error convergence", "CPU time [s]", "Error [\\%]");
  graph_cpu.add_row("exact error", "b", "-", "o");
  graph_cpu.add_row("est. error ", "r", "-", "s");
  graph_cpu.set_log_y();
  graph_cpu.show_legend(); 
  graph_cpu.show_grid();

  // prepare selector
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);
	//selector.set_option(H2D_PREFER_SYMMETRIC_MESH, false);
	//selector.set_error_weights(2.25, 1, sqrt(2.0));


  //////////////////////////////  adaptivity loop  /////////////////////////////
  
	// start time measurement
	cpu_time.tick();
	
	// enumerate DoF
	int ndof = assign_dofs(2, &space1, &space2);
		
	// initial coarse mesh solution
	LinSystem sys(&wf, &umfpack, 2, &space1, &space2);
	sys.assemble();
	
	double error1, error2, err_est;
	double cta;

	for (int iadapt = 0; iadapt < MAX_ADAPT_NUM; iadapt++) {
		cpu_time.tick();	
		info("!---- Adaptivity step %d ---------------------------------------------", iadapt);
    cpu_time.tick(H2D_SKIP);
    
    // solve the fine mesh problem
 
    RefSystem refsys(&sys);
    refsys.assemble();	
      
    cpu_time.tick();
    
    int ndof_ref =	refsys.get_num_dofs(0) + refsys.get_num_dofs(1);  
		info("---------- Reference mesh solution; NDOF=%d ----------------", ndof_ref);	
		
    cpu_time.tick(H2D_SKIP);
    
		refsys.solve(2, &sln1_ref, &sln2_ref);
		
		if (SOLVE_ON_COARSE_MESH) {
			cpu_time.tick();	
			info("----------- Coarse mesh solution; NDOF=%d -----------------", ndof);	  
			cpu_time.tick(H2D_SKIP);
			 
			sys.assemble();	
			sys.solve(2, &sln1, &sln2);
		}	else {
			// project the fine mesh solution on the new coarse mesh
			cpu_time.tick();
	    info("---- Projecting fine mesh solution on new coarse mesh -----------------");
	    cpu_time.tick(H2D_SKIP);
	    sys.project_global(&sln1_ref, &sln2_ref, &sln1, &sln2);
		}
		
		// calculate element errors and total error estimate
    
    H1AdaptCustomNorm hp(NORM_TYPE, Tuple<Space*>(&space1, &space2));
    if (ENERGY_NORM == 2) {
			hp.set_biform(0, 0, callback_egnorm(biform_0_0));
			hp.set_biform(0, 1, callback_egnorm(biform_0_1));
			hp.set_biform(1, 0, callback_egnorm(biform_1_0));
			hp.set_biform(1, 1, callback_egnorm(biform_1_1));
		} else if (ENERGY_NORM == 1) {
			hp.set_biform(0, 0, callback_egnorm(biform_0_0));
			hp.set_biform(1, 1, callback_egnorm(biform_1_1));
		}
			
    hp.set_solutions(Tuple<Solution*>(&sln1, &sln2),
      Tuple<Solution*>(&sln1_ref, &sln2_ref));
    err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
        	
    // report results
    
    cpu_time.tick();            
   	cta = cpu_time.accumulated();
   	
   	info("flux1_dof=%d, flux2_dof=%d", space1.get_num_dofs(), space2.get_num_dofs());

    view1.show(&sln1);
    view2.show(&sln2);
    
    // error w.r.t. the exact solution
    ExactSolution exact1(&mesh1, exact_flux1);
    ExactSolution exact2(&mesh2, exact_flux2);  
	  error1 = h1_error(&sln1, &exact1) * 100;
    error2 = h1_error(&sln2, &exact2) * 100;
		SimpleFilter err1(err_fn, &sln1, &exact1);
		SimpleFilter err2(err_fn, &sln2, &exact2);
		view3.show(&err1);
    view4.show(&err2);
        
    oview1.show(&space1);
    oview2.show(&space2);
    
    info("Exact solution error for phi_1 (H1 norm): %g%%", error1);
    info("Exact solution error for phi_2 (H1 norm): %g%%", error2);
    info("Estimate of error wrt. ref. solution (energy norm): %g%%", err_est);
    
    double error = std::max(error1,error2);
    
    if (ndof > 100) {				
		  // add entry to DOF convergence graphs
		  graph_dof.add_values(ERR_PLOT, ndof, error);  
		  graph_dof.add_values(ERR_EST_PLOT, ndof, err_est);
		  // add entry to DOF evolution graphs
		  graph_dof_evol.add_values(GROUP_1, iadapt, sys.get_num_dofs(0));
 		  graph_dof_evol.add_values(GROUP_2, iadapt, sys.get_num_dofs(1));
		  // add entry to CPU convergence graphs
		  graph_cpu.add_values(ERR_PLOT, cta, error);
		  graph_cpu.add_values(ERR_EST_PLOT, cta, err_est);
		}
           
    cpu_time.tick(H2D_SKIP);   
    	
    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) {
    	break;
    }
    else {
    	hp.adapt(&selector, THRESHOLD, STRATEGY,  MESH_REGULARITY);	
    	ndof = assign_dofs(2, &space1, &space2);
    	
    	if (ndof >= NDOF_STOP)
    		break;
		}
	}
	
	cpu_time.tick();
	cta = cpu_time.accumulated();
  verbose("Total running time: %g s", cta);
  
  FILE* f = fopen ("results.txt","w");

  fprintf(f, "Total running time: %g s\n", cta);
 	fprintf(f, "Exact solution error for phi_1 (H1 norm): %g%%\n", error1);
  fprintf(f, "Exact solution error for phi_2 (H1 norm): %g%%\n", error2);
  fprintf(f, "Estimate of error wrt. ref. solution (energy norm): %g%%\n", err_est);

	fclose(f);
	
  view2.scale(25);
  
  // save convergence graphs
  graph_dof.save("conv_dof.gp");
  graph_cpu.save("conv_cpu.gp");
  graph_dof_evol.save("dof_evol.gp");
  
  // wait for all views to be closed
  View::wait();
  return 0;
};
