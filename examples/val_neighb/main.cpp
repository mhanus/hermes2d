#include "hermes2d.h"
#include "solver_umfpack.h"
#define H2D_REPORT_FILE "test.log"

/*
 * This example is only for showing how to use the values from neighbor in linear surface forms.
 * Solution of the original problem is function x^3 + y^3. As simulated input from previous step (the need to have values from neighbor is originally from time dependent
 * problem) is taken exact solution.
 * The example doesn't have any real meaning.
 */


const int P_INIT = 1;          // Polynomial degree of mesh elements
const int INIT_REF_NUM = 3;    // Number of initial uniform mesh refinements

scalar proj_func(double x, double y, double &dx, double &dy)
{
  return x*x*x + y*y*y;
}

// bilinear and linear form defining the projection
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += wt[i] * ((pow(e->x[i], 1) + pow(e->y[i], 1)) * v->val[i]);
  }
  return result;
}


template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  // here to get values from neighbor you have to use method get_fn_neighbor(). This is for safety.
  for (int i = 0; i < n; i++){
     result += 0.5*wt[i] * (ext->fn[0]->val[i] + ext->get_fn_neighbor(0)->val[i]) * v->val[i];
  }
  return result;
}

// boundary conditions
BCType bc_types(int marker)
{
   return BC_NONE;
}
// function values for Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

int main(int argc, char* argv[])
{
  char *f;
  char *def_f = "square.mesh";

  if (argc < 2) {
    info("Using square.mesh for the meshfile.");
    f = def_f;
  }  else {
    f = argv[1];
  }

  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load(f, &mesh);

  // uniform mesh refinements
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  mesh.refine_element(29);
  mesh.refine_element(87);

  // display the mesh
/*	 MeshView mview("info_neighbor", 100, 100, 500, 500);
	 mview.show(&mesh);
  // wait for keyboard or mouse input
   View::wait("Waiting for keyboard or mouse input.");
*/
  // create the L2 space
  L2Space space(&mesh,P_INIT);
  space.set_bc_types(bc_types);

  BaseView bview;
  bview.show(&space);
  bview.wait_for_close();

  Solution sln;
  Solution xprev;
  xprev.set_exact(&mesh, &proj_func);
  // matrix solver
  UmfpackSolver umfpack;

  // initialize the weak formulation
    WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  // if you want to work in linear surface form with values from neighbors use flag H2D_ANY_INNER_EDGE.
  wf.add_vector_form_surf(callback(linear_form_surf), H2D_ANY_INNER_EDGE, Tuple<MeshFunction* >(&xprev));

  // assemble and solve the finite element problem
  LinSystem sys(&wf, &umfpack, &space);
  sys.assemble();
  sys.solve(&sln);

  // visualize the solution
  ScalarView view1("Solution");
  view1.show(&sln);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

