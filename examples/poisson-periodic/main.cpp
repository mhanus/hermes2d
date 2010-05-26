#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE

#include "hermes2d.h"
#include "solver_umfpack.h"  // defines the class UmfpackSolver
#include "H1PeriodicSpace.h"

// This example shows how to solve a first simple PDE:
//   - load the mesh,
//   - perform initial refinements
//   - create a H1 space over the mesh
//   - define weak formulation
//   - initialize matrix solver
//   - assemble and solve the matrix system
//   - visualize the solution
//
// PDE: Poisson equation -Laplace u = CONST_F with homogeneous (zero)
//      Dirichlet boundary conditions.
//
// You can change the constant right-hand side CONST_F, the
// initial polynomial degree P_INIT, and play with various initial
// mesh refinements at the beginning of the main() function.

double CONST_F = 2.0;   // Constant right-hand side.
int P_INIT = 5;         // Uniform polynomial degree of mesh elements.

// boundary condition types (essential = Dirichlet)
int bc_types(int marker)
{
	if (marker < SRC_MARKERS_START)
	  return BC_ESSENTIAL;
	else
		return BC_NATURAL;
}

// function values for Dirichlet boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return 0;
}

scalar rhs(double x, double y) {
  return exp(-x*x+y*y)/*sin(2*M_PI*x) + sin(2*M_PI*y)*/;
}

// return the value \int \nabla u . \nabla v dx
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

/*
// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F * int_v<Real, Scalar>(n, wt, v);
}
*/

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F * int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

// Integration order for the volumetric linear form
Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // sample "manual" mesh refinement
  mesh.refine_element(0);
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1PeriodicSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  BaseView bview;
  bview.show(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  //wf.add_liform(0, callback(linear_form));
  wf.add_liform(0, linear_form, linear_form_ord);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  // assemble the stiffness matrix and solve the system
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);

  // visualize the solution
  ScalarView view("Solution");
  view.show(&sln);
  
  // wait for a view to be closed
  View::wait();
  return 0;
}

