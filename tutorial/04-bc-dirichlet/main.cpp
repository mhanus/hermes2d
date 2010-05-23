#include "hermes2d.h"
#include "solver_umfpack.h"

// This example illustrates how to use nonhomogeneous (nonzero)
// Dirichlet boundary conditions.
//
// PDE: Poisson equation -Laplace u = CONST_F, where CONST_F is
// a constant right-hand side. It is not difficult to see that
// the function u(x,y) = (-CONST_F/4)*(x^2 + y^2) satisfies the
// above PDE. Since also the Dirichlet boundary conditions
// are chosen to match u(x,y), this function is the exact
// solution.
//
// Note that since the exact solution is a quadratic polynomial,
// Hermes will compute it exactly if all mesh elements are quadratic
// or higher (then the exact solution lies in the finite element space).
// If some elements in the mesh are linear, Hermes will only find
// an approximation, Below you can play with the parameters CONST_F,
// P_INIT, and UNIFORM_REF_LEVEL.

double CONST_F = -4.0;       // constant right-hand side
int P_INIT = 2;              // initial polynomial degree in all elements
int UNIFORM_REF_LEVEL = 3;   // number of initial uniform mesh refinements

// boundary condition type (essential = Dirichlet)
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
{
  return (-CONST_F/4.0)*(x*x + y*y);
}

// return the value \int \nabla u . \nabla v dx
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F*int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Refine mesh uniformly.
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache.
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create an H1 space.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate degrees of freedom.
  int ndof = assign_dofs(&space);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_biform(callback(bilinear_form));
  wf.add_liform(callback(linear_form));

  // Initialize the linear system and solver.
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_space(&space);
  sys.set_pss(&pss);

  // Assemble the stiffness matrix and solve the system.
  Solution sln;
  sys.assemble();
  sys.solve(&sln);

  // Visualize the solution.
  ScalarView view("Solution");
  view.show(&sln);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
