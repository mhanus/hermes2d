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
     //DEPRECATED 
     //result += 0.5*wt[i] * (ext->fn[0]->val[i] + ext->get_fn_neighbor(0)->val[i]) * v->val[i];
     result += 0.5*wt[i] * (ext->fn[0]->get_val_central(i) + ext->fn[0]->get_val_neighbor(i)) * v->val[i];
  }
  return result;
}

// boundary conditions
BCType bc_types(int marker)
{
   return BC_ESSENTIAL;
}
// function values for Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return x*x*x + y*y*y;
}

int main(int argc, char* argv[])
{
  // Decide which mesh to refine and how from the command line argument.

  enum ref_type {NO_REF, DEFAULT_REF, TEST_REF};
  ref_type refinement = NO_REF;
  std::string mesh_file = "square.mesh";

  if (argc < 2) {
    info("Using default refinement of square.mesh.");
    refinement = DEFAULT_REF;
  }  else {
    if (!strcmp(argv[1], "test")) {
      info("Using testing refinement of square.mesh.");
      refinement = TEST_REF;
    }
    else {
      mesh_file = argv[1];
      info("Using mesh %s", argv[1]);
    }
  }

  // Load and refine the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(mesh_file.c_str(), &mesh);

  switch (refinement) {
    case DEFAULT_REF:
      // Initial uniform refinement.
      for (int i=0; i<3; i++) mesh.refine_all_elements();

      // Two additional refinements of selected elements.
      mesh.refine_element(29);
      mesh.refine_element(87);

      break;
    case TEST_REF:
      mesh.refine_element(0);
      mesh.refine_element(1);
      mesh.refine_element(7);
      break;
  }

  // create the L2 space
  L2Space space(&mesh,P_INIT);
  space.set_element_order(3, H2D_MAKE_QUAD_ORDER(1,4));
  space.set_element_order(4, H2D_MAKE_QUAD_ORDER(3,1));
  space.set_essential_bc_values(essential_bc_values);
  space.set_bc_types(bc_types); // Zde se vola assign_dofs(), kde nasledne update_edge_bc_values(); tam je zapotrebi mit essential_bc_values_by_coord, proto se set_bc_types musi volat az za set_essential_bc_values.

  // display the mesh
  OrderView oview("info_neighbor", 100, 100, 500, 500);
  oview.show(&space);
  BaseView bview;
  bview.show(&space);
  
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

