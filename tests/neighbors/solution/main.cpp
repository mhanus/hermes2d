#include "hermes2d.h"
#include "solver_umfpack.h"
#include <iostream>
using namespace std;
// This test is for testing if works in proper way transformation of solution on central or neighbor element and if
// afterward orientation of neighbors array of function values is done right.


#define H2D_ERROR_SUCCESS   0
#define H2D_ERROR_FAILURE  -1

const double NEIGHBOR_TRASHOLD = 1e-010;

bool equal_double(double value, double compare_value)
{
	if(value > (compare_value - NEIGHBOR_TRASHOLD) && value < (compare_value + NEIGHBOR_TRASHOLD))
		return true;
	else
		return false;
}


// The solution is function F bellow. Together with a mesh they were chosen to have wide range of function values,
// to have continuous solution over the mesh and to suppress possibility of getting pass of the test even though it wouldn't.

// projected function
double F(double x, double y)
{
//	return 1;
  return 2*(x - 0.5)*(x - 0.5) + 6*(y + 0.7)*(y + 0.7) - 36;
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
  for (int i = 0; i < n; i++)
    result += wt[i] * (( 2*pow(e->x[i] - 0.5, 2) + 6 * pow(e->y[i] + 0.7, 2) - 25) * v->val[i]);
  return result;
}
/*
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result;
	for(int i = 0; i < n; i++)
		result += ext->fn[0]->val[i] + ext->get_fn_neighbor(0)->val[i] ;
	return result+1;
}

*/

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
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

	// perform refinement to get complex mesh
	mesh.refine_element(0, 1);
	mesh.refine_element(2);
	mesh.refine_element(3, 1);
	mesh.refine_element(8, 2);
	mesh.refine_element(14, 2);
	mesh.refine_element(10);
	mesh.refine_element(16);
	mesh.refine_element(18, 1);
	mesh.refine_element(21);
	mesh.refine_element(24);
	mesh.refine_element(27, 1);
	mesh.refine_element(30);

  MeshView mview("Mesh ", 100, 100, 500, 500);
  mview.show(&mesh);

 // Create an H1 space.
 const int P_INIT = 1;
  H1Space space(&mesh, bc_types, NULL, P_INIT);

  // set uniform polynomial degrees
//	 space.set_uniform_order(P_INIT);
  Element* e = NULL;
  int order = 1;
	 for_all_active_elements(e, &mesh){
	   if ( e->is_quad() )
      space.set_element_order(e->id, H2D_MAKE_QUAD_ORDER(order, order));
     else
      space.set_element_order(e->id, order);
		 order = order < 9 ? order + 1 : 1;
	 }

	// enumerate basis functions
  Solution sln;
  Solution xprev;
  xprev.set_zero(&mesh);
  // matrix solver
  UmfpackSolver umfpack;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
/*
  wf.add_liform_surf(0, callback(linear_form_surf), ANY_INNER_EDGE, 1, &xprev);
*/

  // assemble and solve the finite element problem
  LinSystem sys(&wf, &umfpack, &space);
  sys.assemble();
  sys.solve(&sln);

  // visualize the solution
   ScalarView view1("Solution 1");
   view1.show(&sln);
   View::wait("Waiting for all views to be closed.");


  // begin of test
  NeighborSearch* neighb = NULL;
  int n_neighbors = 0;
  scalar* fn_central = NULL;
  scalar* fn_neighbor = NULL;
  int n_integ_points = 0;
  bool test = false;
	e = NULL;
  for_all_active_elements(e, &mesh)
  {
  	neighb = new NeighborSearch(e, &mesh, &sln, &space);
  	for(int i = 0; i < e->nvert; i++)
  	{
  		if(e->en[i]->bnd == 0){
    		neighb->set_active_edge(i);
  			n_neighbors = neighb->get_number_of_neighbs();
  			for(int j = 0; j < n_neighbors; j++)
  			{
  				fn_central = neighb->get_fn_values_central(j);
  				fn_neighbor = neighb->get_fn_values_neighbor(j);
					n_integ_points = neighb->get_n_integ_points(j);
					for(int k = 0; k < n_integ_points; k++)
					{
						test = equal_double(fn_central[k], fn_neighbor[k]);
						if(test == false)
						{
							for(int l = 0; l < n_integ_points; l++)
								debug_log("fn_central and fn_neighbor[%d]: %f, %f", l, fn_central[l], fn_neighbor[l]);
							printf("failure! \n");
							return H2D_ERROR_FAILURE;
						}
					}
  			}
  		}
  	}
  	delete neighb;
  }

  printf("success! \n");
  return H2D_ERROR_SUCCESS;
}
