// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "space_l2.h"
#include "matrix_old.h"
#include "quad_all.h"
#include "shapeset_l2_all.h"


#include <iostream>
#include <iomanip>
using namespace std;

// Arrays that store precalculated Cholesky decomposition data for each element type and each edge.
double** L2Space::l2_proj_mat[2][4] = { { NULL, NULL, NULL, NULL }, { NULL, NULL, NULL, NULL } };
double*  L2Space::l2_chol_p[2][4]   = { { NULL, NULL, NULL, NULL }, { NULL, NULL, NULL, NULL } };
int      L2Space::l2_proj_ref[2] = { 0, 0 };

L2Space::L2Space(Mesh* mesh, int p_init, Shapeset* shapeset)
  : Space(mesh, shapeset, NULL, NULL, p_init)
{
  if (shapeset == NULL) this->shapeset = new L2Shapeset;
  ldata = NULL;
  lsize = 0;

  // Unlike the H1 edge shape functions, the L2 bubble shape functions used for boundary conditions projection
  // take different values on each edge of either a quadrilateral or a triangle, so we need to precalculate the projection
  // matrix using both shapesets and on all edges. As of now, existence of each edge of each element type on the boundary 
  // is not checked - with current maximum shape function order of 10, doing Cholesky decompositions for both quad and triangle
  // edges takes about (4+3)*222 add/mult's and (4+3)*11 sqrt's. 
  // FIXME: Should the order be radically increased, it would be advantageous to check if there are any triangles on the
  // boundary, even by going through all active boundary elements and testing their type (provided that a quadrilateral 
  // cannot change into triangle during adaptivity, or vice versa).
  
  if (!l2_proj_ref[H2D_MODE_TRIANGLE]++)
    precalculate_projection_matrix(H2D_MODE_TRIANGLE);
  if (!l2_proj_ref[H2D_MODE_QUAD]++)
    precalculate_projection_matrix(H2D_MODE_QUAD);
 
  // set uniform poly order in elements
  if (p_init < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}


L2Space::~L2Space()
{
  if (!--l2_proj_ref[H2D_MODE_TRIANGLE])
  {
    for (int edge = 0; edge < 3; edge++) {
      delete [] l2_proj_mat[H2D_MODE_TRIANGLE][edge];
      delete [] l2_chol_p[H2D_MODE_TRIANGLE][edge];
    }
  }
  if (!--l2_proj_ref[H2D_MODE_QUAD])
  {
    for (int edge = 0; edge < 4; edge++) {
      delete [] l2_proj_mat[H2D_MODE_QUAD];
      delete [] l2_chol_p[H2D_MODE_QUAD];
    }
  }
  
  ::free(ldata);
}


Space* L2Space::dup(Mesh* mesh) const
{
  L2Space* space = new L2Space(mesh, 0, shapeset);
  space->copy_callbacks(this);
  return space;
}


//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void L2Space::resize_tables()
{
  if (lsize < mesh->get_max_element_id())
  {
    if (!lsize) lsize = 1000;
    while (lsize < mesh->get_max_element_id()) lsize = lsize * 3 / 2;
    ldata = (L2Data*) realloc(ldata, sizeof(L2Data) * lsize);
  }
  Space::resize_tables();
}


void L2Space::assign_bubble_dofs()
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    shapeset->set_mode(e->get_mode());
    ElementData* ed = &edata[e->id];
    ed->bdof = next_dof;
    ed->n = shapeset->get_num_bubbles(ed->order); //FIXME: this function might return invalid value because retrieved bubble functions for non-uniform orders might be invalid for the given order.
    next_dof += ed->n * stride;
    
    // Now we also assign dof = -1 to all essential boundary edge nodes (which we identify by NodeData::n of their endpoints as set during Space::reset_dof_assignment),
    // as well as the number of b.c. projection coefficients.
    for (unsigned int ie = 0; ie < e->nvert; ie++) {
      if (ndata[e->vn[ie]->id].n == BC_ESSENTIAL && ndata[e->vn[e->next_vert(ie)]->id].n == BC_ESSENTIAL) {
        ndata[e->en[ie]->id].dof = H2D_CONSTRAINED_DOF;   // In L2, the only constraint may be on the boundary.
        ndata[e->en[ie]->id].n = ed->n;  // We need to project the essential b.c. on all bubble functions.
      }
    }
  }
}

  
  


//// assembly lists ////////////////////////////////////////////////////////////////////////////////

void L2Space::get_element_assembly_list(Element* e, AsmList* al)
{
  // some checks
  if (e->id >= esize || edata[e->id].order < 0)
    error("Uninitialized element order (id = #%d).", e->id);
  if (!is_up_to_date())
    error("The space is out of date. You need to update it with assign_dofs()"
          " any time the mesh changes.");

  // add vertex, edge and bubble functions to the assembly list
  al->clear();
  shapeset->set_mode(e->get_mode());

  get_bubble_assembly_list(e, al);
  
  // If the element lies on boundary, add the boundary condition projection coefficients.
  for (unsigned int ie = 0; ie < e->nvert; ie++)
    if (ndata[e->en[ie]->id].dof == H2D_CONSTRAINED_DOF)
      this->get_boundary_assembly_list(e, ie, al);
}

void L2Space::get_bubble_assembly_list(Element* e, AsmList* al)
{
  ElementData* ed = &edata[e->id];
  if (!ed->n) return;

  int* indices = shapeset->get_bubble_indices(ed->order);
  for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride) {
    //printf("triplet: %d, %d, %f\n", *indices, dof, 1.0);
    al->add_triplet(*indices++, dof, 1.0);
  }
}

void L2Space::get_edge_assembly_list_internal(Element* e, int ie, AsmList* al)
{
  this->get_bubble_assembly_list(e, al);
   
  // If the edge is the boundary edge, add the boundary condition projections
  if (ndata[e->en[ie]->id].dof == H2D_CONSTRAINED_DOF)
    this->get_boundary_assembly_list(e, ie, al);
}

void L2Space::get_boundary_assembly_list(Element* e, int ie, AsmList* al)
{
  Node* en = e->en[ie];
  NodeData* nd = &ndata[en->id];
  ElementData* ed = &edata[e->id];
  
  assert(nd->n == ed->n);
  if (!nd->n) return;
  
  int* indices = shapeset->get_bubble_indices(ed->order);
  for (int i = 0; i < nd->n; i++)
    al->add_triplet(*indices++, -1, nd->edge_bc_proj[i]);
}


scalar* L2Space::get_bc_projection(EdgePos* ep, int order)
{
  assert(order >= 1);  
  
  // Unlike in the H1 space we cannot distinguish the two spanning subsets (of vertex/edge shape functions) 
  // of the target space here, so we define the projection by orthogonalizing the difference between the projected
  // and the original boundary function with respect to all shape functions at once.

  // Allocate the output array to the number of shape functions.
  int mode = ep->base->get_mode();
  shapeset->set_mode(mode); 
  
  //if (mode == H2D_MODE_QUAD) order = H2D_MAKE_QUAD_ORDER(order, order);
  int* bubble_indices = shapeset->get_bubble_indices(order);
  int num_bubbles = shapeset->get_num_bubbles(order);
  scalar* proj = new scalar[num_bubbles];
  
  Quad2DStd quad2d;
  quad2d.set_mode(mode); 
  
  // Obtain integration points on an arbitrary edge using the maximal order quadrature rule (we cannot estimate regularity
  // of the boundary condition).
  int eo = quad2d.get_edge_points(ep->edge);
  int np = quad2d.get_num_points(eo);
  double3* pt = quad2d.get_points(eo);

  // Loop over all shape functions up to requested order and construct the right-hand side of the discrete projection problem.
  for (int i = 0; i < num_bubbles; i++)
  {
    int bubble_ii = bubble_indices[i];
    proj[i] = 0.0;
      
    for (int j = 0; j < np; j++)
    {      
      // Calculate the 1D parameter ep->t going between 0 and 1 within [ep->lo, ep->hi], corresponding to current i.p.
      // TODO: Check that it works with curves.
      double pos;
      if (mode == H2D_MODE_QUAD)
        pos = (ep->edge == 1  || ep->edge == 3) ? pt[j][1] : pt[j][0];
      else
        pos = (ep->edge == 2) ? pt[j][1] : pt[j][0];
    
      double t = (pos + 1) * 0.5, s = 1.0 - t;
      ep->t = ep->lo * s + ep->hi * t;
      
      proj[i] += pt[j][2] * shapeset->get_fn_value(bubble_ii, pt[j][0], pt[j][1], 0)
                          * bc_value_callback_by_edge(ep);
    }
  }

  // solve the system using a precalculated Cholesky decomposed projection matrix
  cholsl(l2_proj_mat[mode][ep->edge], order, l2_chol_p[mode][ep->edge], proj, proj);

  return proj;
}

void L2Space::precalculate_projection_matrix(int mode)
{ 
  shapeset->set_mode(mode);
  
  int order = shapeset->get_max_order();
  if (mode == H2D_MODE_QUAD) order = H2D_MAKE_QUAD_ORDER(order, order);
  int* bubble_indices = shapeset->get_bubble_indices(order);
  
  int num_bubbles = shapeset->get_num_bubbles(order);

  Quad2DStd quad2d;
  quad2d.set_mode(mode);
  for (int edge = 0; edge < (mode == H2D_MODE_QUAD) ? 4 : 3; edge++) 
  {
    l2_proj_mat[mode][edge] = new_matrix<double>(num_bubbles, num_bubbles);

    for (int i = 0; i < num_bubbles; i++)
    {
      int bubble_ii = bubble_indices[i];
      int oi = shapeset->get_order(bubble_ii);
      
      for (int j = i; j < num_bubbles; j++)
      {
        int bubble_ij = bubble_indices[j];  
        int oj = shapeset->get_order(bubble_ij);
        
        if (mode == H2D_MODE_QUAD)
          order = (edge == 1 || edge == 3) ? H2D_GET_V_ORDER(oi) + H2D_GET_V_ORDER(oj) 
                                           : H2D_GET_H_ORDER(oi) + H2D_GET_H_ORDER(oj);
        else
          order = oi + oj;
        
        //order = shapeset.get_order(bubble_ii) + shapeset.get_order(bubble_ij);
        //order = std::max(H2D_GET_H_ORDER(order), H2D_GET_V_ORDER(order));
        
        int eo = quad2d.get_edge_points(edge, order);
        double3* pt = quad2d.get_points(eo);
        int np = quad2d.get_num_points(eo);
        
        for (int k = 0; k < np; k++)
          l2_proj_mat[mode][edge][i][j] += pt[k][2] * shapeset->get_fn_value(bubble_ii, pt[k][0], pt[k][1], 0)
                                                    * shapeset->get_fn_value(bubble_ij, pt[k][0], pt[k][1], 0);
      }
    }
    
    double **a = l2_proj_mat[mode][edge];
    cout << endl;
    for (int i = 0; i < num_bubbles; i++) {
      for (int j = 0; j < num_bubbles; j++) {
        cout << setw(15) << a[i][j];
      }
      cout << endl;
    }
    cout << endl;
    
    l2_chol_p[mode][edge] = new double[num_bubbles];
    choldc(a, num_bubbles, l2_chol_p[mode][edge]);
  }
}
