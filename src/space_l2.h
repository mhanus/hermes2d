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

#ifndef __H2D_SPACE_L2
#define __H2D_SPACE_L2

#include "space.h"


/// L2Space represents a space of scalar functions with discontinuities along
/// mesh edges.
///
///
class H2D_API L2Space : public Space
{
public:

  L2Space(Mesh* mesh = NULL, int p_init = 1, Shapeset* shapeset = NULL);
  virtual ~L2Space();

  virtual Space* dup(Mesh* mesh) const;

  virtual int get_edge_order(Element* e, int edge) { return 0;}

  virtual int get_type() const { return 3; }

  virtual void get_element_assembly_list(Element* e, AsmList* al);
  
protected:

  struct L2Data
  {
    int vdof[4];
    int edof[4];
  };

  L2Data* ldata;
  int lsize;

  virtual void resize_tables();

  virtual void assign_vertex_dofs() {}
  virtual void assign_edge_dofs() {}
  virtual void assign_bubble_dofs();

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al) {}
  virtual void get_edge_assembly_list_internal(Element* e, int ie, AsmList* al);
  virtual void get_boundary_assembly_list(Element* e, int ie, AsmList* al);
  virtual void get_bubble_assembly_list(Element* e, AsmList* al);
  
  /// Project the boundary condition at given edge onto a space spanned by shape functions with non-zero values on that edge
  /// (these are in fact all the bubble functions spanning our discrete L2), using the edge-wise L2 inner product. 
  /// Returns coefficients of the linear combination of the shape functions that defines the projected boundary condition function.
  virtual scalar* get_bc_projection(EdgePos* ep, int order);
      
  /// Precalculates Cholesky decomposition of the projection matrix used for projecting boundary conditions.
  /// Currently, it does that for both the reference triangle and quad (see the comment in the function body).
  virtual void precalculate_projection_matrix(int mode);
  
  /// Arrays storing the Cholesky decomposition data (currently for both element types and every edge).
  static double** l2_proj_mat[2][4];
  static double*  l2_chol_p[2][4];
  
  /// Number of spaces that use the single precomputed Cholesky decomposition for quad / triangle.
  static int      l2_proj_ref[2]; 
};



#endif
