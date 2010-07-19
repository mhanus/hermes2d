// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

// $Id: view.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __H2D_STREAM_VIEW_H
#define __H2D_STREAM_VIEW_H

#include "view.h"

// you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

/// \brief Visualizes streamlines of a vector PDE solution.
///
/// StreamView is a visualization window for all vector-valued PDE solutions (especially for flow problems).
///
class H2D_API StreamView : public View
{
public:

  StreamView(const char* title = "StreamView", DEFAULT_WINDOW_POS);
  StreamView(const char* title = "StreamView", WinGeom* wg = NULL);
  StreamView(char* title, WinGeom* wg = NULL);
  virtual ~StreamView();

  /// Using velocity components (xsln, ysln) it creates streamlines that begin at the boundary with "marker"
  /// and the distance between starting points is "step"
  void show(MeshFunction* xsln, MeshFunction* ysln, int marker, double step, double eps = H2D_EPS_NORMAL);
  void show(MeshFunction* xsln, MeshFunction* ysln, int marker, double step, double eps, int xitem, int yitem);

  /// Creates additional streamline with strarting point (x,y)
  /// Note: Can be called only after StreamView::show
  void add_streamline(double x, double y);

protected:

  struct Node
  {
    bool leaf;
    int level;
    Node* sons[2];
    int elements[100];
    int num_elem;
  };

  Vectorizer vec;
  double max_mag;
  bool lines, pmode;

  double initial_tau;
  double min_tau;
  double max_tau;
  int num_stream;
  double2** streamlines;
  int* streamlength;
  Node* root;
  double root_x_min;
  double root_x_max;
  double root_y_min;
  double root_y_max;

  int find_triangle_in_tree(double x, double y, Node* father, double x_min, double x_max, double y_min, double y_max, double3& bar);
  void add_element_to_tree(Node* father, int e_idx, double x_min, double x_max, double y_min, double y_max);
  void build_tree();
  void delete_tree(Node* father);

  bool is_in_triangle(int idx, double x, double y, double3& bar);
  bool get_solution_values(double x, double y, double& xval, double& yval);

  int create_streamline(double x_start, double y_start, int idx);
  void find_initial_points(int marker, double step, double2*& initial_points);
  int find_initial_edge(int num_edges, int3* edges);

  virtual void on_display();
  virtual void on_mouse_move(int x, int y);
  virtual void on_key_down(unsigned char key, int x, int y);
  virtual void on_left_mouse_down(int x, int y);
  virtual const char* get_help_text() const;

};

#else // NOGLUT

/* Empty dummy implementation in a case GLUT is not used */
class H2D_API StreamView : public View
{
public:

  StreamView(const char* title = "StreamView", DEFAULT_WINDOW_POS) : View(title, x, y, width, height) {};
  virtual ~StreamView() {};

  void show(MeshFunction* xsln, MeshFunction* ysln, int marker, double step, double eps = H2D_EPS_NORMAL) {};
  void show(MeshFunction* xsln, MeshFunction* ysln, int marker, double step, double eps, int xitem, int yitem) {};

  void add_streamline(double x, double y) {};
};

#endif // NOGLUT

#endif
