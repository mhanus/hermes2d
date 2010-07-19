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

#ifndef __H2D_BASE_VIEW_H
#define __H2D_BASE_VIEW_H

#include "../common.h"
#include "scalar_view.h"

// you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

/// \brief Visualizes the basis functions of a space.
///
/// BaseView is a debugging tool for the visualization of the basis functions
/// of a given space.
///
class H2D_API BaseView : public ScalarView
{
public:

  BaseView(const char* title = "BaseView", DEFAULT_WINDOW_POS);
  BaseView(const char* title = "BaseView", WinGeom* wg = NULL);
  BaseView(char* title, WinGeom* wg = NULL);

  void show(Space* space, double eps = H2D_EPS_LOW, int item = H2D_FN_VAL_0);

  virtual ~BaseView() { free(); }

protected:

  Space* space;
  PrecalcShapeset* pss;
  Solution* sln;

  double eps;
  int ndof, item;
  int base_index;

  void free();
  void update_solution();
  void update_title();

  virtual void on_special_key(int key, int x, int y);
  virtual const char* get_help_text() const;

};

#else // NOGLUT

class H2D_API BaseView : public ScalarView
{
public:
  BaseView(const char* title = "BaseView", DEFAULT_WINDOW_POS) {};
  BaseView(const char* title = "BaseView", WinGeom* wg = NULL) {};
  BaseView(char* title, WinGeom* wg = NULL) {};
  virtual ~BaseView() {}
  void show(Space* space, double eps = H2D_EPS_LOW, int item = H2D_FN_VAL_0)
     { verbose("BaseView: Hermes2D compiled without OpenGL support, skipping visualization."); }
};

#endif

#endif
