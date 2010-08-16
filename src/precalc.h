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

#ifndef __H2D_PRECALC_H
#define __H2D_PRECALC_H

#include "function.h"
#include "shapeset.h"


/// \brief Caches precalculated shape function values.
///
/// PrecalcShapeset is a cache of precalculated shape function values.
///
///
class H2D_API PrecalcShapeset : public RealFunction
{
public:

  /// \brief Constructs a standard (master) precalculated shapeset class.
  /// \param shapeset [in] Pointer to the shapeset to be precalculated.
  PrecalcShapeset(Shapeset* shapeset);

  /// \brief Constructs a slave precalculated shapeset class.
  /// \details The slave instance does not hold any precalculated tables.
  /// Instead, it refers to those contained in the master instance. However,
  /// the slave can have different shape function active, different transform
  /// selected, etc. Slave pss's are used for test functions when calling
  /// bilinear forms, inside Solution so as not to disrupt user's pss, etc.
  /// \param master_pss [in] Master precalculated shapeset pointer.
  PrecalcShapeset(PrecalcShapeset* master_pss);

  /// Destructor.
  virtual ~PrecalcShapeset();

  virtual void set_quad_2d(Quad2D* quad_2d);

  /// \brief Frees all precalculated tables.
  virtual void free();

  /// Ensures subsequent calls to get_active_element() will be returning 'e'.
  /// Switches the class to the appropriate mode (triangle, quad).
  virtual void set_active_element(Element* e);

  /// Activates a shape function given by its index. The values of the shape function
  /// can then be obtained by setting the required integration rule order by calling
  /// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
  /// \param index [in] Shape index.
  void set_active_shape(int index);

  /// Returns the index of the active shape (can be negative if the shape is constrained).
  int get_active_shape() const { return index; };

  /// Returns a pointer to the shapeset which is being precalculated.
  Shapeset* get_shapeset() const { return shapeset; }

  /// Returns type of space (0 - H1, 1 - Hcurl, 2 - Hdiv)
  int get_type() const { return shapeset->get_id() / 10;}

  /// Internal. Use set_active_element() instead.
  void set_mode(int mode);

  /// For internal use only.
  void set_master_transform()
  {
    assert(master_pss != NULL);
    sub_idx = master_pss->sub_idx;
    top = master_pss->top;
    stack[top] = *(master_pss->ctm);
    ctm = stack + top;
  }

  void dump_info(int quad, const char* filename); // debug

  int get_edge_fn_order(int edge) { return make_edge_order(shapeset->get_order(index), edge); };
  
protected:

  Shapeset* shapeset;

  void* tables; ///< primary Judy array of shapes

  int mode;
  int index;
  int max_index[2];

  PrecalcShapeset* master_pss;

  bool is_slave() const { return master_pss != NULL; }

  virtual void precalculate(int order, int mask);

  void update_max_index();

  /// Forces a transform without using push_transform() etc.
  /// Used by the Solution class. <b>For internal use only</b>.
  void force_transform(uint64_t sub_idx, Trf* ctm)
  {
    this->sub_idx = sub_idx;
    this->ctm = ctm;
  }

  friend class Solution;
  friend class RefMap;

};



#endif
