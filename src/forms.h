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


#ifndef __H2D_FORMS_H
#define __H2D_FORMS_H

#include "common.h"
#include "quad.h"
#include "function.h"
#include "solution.h"
#include "refmap.h"

#define callback(a)	a<double, scalar>, a<Ord, Ord>

// Base type for orders of functions
//
// We defined a special arithmetics with this type to be able to analyze forms
// and determine the necessary integration order.  This works for forms, but it also
// works for user-defined functions.
class Ord
{
public:

  Ord(): order(0) {}
  explicit Ord(int o): order(o) {}
  Ord(double d): order(0) {}

  int get_order() const { return order; }
  int get_max_order() const {return 30;}

  Ord operator+(const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord operator+(double d) { return *this; }
  Ord operator-(const Ord &o) { return Ord(std::max(this->order, o.order)); }
  Ord operator-(double d) { return *this; }
  Ord operator*(const Ord &o) { return Ord(this->order + o.order); }
  Ord operator*(double d) { return *this; }
  Ord operator/(const Ord &o) { return Ord(this->get_max_order()); }
  Ord operator/(double d) { return *this; }

  Ord operator+=(const Ord &o) { this->order = std::max(this->order, o.order); return *this; }

  Ord operator+=(const double &d) { return *this; }
  Ord operator-=(const double &d) { return *this; }
  Ord operator*=(const double &d) { return *this; }
  Ord operator/=(const double &d) { return *this; }

  bool operator<(double d) { return true; }

protected:
  int order;

};

inline Ord operator/(const scalar &a, const Ord &b) { return Ord(b.get_max_order()); }
inline Ord operator*(const scalar &a, const Ord &b) { return b; }
inline Ord operator+(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const Ord &a) { return a; }

inline Ord pow(const Ord &a, const double &b) { return Ord((int) ceil(fabs(b)) * a.get_order()); }
inline Ord sqrt(const Ord &a) { return a; }
inline Ord sqr(const Ord &a) { return Ord(2 * a.get_order()); }
inline Ord conj(const Ord &a) { return a; }
inline Ord abs(const Ord &a) { return a; }

inline Ord atan2(const Ord &a, const Ord &b) { return Ord(a.get_max_order()); }
inline Ord atan(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord sin(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord cos(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord log(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord exp(const Ord &a) { return Ord(3 * a.get_order()); }

// Function
template<typename T>
class Func
{
protected:
  const int num_gip; ///< A number of integration points used by this intance.
  static const char* ERR_UNDEFINED_NEIGHBORING_ELEMENTS;
  
public:
  const int nc;	///< A number of components. Currently accepted values are 1 (H1, L2 space) and 2 (Hcurl, Hdiv space).
  T *val;					// function values. If orders differ for a diffrent
                                                // direction, this returns max(h_order, v_order).
  T *dx, *dy; 					// derivatives
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  T *laplace;                                   // must be enabled by defining H2D_SECOND_DERIVATIVES_ENABLED
                                                // in common.h. Default is NOT ENABLED.
#endif
  T *val0, *val1;				// components of function values
  T *dx0, *dx1;					// components of derivatives
  T *dy0, *dy1;

  T *curl;					 // components of curl

  /// Constructor.
  /** \param[in] num_gip A number of integration points.
   *  \param[in] num_comps A number of components. */
  explicit Func(int num_gip, int num_comps) : num_gip(num_gip), nc(num_comps) {
    val = val0 = val1 = NULL;
    dx = dx0 = dx1 = NULL;
    dy = dy0 = dy1 = NULL;
    curl = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    laplace = NULL;
#endif
  };
  
/// Subtract arrays stored in a given attribute from the same array in provided function.
#define H2D_SUBTRACT_IF_NOT_NULL(__ATTRIB, __OTHER_FUNC) { if (__ATTRIB != NULL) { \
  assert_msg(__OTHER_FUNC.__ATTRIB != NULL, "Unable to subtract a function expansion " #__ATTRIB " is NULL in the other function."); \
  for(int i = 0; i < num_gip; i++) __ATTRIB[i] -= __OTHER_FUNC.__ATTRIB[i]; } }

  /// Calculate this -= func for each function expations and each integration point.
  /** \param[in] func A function which is subtracted from *this. A number of integratio points and a number of component has to match. */
  virtual void subtract(const Func<T>& func) {
    assert_msg(num_gip == func.num_gip, "Unable to subtract a function due to a different number of integration points (this: %d, other: %d)", num_gip, func.num_gip);
    assert_msg(nc == func.nc, "Unable to subtract a function due to a different number of components (this: %d, other: %d)", nc, func.nc);
    H2D_SUBTRACT_IF_NOT_NULL(val, func)
    H2D_SUBTRACT_IF_NOT_NULL(dx, func)
    H2D_SUBTRACT_IF_NOT_NULL(dy, func)
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    H2D_SUBTRACT_IF_NOT_NULL(laplace, func)
#endif
    if (nc > 1) {
      H2D_SUBTRACT_IF_NOT_NULL(val0, func)
      H2D_SUBTRACT_IF_NOT_NULL(val1, func)
      H2D_SUBTRACT_IF_NOT_NULL(dx0, func)
      H2D_SUBTRACT_IF_NOT_NULL(dx1, func)
      H2D_SUBTRACT_IF_NOT_NULL(dy0, func)
      H2D_SUBTRACT_IF_NOT_NULL(dy1, func)
      H2D_SUBTRACT_IF_NOT_NULL(curl, func)
    }
  };
#undef H2D_SUBTRACT_IF_NOT_NULL

  virtual void free_ord() {
    delete val;
    val = val0 = val1 = NULL;
    dx = dx0 = dx1 = NULL;
    dy = dy0 = dy1 = NULL;
    curl = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    laplace = NULL;
#endif
  }
  virtual void free_fn()
  {
    delete [] val; val = NULL;
    delete [] dx; dx = NULL;
    delete [] dy; dy = NULL;
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    delete [] laplace; laplace = NULL;
#endif

    delete [] val0; delete [] val1; val0 = val1 = NULL;
    delete [] dx0;  delete [] dx1; dx0 = dx1 = NULL;
    delete [] dy0;  delete [] dy1; dy0 = dy1 = NULL;
    delete [] curl; curl = NULL;
  }
  
  virtual void ~Func() { } // All deallocation done via free_fn / free_ord. 
                           // This is to allow proper destruction of DiscontinuousFunc by applying delete on a Func pointer.
  
  virtual T& get_val_central(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); } 
  virtual T& get_val_neighbor(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); }
  virtual T& get_dx_central(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); }
  virtual T& get_dx_neighbor(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); }
  virtual T& get_dy_central(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS);  }
  virtual T& get_dy_neighbor(int k) { error(ERR_UNDEFINED_NEIGHBORING_ELEMENTS); }
};

#include "neighbor.h"

#define GET_CENT(__ATTRIB) return (fn_central != NULL ? fn_central->__ATTRIB[k] : zero;
#define GET_NEIB(__ATTRIB) return (fn_neighbor != NULL ? fn_neighbor->__ATTRIB[ reverse_neighbor_side ? num_gip-k-1 : k ] : zero;
template<typename T>
class DiscontinuousFunc : public Func<T>
{
  private:
    Func<T> *fn_central, *fn_neighbor;
    bool reverse_neighbor_side;
    static const T zero;
    
  public: 
    DiscontinuousFunc(const Func<T>* fn, bool support_on_neighbor = false) : Func(fn->num_gip, fn->nc), 
                                                                             fn_central(NULL), fn_neighbor(NULL),
                                                                             reverse_neighbor_side(false) { 
      assert_msg(fn != NULL, "Invalid arguments to DiscontinuousFunc constructor.");
      if (support_on_neighbor) fn_neighbor = fn; else fn_central = fn;
    }
    
    DiscontinuousFunc(const Func<T>* fn_c, const Func<T>* fn_n, bool reverse = false) : Func(fn_c->num_gip, fn_c->nc) { 
      assert_msg(fn_c != NULL && fn_n != NULL, "Invalid arguments to DiscontinuousFunc constructor.");
      assert_msg(fn_c->num_gip != fn_n->num_gip || fn_c->nc != fn_n->nc,
                 "DiscontinuousFunc must be formed by two Func's with same number of integration points and components.");
      fn_central = fn_c;
      fn_neighbor = fn_n;
      reverse_neighbor_side = reverse;
    }
    
    T& get_val_central(int k) { GET_CENT(val); }
    T& get_val_neighbor(int k) { GET_NEIB(val); }
    T& get_dx_central(int k) { GET_CENT(dx); }
    T& get_dx_neighbor(int k) { GET_NEIB(dx); }
    T& get_dy_central(int k) { GET_CENT(dy); }
    T& get_dy_neighbor(int k) { GET_NEIB(dy); }
    
    #ifdef H2D_SECOND_DERIVATIVES_ENABLED
      T& get_laplace_central(int k) { GET_CENT(laplace); }
      T& get_laplace_neighbor(int k) { GET_NEIB(laplace); }
    #endif
    
    void subtract(const Func<T>& func) { if (fn_central != NULL) fn_central->subtract(func); else fn_neighbor->subtract(func); }
    void subtract(const DiscontinuousFunc<T>& func) { /*TODO*/ }
    
    // Default destructor may be used. Deallocation is done using the following functions.
    // FIXME: This is not safe since it allows calling free_ord in a Func<scalar> object. Template-specialized destructors should be used instead (also in Func).
    void free_fn() { 
      if (fn_central != NULL) {
        fn_central->free_fn(); 
        fn_central = NULL;
      }
      if (fn_neighbor != NULL) {
        fn_neighbor->free_fn();
        fn_neighbor = NULL;
      }
    }
    void free_ord() {
      if (fn_central != NULL) {
        fn_central->free_ord(); 
        fn_central = NULL;
      }
      if (fn_neighbor != NULL) {
        fn_neighbor->free_ord();
        fn_neighbor = NULL;
      }
    }
};


/// Geometry (coordinates, normals, tangents) of either an element or an edge
template<typename T>
class Geom
{
public:
  int marker;       // marker
  int id;
  T diam;           // element diameter
  //Element *element;   // active element. NOTE: We used this for some time but
                        // decided against it because (a) it disables automatic
                        // order parsing and (b) if the form is called with T
                        // == Ord, element is not initialized, so the user has
                        // to be aware of this and test it in his weak form.

  T *x, *y;         // coordinates [in physical domain]
  T *nx, *ny;       // normals [in physical domain] (locally oriented
                    // to point outside the element)
  T *tx, *ty;       // tangents [in physical domain]
  int orientation;  // 0 .... if (nx, ny) is equal to the global normal,
                    // otherwise 1 (each mesh edge has a unique global normal
                    // vector)

  Geom()
  {
    marker = 0;
    id = 0;
    x = y = NULL;
    nx = ny = NULL;
    tx = ty = NULL;
    diam = 0;
  }

  void free()
  {
    delete [] tx;    delete [] ty;
    delete [] nx;    delete [] ny;
  }
};

/// Init element geometry for calculating the integration order
Geom<Ord>* init_geom_ord();
/// Init element geometry for volumetric integrals
Geom<double>* init_geom_vol(RefMap *rm, const int order);
/// Init element geometry for surface integrals
Geom<double>* init_geom_surf(RefMap *rm, EdgePos* ep, const int order);


/// Init the function for calculation the integration order
Func<Ord>* init_fn_ord(const int order);
/// Init the shape function for the evaluation of the volumetric/surface integral (transformation of values)
Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
/// Init the mesh-function for the evaluation of the volumetric/surface integral
Func<scalar>* init_fn(MeshFunction *fu, RefMap *rm, const int order);


/// User defined data that can go to the bilinear and linear forms.
/// It also holds arbitraty number of functions, that user can use.
/// Typically, these functions are solutions from the previous time/iteration levels.
template<typename T>
class ExtData {
public:
	int nf;			  			// number of functions in 'fn' array
	Func<T>** fn;				// array of pointers to functions

	ExtData() {
		nf = 0;
		fn = NULL;
		nf_neighbor = 0;
		fn_neighbor = NULL;
	}

  void free()
  {
    for (int i = 0; i < nf; i++)
    {
      fn[i]->free_fn();
      delete fn[i];
    }
    delete [] fn;
  }

  void free_neighbor()
  {
    for (int i = 0; i < nf_neighbor; i++)
    {
      fn_neighbor[i]->free_fn();
      delete fn_neighbor[i];
    }
    delete [] fn_neighbor;
  }


  void free_ord()
  {
    for (int i = 0; i < nf; i++)
    {
      fn[i]->free_ord();
      delete fn[i];
    }
    delete [] fn;
  }

  void free_neighbor_ord()
  {
    for (int i = 0; i < nf_neighbor; i++)
    {
      fn_neighbor[i]->free_ord();
      delete fn_neighbor[i];
    }
    delete [] fn_neighbor;
  }


  // Used for getting neighbor function values in forms (now only in linear surface form)
  int get_nf_neighbor()
  {
  	return nf_neighbor;
  }

  Func<T>* get_fn_neighbor(int index)
  {
  	return fn_neighbor[index];
  }

	void set_nf_neighbor(int n){ nf_neighbor = n;}

	void set_fn_neighbor(Func<T>** functions){fn_neighbor = functions;}
private:
	int nf_neighbor;			  			// number of functions in 'fn_neighbor' array
	Func<T>** fn_neighbor;				// array of pointers to functions on neighbor element

};

#endif
