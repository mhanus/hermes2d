// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@utep.edu>
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

// $Id: adapt_ortho_h1.cpp 1095 2008-10-24 13:19:45Z jakub $

#include "common.h"
#include "solution.h"
#include "linsystem.h"
#include "refmap.h"
#include "shapeset_h1_all.h"
#include "quad_all.h"
#include "integrals_h1.h"
#include "matrix.h"
#include "adapt_ortho_h1.h"
#include "traverse.h"
#include "norm.h"


H1OrthoHP::H1OrthoHP(int num, ...)
{
  this->num = num;
  
  va_list ap;
  va_start(ap, num);
  for (int i = 0; i < num; i++)
    spaces[i] = va_arg(ap, Space*);
  va_end(ap);


  memset(errors, 0, sizeof(errors));
  esort = NULL;
  have_errors = false;
}


H1OrthoHP::~H1OrthoHP()
{
  for (int i = 0; i < num; i++)
    if (errors[i] != NULL) 
      delete [] errors[i];
    
  if (esort != NULL) delete [] esort;
}


//// orthonormal base construction /////////////////////////////////////////////////////////////////

double3** H1OrthoHP::obase[2][9];
int  H1OrthoHP::basecnt[2][11];
bool H1OrthoHP::obase_ready = false;


void H1OrthoHP::calc_ortho_base()
{
  int i, j, k, l, m, ii, nb, np, o, r;
  int n, idx[121];

  H1Shapeset shapeset;

  // allocate the orthonormal base tables - these are simply the values of the
  // orthonormal functions in integration points; we store the basic functions
  // plus four son cut-outs of them (i.e. 5 times)
  for (i = 0; i < 9; i++)
  {
    if ((i < 4) || (i >= 8)) 
      obase[0][i] = new_matrix<double3>(66, 79); // tri
    obase[1][i] = new_matrix<double3>(121, 121); // quad
  }
  
  // repeat for triangles and quads
  for (m = 0; m <= 1; m++)
  {
    shapeset.set_mode(m);

    // obtain a list of all shape functions up to the order 10, from lowest to highest order
    n = 0;
    int nv = m ? 4 : 3;
    int num_sons = m ? 8 : 4;
    for (i = 0; i < nv; i++)
      idx[n++] = shapeset.get_vertex_index(i);
    basecnt[m][0] = 0;
    basecnt[m][1] = n;

    for (i = 2; i <= 10; i++)
    {
      for (j = 0; j < nv; j++)
        idx[n++] = shapeset.get_edge_index(j, 0, i);

      ii = m ? make_quad_order(i, i) : i;
      nb = shapeset.get_num_bubbles(ii);
      int* bub = shapeset.get_bubble_indices(ii);
      for (j = 0; j < nb; j++)
      {
        o = shapeset.get_order(bub[j]);
        if (get_h_order(o) == i || get_v_order(o) == i)
          idx[n++] = bub[j];
      }
      basecnt[m][i] = n;
    }

    // obtain their values for integration rule 20
    g_quad_2d_std.set_mode(m);
    np = g_quad_2d_std.get_num_points(20);
    double3* pt = g_quad_2d_std.get_points(20);

    for (i = 0; i < n; i++)
      for (j = 0; j < np; j++)
        for (k = 0; k < 3; k++)
          obase[m][8][i][j][k] = shapeset.get_value(k, idx[i], pt[j][0], pt[j][1], 0);
  
    for (l = 0; l < num_sons; l++)
    {
      Trf* tr = (m ? quad_trf : tri_trf) + l;
      for (i = 0; i < n; i++)
        for (j = 0; j < np; j++)
        {
          double x = tr->m[0]*pt[j][0] + tr->t[0],
                 y = tr->m[1]*pt[j][1] + tr->t[1];
          for (k = 0; k < 3; k++)
            obase[m][l][i][j][k] = shapeset.get_value(k, idx[i], x, y, 0);
        }
    }

    // orthonormalize the basis functions
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
      {
        double prod = 0.0;
        for (k = 0; k < np; k++) {
          double sum = 0.0;
          for (r = 0; r < 3; r++)
            sum += obase[m][8][i][k][r] * obase[m][8][j][k][r];
          prod += pt[k][2] * sum;
        }
        
        for (l = 0; l < 9; l++)
          if (m || l < 4 || l >= 8)
            for (k = 0; k < np; k++)
              for (r = 0; r < 3; r++)
                obase[m][l][i][k][r] -= prod * obase[m][l][j][k][r];
      }

      double norm = 0.0;
      for (k = 0; k < np; k++) {
        double sum = 0.0;
        for (r = 0; r < 3; r++)
          sum += sqr(obase[m][8][i][k][r]);
        norm += pt[k][2] * sum;
      }
      norm = sqrt(norm);

      for (l = 0; l < 9; l++)
        if (m || l < 4 || l >= 8)
          for (k = 0; k < np; k++)
            for (r = 0; r < 3; r++)
              obase[m][l][i][k][r] /= norm;
    }
    
    // check the orthonormal base
/*    if (m) {
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
        double check = 0.0;
        for(int son = 4; son < 6; son++ )
          for (k = 0; k < np; k++)
            check += pt[k][2] * (obase[m][son][i][k][0] * obase[m][son][j][k][0] + 
                                 obase[m][son][i][k][1] * obase[m][son][j][k][1] +
                                 obase[m][son][i][k][2] * obase[m][son][j][k][2]);
        check *= 0.5;
        if ((i == j && fabs(check - 1.0) > 1e-8) || (i != j && fabs(check) > 1e-8))
          warn("Not orthonormal: base %d times base %d = %g", i, j , check);  
      }
    }*/
  }
  obase_ready = true;
}


void H1OrthoHP::free_ortho_base()
{
  if (!obase_ready) return;
    
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 2; j++)
      delete [] obase[j][i];
    
  obase_ready = false;
}


//// optimal refinement search /////////////////////////////////////////////////////////////////////

void H1OrthoHP::calc_projection_errors(Element* e, int order, Solution* rsln,
                                       double herr[8][11], double perr[11])
{
  int i, j, s, k, r, son;
  int m = e->get_mode();
  double error;
  scalar prod;
  
  if (!obase_ready) calc_ortho_base();
    
  // select quadrature, obtain integration points and weights
  Quad2D* quad = &g_quad_2d_std;
  quad->set_mode(m);
  rsln->set_quad_2d(quad);
  double3* pt = quad->get_points(20);
  int np = quad->get_num_points(20);

  // everything is done on the reference domain
  // -- no reference mapping, no transformations
  rsln->enable_transform(false);
  
  // obtain reference solution values on all four refined sons
  scalar* rval[4][3];
  Element* base = rsln->get_mesh()->get_element(e->id);
  assert(!base->active);
  for (son = 0; son < 4; son++)
  {
    Element* e = base->sons[son];
    assert(e != NULL);
    rsln->set_active_element(e);
    rsln->set_quad_order(20);
    rval[son][0] = rsln->get_fn_values();
    rval[son][1] = rsln->get_dx_values();
    rval[son][2] = rsln->get_dy_values();
  }

  // h-cadidates: calculate products of the reference solution with orthonormal basis
  // functions on son elements, obtaining (partial) projections and their errors
  scalar3 proj[4][121];
  for (son = 0; son < 4; son++)
  {
    memset(proj[0], 0, sizeof(proj[0]));
    for (i = 1; i <= order; i++)
    {
      // update the projection to the current order
      for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
      {
        for (k = 0, prod = 0.0; k < np; k++)
          prod += pt[k][2] * (rval[son][0][k] * obase[m][8][j][k][0] +
                              rval[son][1][k] * obase[m][8][j][k][1] +
                              rval[son][2][k] * obase[m][8][j][k][2]);  
                                                                        
        for (k = 0; k < np; k++)                                        
          for (r = 0; r < 3; r++)
            proj[0][k][r] += obase[m][8][j][k][r] * prod;
      }

      // calculate the error of the projection
      for (k = 0, error = 0.0; k < np; k++)
        error += pt[k][2] * (sqr(rval[son][0][k] - proj[0][k][0]) +
                             sqr(rval[son][1][k] - proj[0][k][1]) +
                             sqr(rval[son][2][k] - proj[0][k][2]));
      herr[son][i] = error;
    }
  }

  // aniso-candidates: calculate projections and their errors (only quadrilaterals)
  if (m)
  {
    const double mx[4] = { 2.0, 2.0, 1.0, 1.0};
    const double my[4] = { 1.0, 1.0, 2.0, 2.0};
    const int sons[4][2] = { {0,1}, {3,2}, {0,3}, {1,2} };
    const int tr[4][2]   = { {6,7}, {6,7}, {4,5}, {4,5} };
  
    for (son = 0; son < 4; son++) // 2 sons for vertical split, 2 sons for horizontal split
    {
      memset(proj, 0, sizeof(proj));
      for (i = 1; i <= order+1; i++)  // h-candidates: max order equals to original element order+1
      {
        // update the projection to the current order
        for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
        {
          for (s = 0, prod = 0.0; s < 2; s++) // each son has 2 subsons (regular square sons)
            for (k = 0; k < np; k++)
              prod += pt[k][2] * ( rval[sons[son][s]][0][k]           * obase[m][tr[son][s]][j][k][0] +
                                   rval[sons[son][s]][1][k] * mx[son] * obase[m][tr[son][s]][j][k][1] +
                                   rval[sons[son][s]][2][k] * my[son] * obase[m][tr[son][s]][j][k][2]);
          prod *= 0.5;
    
          for (s = 0; s < 2; s++)
            for (k = 0; k < np; k++)
              for (r = 0; r < 3; r++)
                proj[s][k][r] += prod * obase[m][tr[son][s]][j][k][r];
        }
      
        // calculate the error of the projection
        for (s = 0, error = 0.0; s < 2; s++)
          for (k = 0; k < np; k++)
            error += pt[k][2] * (sqr(rval[sons[son][s]][0][k]           - proj[s][k][0]) +
                                 sqr(rval[sons[son][s]][1][k] * mx[son] - proj[s][k][1]) + 
                                 sqr(rval[sons[son][s]][2][k] * my[son] - proj[s][k][2]));
        herr[4 + son][i] = error * 0.5;
      }
    }
  }

  // p-candidates: calculate projections and their errors
  memset(proj, 0, sizeof(proj));
  for (i = 1; i <= std::min(order+2, 10); i++)
  {
    // update the projection to the current order
    for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
    {
      for (son = 0, prod = 0.0; son < 4; son++)
      {
        // (transforming to the quarter of the reference element)
        double mm = (e->is_triangle() && son == 3) ? -2.0 : 2.0;  
        
        for (k = 0; k < np; k++)
        {
          prod += pt[k][2] * (rval[son][0][k] *      obase[m][son][j][k][0] +
                              rval[son][1][k] * mm * obase[m][son][j][k][1] +
                              rval[son][2][k] * mm * obase[m][son][j][k][2]);
        }
      }
      prod *= 0.25;

      for (son = 0; son < 4; son++)
        for (k = 0; k < np; k++)
          for (r = 0; r < 3; r++)
            proj[son][k][r] += prod * obase[m][son][j][k][r];
    }

    // calculate the error of the projection
    for (son = 0, error = 0.0; son < 4; son++)
    {
      double mm = (e->is_triangle() && son == 3) ? -2.0 : 2.0; 
      
      for (k = 0; k < np; k++)
        error += pt[k][2] * (sqr(rval[son][0][k]      - proj[son][k][0]) +
                             sqr(rval[son][1][k] * mm - proj[son][k][1]) + 
                             sqr(rval[son][2][k] * mm - proj[son][k][2]));
    }
    perr[i] = error * 0.25;
  }
}


void H1OrthoHP::get_optimal_refinement(Element* e, int order, Solution* rsln,
                                       int& split, int p[4], int q[4], bool aniso, bool h_adapt)
{
  int i, j, k, n = 0;
  const int maxcand = 300;
  
  order = std::max(get_h_order(order), get_v_order(order));
  bool tri = e->is_triangle();

  // calculate maximal order of elements
  // linear elements = 9
  // curvilinear elements = depends on iro_cache (how curved they are)
  int max_order = (20 - e->iro_cache)/2 - 1;

 
  
  struct Cand
  {
    double error;
    int dofs, split, p[4];
  };
  Cand cand[maxcand];
  
  #define make_p_cand(q) { \
    assert(n < maxcand);   \
    cand[n].split = -1; \
    cand[n].p[1] = cand[n].p[2] = cand[n].p[3] = 0; \
    cand[n++].p[0] = (q); }
  
  #define make_hp_cand(q0, q1, q2, q3) { \
    assert(n < maxcand);  \
    cand[n].split = 0; \
    cand[n].p[0] = (q0); \
    cand[n].p[1] = (q1); \
    cand[n].p[2] = (q2); \
    cand[n++].p[3] = (q3); }
    
  #define make_ani_cand(q0, q1, iso) { \
    assert(n < maxcand);  \
    cand[n].split = iso; \
    cand[n].p[2] = cand[n].p[3] = 0; \
    cand[n].p[0] = (q0); \
    cand[n++].p[1] = (q1); }

  if (h_adapt)
  {
    make_p_cand(order);
    make_hp_cand(order, order, order, order);
    if ((!tri) && (e->iro_cache < 8) && aniso) {
      make_ani_cand(order, order, 1);
      make_ani_cand(order, order, 2);
    }
  }
  else {
    // prepare p-candidates
    int p0, p1 = std::min(max_order, order+1);
    for (p0 = order; p0 <= p1; p0++)
      make_p_cand(p0);
    
    //prepare hp-candidates
    p0 = (order+1) / 2;
    p1 = std::min(p0 + 3, order);
    int q0, q1, q2, q3;
    for (q0 = p0; q0 <= p1; q0++)
      for (q1 = p0; q1 <= p1; q1++)
        for (q2 = p0; q2 <= p1; q2++)
          for (q3 = p0; q3 <= p1; q3++)
            make_hp_cand(q0, q1, q2, q3);
  
    //prepare anisotropic candidates
    //only for quadrilaterals
    //too distorted (curved) elements cannot have aniso refinement (produces even worse elements)
    if ((!tri) && (e->iro_cache < 8) && aniso) {
      p0 = 2 * (order+1) / 3;
      int p_max = std::min(max_order, order+1);
      p1 = std::min(p0 + 3, p_max);    
      for (q0 = p0; q0 <= p1; q0++)
        for (q1 = p0; q1 <= p1; q1++) {
          if ((q0 < order+1) || (q1 < order+1)) {      
            make_ani_cand(q0, q1, 1);
            make_ani_cand(q0, q1, 2);
          }
        }
    }
  }
  // calculate (partial) projection errors
  double herr[8][11], perr[11];
  calc_projection_errors(e, order, rsln, herr, perr);
        
  // evaluate candidates (sum partial projection errors, calculate dofs)
  double avg = 0.0;
  double dev = 0.0;
  for (i = k = 0; i < n; i++)
  {
    Cand* c = cand + i;
    if (c->split == 0)
    {
      c->error = 0.0;
      c->dofs = tri ? 6 : 9;
      for (j = 0; j < 4; j++)
      {
        int o = c->p[j];
        c->error += herr[j][o];// * 0.25; // spravny vypocet chyby
        if (tri) {
          c->dofs += (o-2)*(o-1)/2;
          if (j < 3) c->dofs += std::min(o, c->p[3])-1 + 2*(o-1);
        }
        else {
          c->dofs += sqr(o)-1;
          c->dofs += /*2 * */std::min(o, c->p[j>0 ? j-1 : 3]) - 1;
        }
      }
    }
    else if (c->split == 1 || c->split == 2)  // aniso splits
    {
      c->dofs  = 6 /* vertex */ + 3*(c->p[0] - 1 + c->p[1] - 1); // edge fns
      c->dofs += std::min(c->p[0], c->p[1]) - 1; // common edge
      c->dofs += sqr(c->p[0] - 1) + sqr(c->p[1] - 1); // bubbles
      for (c->error = 0.0, j = 0; j < 2; j++)
        c->error += herr[(c->split == 1) ? j+4 : j+6][c->p[j]];// * 0.5;  // spravny vypocet chyby

    }
    else
    {
      int o = c->p[0];
      c->error = perr[o];
      c->dofs  = tri ? (o+1)*(o+2)/2 : sqr(o+1);
    }
    c->error = sqrt(c->error);

    if (!i || c->error <= cand[0].error)
    { 
      avg += log10(c->error); 
      dev += sqr(log10(c->error));
      k++; 
    }
  }
  avg /= k;  // mean
  dev /= k;  // second moment
  dev = sqrt(dev - sqr(avg));  // deviation is square root of variance   

  // select an above-average candidate with the steepest error decrease
  int imax = 0, h_imax = 0;
  double score, maxscore = 0.0, h_maxscore = 0.0;
  for (i = 1; i < n; i++)
  {
    if ((log10(cand[i].error) < avg + dev) && (cand[i].dofs > cand[0].dofs))
    {
      score = (log10(cand[0].error) - log10(cand[i].error)) / (cand[i].dofs - cand[0].dofs);
      if (score > maxscore) { maxscore = score; imax = i; }
      if ((cand[i].split == 0) && (score > h_maxscore)) { h_maxscore = score; h_imax = i; }
    }
  }

  // return result
  split = cand[imax].split;
  memcpy(p, cand[imax].p, 4*sizeof(int));
  memcpy(q, cand[h_imax].p, 4*sizeof(int));

}


//// adapt /////////////////////////////////////////////////////////////////////////////////////////

void H1OrthoHP::adapt(double thr, bool h_only, int strat)
{
  
  if (!have_errors)
    error("Element errors have to be calculated first, see calc_error().");

  bool aniso = true;
  
  int i, j, l;
  int max_id = -1;  
  Mesh* mesh[10]; 
  for (j = 0; j < num; j++) {
    mesh[j] = spaces[j]->get_mesh();
    rsln[j]->set_quad_2d(&g_quad_2d_std);  
    rsln[j]->enable_transform(false); 
    if (mesh[j]->get_max_element_id() > max_id)
      max_id = mesh[j]->get_max_element_id();
  }


  int split[nact];
  memset(split, 0, nact*sizeof(int));
  int p[nact][4], q[nact][4];
  int idx[max_id + 1][num + 1];
  for(j = 0; j < max_id; j++)
    for(l = 0; l < num; l++)
      idx[j][l] = -1; // element not refined

  int nref = nact;
  double err0 = 1000.0;
  double processed_error = 0.0;

  for (i = 0; i < nact; i++)
  {
    int comp = esort[i][1];
    int id = esort[i][0];
    double err = errors[comp][id];

    // first refinement strategy:
    // refine elements until prescribed amount of error is processed
    // if more elements have similar error refine all to keep the mesh symmetric
    if ((strat == 0) && (processed_error > sqrt(thr) * total_err) && fabs((err - err0)/err0) > 1e-3) { nref = i; break; }

    // second refinement strategy:
    // refine all elements whose error is bigger than some portion of maximal error
    if ((strat == 1) && (err < thr * errors[esort[0][1]][esort[0][0]])) { nref = i; break; }

    
    Element* e;
    e = mesh[comp]->get_element(id);
    int current = spaces[comp]->get_element_order(id);
    
    if (!h_only || aniso)
      get_optimal_refinement(e, current, rsln[comp], split[i], p[i], q[i], aniso, h_only);
    else
    {
      p[i][0] = p[i][1] = p[i][2] = p[i][3] = current;
      q[i][0] = q[i][1] = q[i][2] = q[i][3] = current;
    }

    idx[id][comp] = i;
    err0 = err;
    processed_error += err; 
  }

  int k = nref;
  for (i = 0; i < nref; i++)
  {
    int comp = esort[i][1];
    int id = esort[i][0];
    int current = get_h_order(spaces[comp]->get_element_order(id));

    int max_ref = split[i]; // how all the other elements will be refined
    for (j = 0; j < num; j++)
    {
      if (max_ref == 0) break; // iso refinement is max what can be recieved
      if ((j != comp) && (mesh[j] == mesh[comp])) // components share the mesh
      {
        int ii = idx[id][j];
        if ((ii >= 0) && (split[ii] != max_ref) && (split[ii] >= 0)) // ii element refined, refinement differs from max_ref, ii element split
        {
          if (((split[ii] == 1) || (split[ii] == 2)) && (max_ref == -1)) // the only case when aniso refinement
            max_ref = split[ii];
          else // otherwise isotropic refinement
            max_ref = 0;
        }
      }
    }    

    if (max_ref >= 0) 
    {
      for (j = 0; j < num; j++)
      {
        if ((j != comp) && (mesh[j] == mesh[comp])) // components share the mesh
        {
          // change appropriately original element
          if (split[i] != max_ref)
          {
            split[i] = max_ref;
            if (split[i] == 0) memcpy(p[i], q[i], 4*sizeof(int));
            else { // aniso refinements
              p[i][0] = h_only ? current : std::max(1, 2*(current+1)/3);
              p[i][1] = h_only ? current : std::max(1, 2*(current+1)/3);
            }
          }  
          int ii = idx[id][j];
          current = get_h_order(spaces[j]->get_element_order(id));
          if (ii >= 0) 
          {
            if (split[ii] != max_ref)
            {
              split[ii] = max_ref;
              if (split[ii] == 0) memcpy(p[ii], q[ii], 4*sizeof(int));
              else { // aniso refinements
                p[ii][0] = h_only ? current : std::max(1, 2*(current+1)/3);
                p[ii][1] = h_only ? current : std::max(1, 2*(current+1)/3);
              }
            }
          }
          if (ii < 0) // element not refined at all
          {
            split[k] = max_ref;
            if (split[k] == 0)
              for (int r = 0; r < 4; r++)
                p[k][r] = h_only ? current : std::max(1, (current+1)/2);
            else { // aniso refinements
              p[k][0] = h_only ? current : std::max(1, 2*(current+1)/3);
              p[k][1] = h_only ? current : std::max(1, 2*(current+1)/3);
            }
            esort[k][0] = id;
            esort[k][1] = j;
            k++;
          }
        }
      }    
    } 
  }

  for (i = 0; i < k; i++) // go over elements to be refined
  {
    int comp = esort[i][1];    
    int id = esort[i][0];
    Element* e;
    e = mesh[comp]->get_element(id);

    if (split[i] < 0)
      spaces[comp]->set_element_order(id, p[i][0]);
    else if (split[i] == 0) 
    {
      if (e->active)
        mesh[comp]->refine_element(id);
      for (j = 0; j < 4; j++)
        spaces[comp]->set_element_order(e->sons[j]->id, p[i][j]);
    }   
    else {
      if (e->active)
        mesh[comp]->refine_element(id, split[i]);
      for (j = 0; j < 2; j++)
        spaces[comp]->set_element_order(e->sons[ (split[i] == 1) ? j : j+2 ]->id, p[i][j]);
    }  

  }
  
  for (j = 0; j < num; j++)
    rsln[j]->enable_transform(true);


  verbose("Refined %d elements.", i);
  have_errors = false;
}


///// Unrefinements /////////////////////////////////////////////////////////////////////////////////

void H1OrthoHP::unrefine(Solution* sln1, Solution* sln2, double thr)
{
  
  if (!have_errors)
    error("Element errors have to be calculated first, see calc_error().");

  Mesh* mesh[2];
  mesh[0] = spaces[0]->get_mesh();
  mesh[1] = spaces[1]->get_mesh();


  int k = 0;
  if (mesh[0] == mesh[1]) // single mesh
  {
    Element* e;
    for_all_inactive_elements(e, mesh[0])
    {
      bool found = true;
      for (int i = 0; i < 4; i++)
        if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
          { found = false;  break; }
  
      if (found)
      {
        double sum1 = 0.0, sum2 = 0.0;
        int max1 = 0, max2 = 0;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL) 
          {
            sum1 += errors[0][e->sons[i]->id];
            sum2 += errors[1][e->sons[i]->id];
            int oo = spaces[0]->get_element_order(e->sons[i]->id);
            if (oo > max1) max1 = oo;
            oo = spaces[1]->get_element_order(e->sons[i]->id);
            if (oo > max2) max2 = oo;
          }
        if ((sum1 < thr * errors[esort[0][1]][esort[0][0]]) && 
            (sum2 < thr * errors[esort[0][1]][esort[0][0]]))
        {
          mesh[0]->unrefine_element(e->id);
          mesh[1]->unrefine_element(e->id);
          errors[0][e->id] = sum1;
          errors[1][e->id] = sum2;
          spaces[0]->set_element_order(e->id, max1);
          spaces[1]->set_element_order(e->id, max2);
          k++; // number of unrefined elements
        }
      }
    }
    for_all_active_elements(e, mesh[0])
    {
      for (int i = 0; i < 2; i++)
        if (errors[i][e->id] < thr/4 * errors[esort[0][1]][esort[0][0]])
        {
          int oo = get_h_order(spaces[i]->get_element_order(e->id));
          spaces[i]->set_element_order(e->id, std::max(oo - 1, 1));
          k++;
        }
    }
  }
  else // multimesh
  {
    for (int m = 0; m < 2; m++)
    {
      Element* e;
      for_all_inactive_elements(e, mesh[m])
      {
        bool found = true;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
            { found = false;  break; }
    
        if (found)
        {
          double sum = 0.0;
          int max = 0;
          for (int i = 0; i < 4; i++)
            if (e->sons[i] != NULL) 
            {
              sum += errors[m][e->sons[i]->id];
              int oo = spaces[m]->get_element_order(e->sons[i]->id);
              if (oo > max) max = oo;
            }
          if ((sum < thr * errors[esort[0][1]][esort[0][0]]))
          {
            mesh[m]->unrefine_element(e->id);
            errors[m][e->id] = sum;
            spaces[m]->set_element_order(e->id, max);
            k++; // number of unrefined elements
          }
        }
      }
      for_all_active_elements(e, mesh[m])
      {
        if (errors[m][e->id] < thr/4 * errors[esort[0][1]][esort[0][0]])
        {
          int oo = get_h_order(spaces[m]->get_element_order(e->id));
          spaces[m]->set_element_order(e->id, std::max(oo - 1, 1));
          k++;
        }
      }
    }
  }
  verbose("Unrefined %d elements.", k);
  have_errors = false;
}


//// error calculation /////////////////////////////////////////////////////////////////////////////

static double** cmp_err;
static int compare(const void* p1, const void* p2)
{ 
  const int2 (*e1) = ((const int2*) p1);
  const int2 (*e2) = ((const int2*) p2);
  return cmp_err[(*e1)[1]][(*e1)[0]] < cmp_err[(*e2)[1]][(*e2)[0]] ? 1 : -1; 
}


double H1OrthoHP::calc_error_n(int n, ...)
{
  int i, j, k;

  if (n != num) error("Wrong number of solutions.");
  
  va_list ap;
  va_start(ap, n);
  for (i = 0; i < n; i++) {
    sln[i] = va_arg(ap, Solution*);
    sln[i]->enable_transform(true);
  }
  for (i = 0; i < n; i++) {
    rsln[i] = va_arg(ap, Solution*);
    rsln[i]->enable_transform(true);
  }
  va_end(ap);
  
  nact = 0;
  for (j = 0; j < num; j++)
    nact += sln[j]->get_mesh()->get_num_active_elements();
  if (esort != NULL) delete [] esort;
  esort = new int2[nact];

  double total_error = 0.0,
         total_norm  = 0.0;

  double norms[n];
  memset(norms, 0, n*sizeof(double));
  for (i = 0; i < n; i++)
    norms[i] = sqr(h1_norm(rsln[i]));

  for (j = k = 0; j < num; j++)
  {
    Mesh* cmesh = sln[j]->get_mesh();
    Mesh* fmesh = rsln[j]->get_mesh();
    
    int max = cmesh->get_max_element_id();
    if (errors[j] != NULL) delete [] errors[j];
    errors[j] = new double[max];
    memset(errors[j], 0, sizeof(double) * max);
    
    Element* e;
    for_all_active_elements(e, cmesh)
    {
      sln[j]->set_active_element(e);
      update_limit_table(e->get_mode());
      for (i = 0; i < 4; i++)
      {
        sln[j]->push_transform(i);
  
        Element* fe = fmesh->get_element(e->id);
        if (fe->active || fe->sons[i] == NULL || !fe->sons[i]->active)
          error("Bad reference solution.");
        rsln[j]->set_active_element(fe->sons[i]);
  
        RefMap* crm = sln[j]->get_refmap();
        RefMap* frm = rsln[j]->get_refmap();
  
        double err  = int_h1_error(sln[j], rsln[j], crm, frm);
        total_norm += int_h1_norm (rsln[j], frm);
  
        errors[j][e->id] += err / norms[j];
        total_error += err;
  
        sln[j]->pop_transform();
      }
      esort[k][0] = e->id;
      esort[k++][1] = j;
    }
    
    for_all_inactive_elements(e, cmesh)
      errors[j][e->id] = -1.0;
  }
  
  assert(k == nact);
  cmp_err = errors;
  qsort(esort, nact, sizeof(int2), compare);
  
  have_errors = true;
  total_err = total_error /total_norm;
  
  return sqrt(total_error / total_norm);
}


//// energy errors /////////////////////////////////////////////////////////////////////////////////

class ErrorFn : public ScalarFunction
{
public:

  // ErrorFn subtracts two solutions. Used as input to the bilinear (energy) forms.
  ErrorFn(Solution* sln1, Solution* sln2)
  {
    this->sln1 = sln1;
    this->sln2 = sln2;

    element = sln1->get_active_element();
    order = std::max(sln1->get_fn_order(), sln2->get_fn_order());
    num_components = 1;
    assert(sln1->get_num_components() == 1);

    tables = NULL;
    sub_tables = &tables;
    update_nodes_ptr();
    
    set_quad_2d(sln1->get_quad_2d());
  }
  
  ~ErrorFn() { free(); }
  
  void free()
  {
    free_sub_tables(&tables);
  }

protected:

  Solution *sln1, *sln2;
  void* tables;

  virtual void precalculate(int order, int mask)
  {
    Quad2D* quad = quads[cur_quad];
    int np = quad->get_num_points(order);
    assert(!(mask & ~FN_DEFAULT));
    Node* node = new_node(FN_DEFAULT, np);

    sln1->set_quad_order(order, FN_DEFAULT);
    sln2->set_quad_order(order, FN_DEFAULT);
    
    scalar *val1 = sln1->get_fn_values();
    scalar *val2 = sln2->get_fn_values();
    
    scalar *d1dx, *d1dy, *d2dx, *d2dy;
    sln1->get_dx_dy_values(d1dx, d1dy);
    sln2->get_dx_dy_values(d2dx, d2dy);
    
    Trf *ctm1 = sln1->get_ctm(), *ctm2 = sln2->get_ctm();
    double mm = (ctm1->m[0] * ctm2->m[0] < 0) ? -2.0 : 2.0;
    
    for (int i = 0; i < np; i++)
    {
      node->values[0][0][i] = val1[i] - val2[i];
      node->values[0][1][i] = d1dx[i] - mm*d2dx[i];
      node->values[0][2][i] = d1dy[i] - mm*d2dy[i];
    }
    
    replace_cur_node(node);
  }

};


double H1OrthoHP::calc_energy_error_n(int n, ...)
{
  int i, j, k;

  if (n != num) error("Wrong number of solutions.");
  
  // obtain solutions and bilinear forms
  va_list ap;
  va_start(ap, n);
  for (i = 0; i < n; i++) {
    sln[i] = va_arg(ap, Solution*);
    sln[i]->enable_transform(false);
    sln[i]->set_quad_2d(&g_quad_2d_std);
  }
  for (i = 0; i < n; i++) {
    rsln[i] = va_arg(ap, Solution*);
    rsln[i]->enable_transform(false);
    rsln[i]->set_quad_2d(&g_quad_2d_std);
  }
  biform_t bi[10][10];
  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++)
      bi[i][j] = va_arg(ap, biform_t);
  va_end(ap);

  // prepare multi-mesh traversal and error arrays
  Mesh* meshes[2*num];
  Transformable* tr[2*num];
  Traverse trav;
  nact = 0;
  for (i = 0; i < num; i++)
  {
    meshes[i] = sln[i]->get_mesh();
    meshes[i+num] = rsln[i]->get_mesh();
    tr[i] = sln[i];
    tr[i+num] = rsln[i];

    nact += sln[i]->get_mesh()->get_num_active_elements();

    int max = meshes[i]->get_max_element_id();
    if (errors[i] != NULL) delete [] errors[i];
    errors[i] = new double[max];
    memset(errors[i], 0, sizeof(double) * max);
  }
  
  double total_norm = 0.0;
  double norms[num];
  memset(norms, 0, num*sizeof(double));
  double total_error = 0.0;
  if (esort != NULL) delete [] esort;
  esort = new int2[nact];

  // 
  ErrorFn* err[num];
  Element** ee;
  trav.begin(2*num, meshes, tr);
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    for (i = 0; i < num; i++)
      err[i] = new ErrorFn(sln[i], rsln[i]);

    for (i = 0; i < num; i++)
    {
      RefMap* rm1 = rsln[i]->get_refmap();
      for (j = 0; j < num; j++)
      {
        RefMap* rm2 = rsln[j]->get_refmap();
        double e, t;
        if (bi[i][j] != NULL)
        {
          #ifndef COMPLEX
          e = fabs(bi[i][j]( err[i],  err[j], rm1, rm2)) * 0.25;
          t = fabs(bi[i][j](rsln[i], rsln[j], rm1, rm2));
          #else
          e = std::abs(bi[i][j]( err[i],  err[j], rm1, rm2)) * 0.25;
          t = std::abs(bi[i][j](rsln[i], rsln[j], rm1, rm2));
          #endif
        }
  
        norms[i] += t;
        total_norm  += t;
        total_error += e;
        errors[i][ee[i]->id] += e;
      }
    }
    for (i = 0; i < num; i++)
      delete err[i];
  }
  trav.finish();

  Element* e;
  k = 0;
  for (i = 0; i < num; i++)
    for_all_active_elements(e, meshes[i]) {
      esort[k][0] = e->id;
      esort[k++][1] = i;
      errors[i][e->id] /= norms[i];
    }

  assert(k == nact);
  cmp_err = errors;
  qsort(esort, nact, sizeof(int2), compare);
    
  for (i = 0; i < num; i++)
  {
    sln[i]->enable_transform(true);
    rsln[i]->enable_transform(true);
  }
  
  have_errors = true;
  total_err = total_error / total_norm;
  return sqrt(total_error / total_norm);
}