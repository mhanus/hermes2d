template<typename Real, typename Scalar>
Scalar projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_x_u_v<Real, Scalar>(n, wt, u, v, e) + int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_x_u_v<Real, Scalar>(n, wt, ext->fn[0], v, e) + int_x_grad_u_grad_v<Real, Scalar>(n, wt, ext->fn[0], v, e);
}
/*
//////////   Eq 1   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar projection_liform_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return biform_0_0<Real, Scalar>(n, wt, ext->fn[0], v, e, ext); 
}
template<typename Real, typename Scalar>
Scalar projection_liform_surf_0(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return biform_surf_0_0<Real, Scalar>(n, wt, ext->fn[0], v, e, ext);
}

//////////   Eq 2   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar projection_liform_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	return biform_1_0<Real, Scalar>(n, wt, ext->fn[0], v, e, ext) + biform_1_1<Real, Scalar>(n, wt, ext->fn[1], v, e, ext);
}
template<typename Real, typename Scalar>
Scalar projection_liform_surf_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return biform_surf_1_1<Real, Scalar>(n, wt, ext->fn[0], v, e, ext);
}

//////////   Eq 3   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar projection_liform_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	return biform_2_1<Real, Scalar>(n, wt, ext->fn[0], v, e, ext) + biform_2_2<Real, Scalar>(n, wt, ext->fn[1], v, e, ext);
}
template<typename Real, typename Scalar>
Scalar projection_liform_surf_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return biform_surf_2_2<Real, Scalar>(n, wt, ext->fn[0], v, e, ext);
}

//////////   Eq 4   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar projection_liform_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	return biform_3_2<Real, Scalar>(n, wt, ext->fn[0], v, e, ext) + biform_3_3<Real, Scalar>(n, wt, ext->fn[1], v, e, ext);
}
template<typename Real, typename Scalar>
Scalar projection_liform_surf_3(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return biform_surf_3_3<Real, Scalar>(n, wt, ext->fn[0], v, e, ext);
}
*/

class CustomProjector : public Projection
{
	public:
		CustomProjector(Solution *sol1, Solution *sol2, Solution *sol3, Solution *sol4, 
										Space	*space1, Space *space2, Space *space3, Space *space4,  
										PrecalcShapeset* pss1, PrecalcShapeset* pss2, PrecalcShapeset* pss3, PrecalcShapeset* pss4 ) : 
					Projection(	4, sol1, sol2, sol3, sol4, 
											space1, space2, space3, space4,  
											pss1, pss2, pss3, pss4) { };
										
		void project(Solution *res1, Solution *res2, Solution *res3, Solution *res4);
		scalar* project();
};

void CustomProjector::project(Solution *res1, Solution *res2, Solution *res3, Solution *res4)
{
  WeakForm wf(num);
  
/*  wf.add_biform(0, 0, callback(biform_0_0));
  wf.add_biform(1, 1, callback(biform_1_1));
  wf.add_biform(1, 0, callback(biform_1_0));
  wf.add_biform(2, 2, callback(biform_2_2));
  wf.add_biform(2, 1, callback(biform_2_1));
  wf.add_biform(3, 3, callback(biform_3_3));
  wf.add_biform(3, 2, callback(biform_3_2));
  wf.add_liform(0, callback(projection_liform_0), marker_core, 1, slns[0]);
  wf.add_liform(1, callback(projection_liform_1), marker_core, 2, slns[0], slns[1]);
  wf.add_liform(2, callback(projection_liform_2), marker_core, 2, slns[1], slns[2]);
  wf.add_liform(3, callback(projection_liform_3), marker_core, 2, slns[2], slns[3]);
  wf.add_biform_surf(0, 0, callback(biform_surf_0_0), bc_vacuum);
  wf.add_biform_surf(1, 1, callback(biform_surf_1_1), bc_vacuum);
  wf.add_biform_surf(2, 2, callback(biform_surf_2_2), bc_vacuum);
  wf.add_biform_surf(3, 3, callback(biform_surf_3_3), bc_vacuum);
  wf.add_liform_surf(0, callback(projection_liform_surf_0), bc_vacuum, 1, slns[0]);
  wf.add_liform_surf(1, callback(projection_liform_surf_1), bc_vacuum, 1, slns[1]);
  wf.add_liform_surf(2, callback(projection_liform_surf_2), bc_vacuum, 1, slns[2]);
  wf.add_liform_surf(3, callback(projection_liform_surf_3), bc_vacuum, 1, slns[3]);
*/  
  for (int i = 0; i < num; i++)
  {
    wf.add_biform(i, i, callback(projection_biform));
    wf.add_liform(i, callback(projection_liform), H2D_ANY, 1, slns[i]);
  }
	
  LinSystem ps(&wf, solver);
  ps.set_spaces(num, spaces[0], spaces[1], spaces[2], spaces[3]);
  ps.set_pss(num, pss[0], pss[1], pss[2], pss[3]);
  ps.assemble();
  
  ps.solve(4, res1, res2, res3, res4);
}

scalar* CustomProjector::project()
{
  WeakForm wf(num);
  
/*  wf.add_biform(0, 0, callback(biform_0_0));
  wf.add_biform(1, 1, callback(biform_1_1));
  wf.add_biform(1, 0, callback(biform_1_0));
  wf.add_biform(2, 2, callback(biform_2_2));
  wf.add_biform(2, 1, callback(biform_2_1));
  wf.add_biform(3, 3, callback(biform_3_3));
  wf.add_biform(3, 2, callback(biform_3_2));
  wf.add_liform(0, callback(projection_liform_0), marker_core, 1, slns[0]);
  wf.add_liform(1, callback(projection_liform_1), marker_core, 2, slns[0], slns[1]);
  wf.add_liform(2, callback(projection_liform_2), marker_core, 2, slns[1], slns[2]);
  wf.add_liform(3, callback(projection_liform_3), marker_core, 2, slns[2], slns[3]);
  wf.add_biform_surf(0, 0, callback(biform_surf_0_0), bc_vacuum);
  wf.add_biform_surf(1, 1, callback(biform_surf_1_1), bc_vacuum);
  wf.add_biform_surf(2, 2, callback(biform_surf_2_2), bc_vacuum);
  wf.add_biform_surf(3, 3, callback(biform_surf_3_3), bc_vacuum);
  wf.add_liform_surf(0, callback(projection_liform_surf_0), bc_vacuum, 1, slns[0]);
  wf.add_liform_surf(1, callback(projection_liform_surf_1), bc_vacuum, 1, slns[1]);
  wf.add_liform_surf(2, callback(projection_liform_surf_2), bc_vacuum, 1, slns[2]);
  wf.add_liform_surf(3, callback(projection_liform_surf_3), bc_vacuum, 1, slns[3]);
*/  
  for (int i = 0; i < num; i++)
  {
    wf.add_biform(i, i, callback(projection_biform));
    wf.add_liform(i, callback(projection_liform), H2D_ANY, 1, slns[i]);
  }

  LinSystem ps(&wf, solver);
  ps.set_spaces(num, spaces[0], spaces[1], spaces[2], spaces[3]);
  ps.set_pss(num, pss[0], pss[1], pss[2], pss[3]);
  ps.assemble();

  ps.solve(0);
  const scalar* sln_vec = ps.get_solution_vec();
  int ndofs = ps.get_num_dofs();
  vec = new scalar[ndofs];
  memcpy(vec, sln_vec, ndofs * sizeof(scalar));

  return vec;
}

