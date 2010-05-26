//////  Bilinear and linear forms - axisymmetric arrangement  ////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar int_x_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (e->x[i] * u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_x_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

//////////   Eq 1   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar jacobian_0_0(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][0]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, vj, vi, e) +
         (Sr[e->marker - 1][0]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_surf_0_0(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar residual_0(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][0] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * vj->val[i];
  return result;
}

//////////   Eq 2   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar jacobian_1_1(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][1]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, vj, vi, e) +
         (Sr[e->marker - 1][1]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_surf_1_1(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_1_0(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][1][0]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar residual_1(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][1] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * vj->val[i];
  return result;
}

//////////   Eq 3   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar jacobian_2_2(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][2]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, vj, vi, e) +
         (Sr[e->marker - 1][2]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_surf_2_2(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_2_1(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][2][1]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar residual_2(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][2] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * vj->val[i];
  return result;
}

//////////   Eq 4   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar jacobian_3_3(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->marker - 1][3]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, vj, vi, e) +
         (Sr[e->marker - 1][3]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_surf_3_3(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar jacobian_3_2(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->marker - 1][3][2]) * int_x_u_v<Real, Scalar>(n, wt, vj, vi, e);
}

template<typename Real, typename Scalar>
Scalar residual_3(int n, double *wt, Func<Scalar> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->marker - 1][3] / k_eff) * (nu[e->marker - 1][0] * Sf[e->marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->marker - 1][1] * Sf[e->marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->marker - 1][2] * Sf[e->marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->marker - 1][3] * Sf[e->marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * vj->val[i];
  return result;
}
