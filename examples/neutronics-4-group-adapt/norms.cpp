double error_fn_l2_axisym(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o, H2D_FN_VAL);
  sln2->set_quad_order(o, H2D_FN_VAL);

  scalar *uval, *vval;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();

	double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i]*sqr(uval[i] - vval[i]));
  return 2*M_PI*result;
}


// function used to calculate L2 norm of the solution
double norm_fn_l2_axisym(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 *sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o, H2D_FN_VAL);

  scalar* uval = sln->get_fn_values();

	double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i]*sqr(uval[i]));
  return 2*M_PI*result;
}


double l2_error_rel_sq(MeshFunction* sln1, MeshFunction* sln2)
{
  double error = calc_error(error_fn_l2_axisym, sln1, sln2);
  double norm = calc_norm(norm_fn_l2_axisym, sln2);
  return error/norm;
}

double l2_norm_rel_sq(MeshFunction* sln)
{
  return calc_norm(norm_fn_l2_axisym, sln);
}

double l2_error_total_rel_sq(Tuple<Solution*>& slns1, Tuple<Solution*>& slns2)
{
	Tuple<Solution*>::iterator it1, it2;
	double error = 0.0, norm = 0.0;
	
	for (it1=slns1.begin(), it2=slns2.begin(); it1 < slns1.end(); it1++, it2++) {
		assert(it2 < slns2.end());
		error += calc_error(error_fn_l2_axisym, *it1, *it2);
		norm += calc_norm(norm_fn_l2_axisym, *it2);
	}
	
	return error/norm;
}
		  
