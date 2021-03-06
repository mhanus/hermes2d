// This file is part of Hermes2D
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

#ifndef __H2D_SOLVER_NOX_H_
#define __H2D_SOLVER_NOX_H_

#include "itersolver.h"
#include "feproblem.h"
#include "precond.h"
#include "solver_epetra.h"

#ifdef HAVE_NOX
#include <NOX.H>
#include <NOX_Epetra.H>
#endif

class NoxProblemInterface;

/// Encapsulation of NOX nonlinear solver
///
/// @ingroup solvers
class NoxSolver : public IterSolver
{
public:
	NoxSolver(FeProblem *problem);
	virtual ~NoxSolver();

	bool set_init_sln(double *ic);
	bool set_init_sln(EpetraVector *ic);
	virtual bool solve();

	int get_num_iters() { return num_iters; }
  double get_residual()  { return residual;}
  int get_num_lin_iters() { return num_lin_iters; }
  double get_achieved_tol()  { return achieved_tol;}

	// settings for the solver
	void set_nl_method(const char *par);
  void set_output_flags(int flags) { output_flags = flags; }

	// linear solver setters
	void set_ls_type(const char *type) { ls_type = type; }
	void set_ls_max_iters(int iters) { ls_max_iters = iters; }
	void set_ls_tolerance(double tolerance) { ls_tolerance = tolerance; }
	void set_ls_sizeof_krylov_subspace(int size) { ls_sizeof_krylov_subspace = size; }

	// convergence params
#ifdef HAVE_NOX
  void set_norm_type(NOX::Abstract::Vector::NormType type)  { conv.norm_type = type; }
  void set_scale_type(NOX::StatusTest::NormF::ScaleType type)  { conv.stype = type; }
#endif
  void set_conv_iters(int iters)        { conv.max_iters = iters; }
	void set_conv_abs_resid(double resid) { conv_flag.absresid = 1; conv.abs_resid = resid; }
	void set_conv_rel_resid(double resid) { conv_flag.relresid = 1; conv.rel_resid = resid; }
	void set_conv_update(double update)   { conv_flag.update = 1; conv.update = update; }
	void set_conv_wrms(double rtol, double atol) {
		conv_flag.wrms = 1;
		conv.wrms_rtol = rtol;
		conv.wrms_atol = atol;
	}

	void set_precond(Precond *pc);
  void set_precond(const char *pc);

protected:
#ifdef HAVE_NOX
	Teuchos::RCP<NoxProblemInterface> interface;
#endif
	int num_iters;
  double residual;
  int num_lin_iters;
  double achieved_tol;

	const char *nl_dir;

  int output_flags;
	const char *ls_type;
	int ls_max_iters;
	double ls_tolerance;
	int ls_sizeof_krylov_subspace;
  bool precond_yes;
  const char* precond_type;
	// convergence params
	struct conv_t {
		int max_iters;
		double abs_resid;
		double rel_resid;
#ifdef HAVE_NOX
    NOX::Abstract::Vector::NormType norm_type;
    NOX::StatusTest::NormF::ScaleType stype;
#endif
		double update;
		double wrms_rtol;
		double wrms_atol;
	} conv;

	struct conv_flag_t {
		unsigned absresid:1;
		unsigned relresid:1;
		unsigned wrms:1;
		unsigned update:1;
	} conv_flag;
};

#endif
