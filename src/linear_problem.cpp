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

#include "common.h"
#include "limit_order.h"
#include "discrete_problem.h"
#include "linear_problem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"
#include "views/view.h"
#include "views/vector_view.h"
#include "tuple.h"
#include "norm.h"


LinearProblem::LinearProblem() : DiscreteProblem() {};
LinearProblem::LinearProblem(WeakForm* wf_) : DiscreteProblem(wf_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Space* s_) : DiscreteProblem(wf_, s_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Tuple<Space*> spaces_) : DiscreteProblem(wf_, spaces_) {};
LinearProblem::~LinearProblem() {};

void LinearProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly, bool is_complex)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof == 0 in LinearProblem::assemble().");
  Vector* dir_ext = new AVector(ndof, is_complex);
  // The vector dir represents the contribution of the Dirichlet lift, 
  // and for linear problems it has to be subtracted from the right hand side.
  // The NULL stands for the initial coefficient vector that is not used.
  DiscreteProblem::assemble(NULL, mat_ext, dir_ext, rhs_ext, rhsonly);
  // FIXME: Do we really need to handle the real and complex cases separately?
  if (is_complex) for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get_cplx(i));
  else for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get(i));
  delete dir_ext;
}

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> solutions, 
                  MatrixSolverType matrix_solver, bool is_complex) 
{
  // Initialize the linear problem.
  LinearProblem lp(wf, spaces);
  //info("ndof = %d", lp.get_num_dofs());

  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, get_num_dofs(spaces), mat, rhs, solver, is_complex);

  // Assemble stiffness matrix and rhs.
  bool rhsonly = false;
  lp.assemble(mat, rhs, rhsonly, is_complex);

  //mat->print();

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  for (int i=0; i<solutions.size(); i++) {
    solutions[i]->set_fe_solution(spaces[i], rhs);
  }
}

// Solve a typical linear problem using automatic adaptivity.
// Feel free to adjust this function for more advanced applications.
bool solve_linear_adapt(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> slns, 
                        MatrixSolverType matrix_solver, Tuple<Solution *> ref_slns, Tuple<int> proj_norms, 
                        RefinementSelectors::Selector* selector, AdaptivityParamType* apt,
                        Tuple<WinGeom *> sln_win_geom, Tuple<WinGeom *> mesh_win_geom, 
                        bool verbose, Tuple<ExactSolution *> exact_slns, bool is_complex) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity parameters.
  double err_stop = apt->err_stop; 
  int ndof_stop = apt->ndof_stop;
  double threshold = apt->threshold;
  int strategy = apt->strategy; 
  int mesh_regularity = apt->mesh_regularity;
  double to_be_processed = apt->to_be_processed;

  // Number of physical fields in the problem.
  int num_comps = spaces.size();

  // Number of degreeso of freedom 
  int ndof = get_num_dofs(spaces);

  // Number of exact solutions given.
  if (exact_slns.size() != 0 && exact_slns.size() != num_comps) 
    error("Number of exact solutions does not match number of equations.");
  bool is_exact_solution;
  if (exact_slns.size() == num_comps) is_exact_solution = true;
  else is_exact_solution = false;

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, is_complex);

  // Initialize views.
  ScalarView* sview[H2D_MAX_COMPONENTS];
  OrderView*  oview[H2D_MAX_COMPONENTS];
  for (int i = 0; i < num_comps; i++) {
    char* title = (char*)malloc(100*sizeof(char));
    if (sln_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Solution", i); 
      else sprintf(title, "Solution %d", i); 
      sview[i] = new ScalarView(title, sln_win_geom[i]);
      sview[i]->show_mesh(false);
    }
    else sview[i] = NULL;
    if (mesh_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Mesh", i); 
      else sprintf(title, "Mesh %d", i); 
      oview[i] = new OrderView(title, mesh_win_geom[i]);
    }
    else oview[i] = NULL;
  }

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Conversion from Tuple<Solution *> to Tuple<MeshFunction *>
  // so that project_global() below compiles. 
  Tuple<MeshFunction *> ref_slns_mf;
  for (int i = 0; i < num_comps; i++) {
    MeshFunction *s = (MeshFunction*)ref_slns[i];
    ref_slns_mf.push_back(s);
  }

  int as = 1; bool done = false;
  do
  {
    if (verbose) {
      info("---- Adaptivity step %d:", as);
      info("Solving on reference mesh.");
    }

    // Construct globally refined reference mesh(es)
    // and setup reference space(s).
    Tuple<Space *> ref_spaces;
    for (int i = 0; i < num_comps; i++) {
      Mesh *ref_mesh = new Mesh();
      ref_mesh->copy(spaces[i]->get_mesh());
      ref_mesh->refine_all_elements();
      ref_spaces.push_back(spaces[i]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[i]->copy_orders(spaces[i], order_increase);
    }

    // Solve the reference problem.
    solve_linear(ref_spaces, wf, ref_slns, matrix_solver);

    // Project the reference solution on the coarse mesh.
    if (verbose) info("Projecting reference solution on coarse mesh.");
    project_global(spaces, ref_slns_mf, slns, proj_norms, NULL, is_complex); // NULL means that we do not want to know the resulting coefficient vector.

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution (first component only).
    for (int i = 0; i < num_comps; i++) {
      if (sview[i] != NULL) sview[i]->show(slns[i]);
      if (oview[i] != NULL) oview[i]->show(spaces[i]);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    if (verbose) info("Calculating error (est).");
    H1Adapt hp(spaces);
    hp.set_solutions(slns, ref_slns);
    hp.calc_error();
 
    // Calculate error estimate for each solution component.   
    double err_est[H2D_MAX_COMPONENTS];
    double err_est_total = 0;
    for (int i = 0; i < num_comps; i++) {
      err_est[i] = h1_error(slns[i], ref_slns[i]) * 100;
      err_est_total += err_est[i]*err_est[i];
    }
    err_est_total = sqrt(err_est_total);

    // Calculate exact error for each solution component.   
    double err_exact[H2D_MAX_COMPONENTS];
    double err_exact_total = 0;
    if (is_exact_solution == true) {
      if (verbose) info("Calculating error (exact).");
      for (int i = 0; i < num_comps; i++) {
        err_exact[i] = h1_error(slns[i], exact_slns[i]) * 100;
        err_exact_total += err_exact[i]*err_exact[i];
      }
      err_exact_total = sqrt(err_exact_total);
    }

    // Report results.
    if (verbose) {
      if (num_comps == 1) {
        info("ndof: %d, ref_ndof: %d, err_est: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_total);
        if (is_exact_solution == true) info("err_exact: %g%%", err_exact[0]);
      }
      else {
        for (int i = 0; i < num_comps; i++) {
          info("ndof[%d]: %d, ref_ndof[%d]: %d, err_est[%d]: %g%%", 
               i, spaces[i]->get_num_dofs(), i, ref_spaces[i]->get_num_dofs(),
               i, err_est[i]);
          if (is_exact_solution == true) info("err_exact[%d]: %g%%", i, err_exact[i]);
        }
        info("ndof: %d, ref_ndof: %d, err_est: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_total);
      }
    }

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(spaces), err_est_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    if (is_exact_solution == true) {
      graph_dof_exact.add_values(get_num_dofs(spaces), err_exact_total);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_total);
      graph_cpu_exact.save("conv_cpu_exact.dat");
    }

    // If err_est too large, adapt the mesh.
    if (err_est_total < err_stop) done = true;
    else {
      if (verbose) info("Adapting the coarse mesh.");
      done = hp.adapt(selector, threshold, strategy, mesh_regularity, to_be_processed);

      if (get_num_dofs(spaces) >= ndof_stop) done = true;
    }

    as++;
  }
  while (done == false);

  // Close visualization windows.
  /* FIXME: this should be uncommented but it causes segfaults
  for (int i = 0; i < num_comps; i++) {
    if (sview[i] != NULL) sview[i]->close(); 
    if (oview[i] != NULL) oview[i]->close(); 
  }
  */

  if (verbose) info("Total running time: %g s", cpu_time.accumulated());
}


