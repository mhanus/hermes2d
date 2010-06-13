#include "hermes2d.h"

// This test makes sure that example 02-space works correctly.

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);            // original L-shape domain
  //mloader.load("domain_quad.mesh", &mesh);     // reference square
  //mloader.load("domain_tri.mesh", &mesh);      // reference triangle

  // sample element refinement, to see more basis functions
  //mesh.refine_all_elements();

  // create a shapeset and an H1 space
  H1Shapeset shapeset;
  H1Space space(&mesh, &shapeset);

  // new code for the test
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof[10], dof_max[10];
  int success = 1;
  // testing all poly degrees between 1 and 10
  for (int i=1; i <= 10; i++) {
    space.set_uniform_order(i);
    assign_dofs(&space);
    n_dof[i-1] = space.get_num_dofs();
    dof_max[i-1] = space.get_max_dof();

    printf("n_dof = %d\n", n_dof[i-1]);
    printf("dof_max = %d\n", dof_max[i-1]);
  }
  if (n_dof[0] != 8 || dof_max[0] != 7) success = 0;
  if (n_dof[1] != 21 || dof_max[1] != 20) success = 0;
  if (n_dof[2] != 40 || dof_max[2] != 39) success = 0;
  if (n_dof[3] != 65 || dof_max[3] != 64) success = 0;
  if (n_dof[4] != 96 || dof_max[4] != 95) success = 0;
  if (n_dof[5] != 133 || dof_max[5] != 132) success = 0;
  if (n_dof[6] != 176 || dof_max[6] != 175) success = 0;
  if (n_dof[7] != 225 || dof_max[7] != 224) success = 0;
  if (n_dof[8] != 280 || dof_max[8] != 279) success = 0;
  if (n_dof[9] != 341 || dof_max[9] != 340) success = 0;

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

