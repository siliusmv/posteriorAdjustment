#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "cgeneric.h"
#include "smat-operations.h"
#include "spde-precision.h"
#include "test-funcs.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

int main(int argc, char * argv[]) {

  assert(argc == 3);
  double log_rho = atof(argv[1]);
  double log_sigma = atof(argv[2]);

  inla_cgeneric_data_tp data = read_cgeneric_data_from_dir("../cgeneric-data");

  inla_cgeneric_smat_tp Q = spde_precision(log_rho, log_sigma, data.mats, data.smats);
  Q.name = Calloc(1, char);
  strcpy(Q.name, "Q");
  print_smat_as_mat(Q);

  free_smat(&Q);

  return 0;
}
