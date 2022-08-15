

#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <omp.h>

// #include <gsl/gsl_spmatrix.h>

#include "cgeneric.h"
#include "read-write.h"
#include "spde-funcs.h"
#include "unused-code.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))


int main() {

  inla_cgeneric_mat_tp A;
  A.name = NULL;
  A.nrow = 2;
  A.ncol = 3;
  A.x = Calloc(6, double);
  for (int i = 0; i < 6; i++) {
    A.x[i] = i + 1;
  }
  print_mat(A);

  inla_cgeneric_smat_tp D;
  D.name = NULL;
  D.nrow = 2;
  D.ncol = 3;
  D.n = 4;
  int i[] = {0, 1, 0, 1};
  int j[] = {0, 0, 1, 2};
  double x[] = {1, 2, 3, 6};
  D.i = i;
  D.j = j;
  D.x = x;
  print_smat(D);
  print_mat(smat_to_mat(D));

  // gsl_spmatrix *D2 = inla_smat_to_gsl(&D);
  // printf("D2: Type: %s\n", gsl_spmatrix_type(D2));
  // gsl_spmatrix_fprintf(stdout, D2, "%.1f");

  // gsl_spmatrix_free(D2);

  // D2 = inla_smat_to_gsl(&D);
  // printf("D2 second time: Type: %s\n", gsl_spmatrix_type(D2));
  // gsl_spmatrix_fprintf(stdout, D2, "%.1f");

  // gsl_matrix *D3 = gsl_matrix_alloc(2, 3);
  // gsl_spmatrix_sp2d(D3, D2);
  // gsl_matrix_fprintf(stdout, D3, "%.1f");

  // gsl_spmatrix *D4 = gsl_spmatrix_alloc_nzmax(2, 3, 5, GSL_SPMATRIX_CSC);
  // gsl_spmatrix_add(D4, D2, D2);
  // printf("D4: Type: %s\n", gsl_spmatrix_type(D4));
  // gsl_spmatrix_fprintf(stdout, D4, "%.1f");

  // gsl_spmatrix *D5 = csc_to_coo(D4);
  // printf("D5: Type: %s\n", gsl_spmatrix_type(D5));
  // gsl_spmatrix_fprintf(stdout, D5, "%.1f");
  // D5->sptype = GSL_SPMATRIX_CSC;
  // printf("D5: Type: %s\n", gsl_spmatrix_type(D5));
  // D5->sptype = GSL_SPMATIX_COO;
  // printf("D5: Type: %s\n", gsl_spmatrix_type(D5));

  // gsl_spmatrix *D6 = gsl_spmatrix_ccs(D5);
  // printf("D6: Type: %s\n", gsl_spmatrix_type(D6));
  // gsl_spmatrix_fprintf(stdout, D6, "%.1f");
  // gsl_spmatrix_add(D6, D4, D6);
  // printf("D6: Type: %s\n", gsl_spmatrix_type(D6));
  // gsl_spmatrix_fprintf(stdout, D6, "%.1f");

  // gsl_spmatrix_free(D2);
  // gsl_matrix_free(D3);
  // gsl_spmatrix_free(D4);
  // gsl_spmatrix_free(D5);
  // gsl_spmatrix_free(D6);

  double start_time = omp_get_wtime();

  inla_cgeneric_data_tp data = read_cgeneric_data_from_dir("../cgeneric-data");
  printf("n_mat: %d\n", data.n_mat);
  printf("n_smat: %d\n", data.n_smat);

  printf("2\n");

  // Order the matrices correctly
  inla_cgeneric_mat_tp *tmp_mat_pt = data.mats[0];
  data.mats[0] = data.mats[2];
  data.mats[2] = tmp_mat_pt;
  inla_cgeneric_smat_tp *tmp_smat_pt = data.smats[0];
  data.smats[0] = data.smats[2];
  data.smats[2] = tmp_smat_pt;

  inla_cgeneric_smat_tp Q = spde_precision(log(5), log(2), data.mats, data.smats);
  Q.name = NULL;
  // print_smat(Q);
  print_smat_as_mat(Q);
  printf("%d %d %d\n", Q.n, Q.nrow, Q.ncol);

  upper_diag(&Q);
  sort_smat(&Q);
  print_smat(Q);
  printf("upper diag only: %d %d %d\n", Q.n, Q.nrow, Q.ncol);

  //inla_cgeneric_smat_tp Q2 = copy_smat(&Q);

  // block_diag_smat(&Q, 2);
  // print_smat(Q);
  // printf("block diag: %d %d %d\n", Q.n, Q.nrow, Q.ncol);

  print_mat(smat_to_mat(Q));

  printf("Time passed: %f\n", omp_get_wtime() - start_time);

  sleep(2);

  printf("Time passed: %f\n", omp_get_wtime() - start_time);


  return 0;
}
