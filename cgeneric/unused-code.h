#pragma once

#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#include "cgeneric.h"
#include "spde-funcs.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))


double *matrix_multiplication(double *A, double *B, int nrow_A, int ncol_A, int nrow_B, int ncol_B) {
  assert(ncol_A == nrow_B);
  double *res = Calloc(nrow_A * ncol_B, double);
  double tmp;
  for (int i = 0; i < ncol_B; i++) {
    for (int j = 0; j < nrow_A; j++) {
      tmp = 0;
      for (int k = 0; k < ncol_A; k++) {
	tmp += A[j + k * nrow_A] * B[k + i * nrow_B];
      }
      res[j + i * ncol_B] = tmp;
    }
  }
  return res;
}

// Matrices have been changed from column major to row major. Not sure this works anymore
// inla_cgeneric_vec_tp mat_vec_mult(inla_cgeneric_mat_tp A, inla_cgeneric_vec_tp x) {
//   assert(A.ncol == x.len);
//   inla_cgeneric_vec_tp res;
//   res.name = NULL;
//   res.len = A.nrow;
//   res.doubles = matrix_multiplication(A.x, x.doubles, A.nrow, A.ncol, x.len, 1);
//   return res;
// }
// 
// inla_cgeneric_mat_tp mat_mat_mult(inla_cgeneric_mat_tp A, inla_cgeneric_mat_tp B) {
//   assert(A.ncol == B.nrow);
//   inla_cgeneric_mat_tp res;
//   res.name = NULL;
//   res.nrow = A.nrow;
//   res.ncol = B.ncol;
//   res.x = matrix_multiplication(A.x, B.x, A.nrow, A.ncol, B.nrow, B.ncol);
//   return res;
// }
// 
// inla_cgeneric_mat_tp smat_to_mat(inla_cgeneric_smat_tp A) {
//   inla_cgeneric_mat_tp res;
//   res.nrow = A.nrow;
//   res.ncol = A.ncol;
//   res.name = A.name;
//   res.x = Calloc(res.nrow * res.ncol, double);
//   for (int i = 0; i < A.n; i++) {
//     res.x[A.i[i] + A.j[i] * res.nrow] = A.x[i];
//   }
//   return res;
// }
// 
// inla_cgeneric_smat_tp mat_to_smat(inla_cgeneric_mat_tp A) {
//   inla_cgeneric_smat_tp res;
//   res.nrow = A.nrow;
//   res.ncol = A.ncol;
//   res.name = A.name;
//   int n = 0;
//   for (int i = 0; i < res.nrow * res.ncol; i++) {
//     if (A.x[i] != 0) {
//       n += 1;
//     }
//   }
//   res.n = n;
//   res.x = Calloc(n, double);
//   res.i = Calloc(n, int);
//   res.j = Calloc(n, int);
//   int counter = 0;
//   for (int k = 0; k < res.nrow * res.ncol; k++) {
//     if (A.x[k] != 0) {
//       res.x[counter] = A.x[k];
//       res.i[counter] = k%res.nrow;
//       res.j[counter] = k /res.nrow;
//       counter += 1;
//     }
//   }
//   return res;
// }

void print_int_locs(int_loc_tp x) {
  for (int i = 0; i < x.num_ints; i++) {
    printf("%d: ", i);
    for (int j = 0; j < x.n[i]; j++) {
      printf("%d ", x.locs[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_smat(inla_cgeneric_smat_tp A) {
  printf("%s: (n: %d, nrow: %d, ncol: %d):\n", A.name, A.n, A.nrow, A.ncol);
  printf("i j x\n");
  for (int i = 0; i < A.n; i++) {
    printf("%d %d %f\n", A.i[i], A.j[i], A.x[i]);
  }
  printf("\n");
}

// Matrices have been changed from column major to row major. Not sure this works anymore
// void print_mat(inla_cgeneric_mat_tp A) {
//   printf("%s: (nrow: %d, ncol: %d):\n", A.name, A.nrow, A.ncol);
//   for (int i = 0; i < A.nrow; i++) {
//     for (int j = 0; j < A.ncol; j++) {
//       printf("%f ", A.x[i + j * A.nrow]);
//     }
//     printf("\n");
//   }
//   printf("\n");
// }
// 
// void print_smat_as_mat(inla_cgeneric_smat_tp A) {
//   inla_cgeneric_mat_tp tmp = smat_to_mat(A);
//   printf("%s: (nrow: %d, ncol: %d):\n", tmp.name, tmp.nrow, tmp.ncol);
//   for (int i = 0; i < tmp.nrow; i++) {
//     for (int j = 0; j < tmp.ncol; j++) {
//       if (tmp.x[i + j * tmp.nrow] == 0) {
// 	printf("* ");
//       } else {
// 	printf("%f ", tmp.x[i + j * tmp.nrow]);
//       }
//     }
//     printf("\n");
//   }
//   printf("\n");
// }

// We assume that there are no duplicate elements inside each of the index-structs
double sum_of_multiplied_intersects(nonzero_index_tp x, nonzero_index_tp y) {
  double res = 0;
  for (int i = 0; i < x.n; i++) {
    for (int j = 0; j < y.n; j++) {
      if (x.index[i] == y.index[j]) {
	res += x.val[i] * y.val[j];
	break;
      }
    }
  }
  return res;
}

inla_cgeneric_smat_tp smat_smat_mult(const inla_cgeneric_smat_tp *A, const inla_cgeneric_smat_tp *B) {
  assert(A->ncol == B->nrow);
  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.nrow = A->nrow;
  res.ncol = B->ncol;

  int max_alloc_size = 5;
  res.x = Calloc(max_alloc_size, double);
  res.i = Calloc(max_alloc_size, int);
  res.j = Calloc(max_alloc_size, int);
  int n = 0;
  
  int_loc_tp A_row_locs = locate_smat_rows(A); // Where can I find the location of row 0, 1, 2, ...?
  int_loc_tp B_col_locs = locate_smat_cols(B); // Where can I find the location of column 0, 1, 2, ...?

  nonzero_index_tp A_nonzero_cols;
  nonzero_index_tp B_nonzero_rows;
  A_nonzero_cols.index = Calloc(A->ncol, int);
  A_nonzero_cols.val = Calloc(A->ncol, double);
  B_nonzero_rows.index = Calloc(B->nrow, int);
  B_nonzero_rows.val = Calloc(B->nrow, double);

  double tmp;
  for (int i = 0; i < res.nrow; i++) {
    if (A_row_locs.n[i] >= 1) { // Are there any nonzero entries in row nr. i of A?
      // Find all the nonzero columns in row nr. i of A, and their values
      // printf("Row number %d: nonzero cols in A:\n", i);
      for (int k = 0; k < A_row_locs.n[i]; k++) {
	A_nonzero_cols.index[k] = A->j[A_row_locs.locs[i][k]];
	A_nonzero_cols.val[k] = A->x[A_row_locs.locs[i][k]];
      }
      A_nonzero_cols.n = A_row_locs.n[i];

      for (int j = 0; j < res.ncol; j++) {
  	if (B_col_locs.n[j] >= 1) { // Are there any nonzero entries in column nr. j of B?
	  // Find all the nonzero rows in column nr. j of B, and their values
	  for (int k = 0; k < B_col_locs.n[j]; k++) {
	    B_nonzero_rows.index[k] = B->i[B_col_locs.locs[j][k]];
	    B_nonzero_rows.val[k] = B->x[B_col_locs.locs[j][k]];
	  }
	  B_nonzero_rows.n = B_col_locs.n[j];

	  // Find all the values k s.t. A[i, k] != 0 and B[k, j] != 0,
	  // then multiply A[i, k] * B[k, j] and sum over k
	  tmp = sum_of_multiplied_intersects(A_nonzero_cols, B_nonzero_rows);

	  // If tmp != 0, then add it to the results
	  if (tmp != 0) {
	    if (n == max_alloc_size) {
	      max_alloc_size *= 2;
	      res.x = realloc(res.x, max_alloc_size * sizeof(double));
	      res.i = realloc(res.i, max_alloc_size * sizeof(int));
	      res.j = realloc(res.j, max_alloc_size * sizeof(int));
	    }
	    res.x[n] = tmp;
	    res.i[n] = i;
	    res.j[n] = j;
	    n += 1;
	  }
  	}
      }
    }
  }
  res.n = n;
  return res;
}

