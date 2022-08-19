#pragma once

#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

// Create a deep copy of an inla_cgeneric_mat_tp
inla_cgeneric_mat_tp copy_mat(inla_cgeneric_mat_tp const * A) {
  inla_cgeneric_mat_tp res;
  res.name = NULL;
  res.nrow = A->nrow;
  res.ncol = A->ncol;
  res.x = Calloc(A->nrow * A->ncol, double);
  for (int i = 0; i < A->nrow * A->ncol; ++i) {
    res.x[i] = A->x[i];
  }
  return res;
}

// Free memory from an inla_cgeneric_smat_tp
void free_smat(inla_cgeneric_smat_tp * A) {
  free(A->x);
  free(A->i);
  free(A->j);
  free(A->name);
}

// Free memory from an inla_cgeneric_mat_tp
void free_mat(inla_cgeneric_mat_tp * A) {
  free(A->x);
  free(A->name);
}

// This is a struct that describes a diagonal matrix of dimension (dim x dim),
// and with diagonal elements x
typedef struct {
  int dim;
  double *x;
} diag_mat_tp;

// Given an inla_cgeneric_smat_tp A and a diagonal matrix D,
// Compute D * A * D
void diag_smat_diag_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D) {
  assert(D->dim == A->nrow);
  assert(D->dim == A->ncol);
  for (int i = 0; i < A->n; ++i) {
    A->x[i] *= (D->x[A->i[i]] * D->x[A->j[i]]);
  }
}

// Remove all the lower diagonal elements of an inla_cgeneric_smat_tp
void upper_diag(inla_cgeneric_smat_tp * A) {
  int count = 0;
  for (int i = 0; i < A->n; ++i) {
    if (A->i[i] <= A->j[i]) {
      if (count < i) {
	A->i[count] = A->i[i];
	A->j[count] = A->j[i];
	A->x[count] = A->x[i];
      }
      count += 1;
    }
  }
  A->n = count;
}

// This is a struct for containing one element of an inla_cgeneric_smat_tp.
// This struct is used for sorting the elements of an inla_cgeneric_smat_tp,
// in the function sort_smat()
typedef struct {
  int i;
  int j;
  double x;
} smat_element;

// This is a function used as input to the qsort() function, which we use
// for sorting the elements of an inla_cgeneric_smat_tp, in the function sort_smat()
int smat_cmp_func_i_then_j(const void * a, const void * b) {
  if (((smat_element*)a)->i > ((smat_element*)b)->i) {
    return 1;
  } else if (((smat_element*)a)->i < ((smat_element*)b)->i) {
    return -1;
  } else if (((smat_element*)a)->j > ((smat_element*)b)->j) {
    return 1;
  } else if (((smat_element*)a)->j < ((smat_element*)b)->j) {
    return -1;
  } else {
    return 0;
  }
}

// Sort the elements of an inla_cgeneric_smat_tp, such that i is strictly increasing,
// while j is strictly increasing within each value of i.
void sort_smat(inla_cgeneric_smat_tp * A) {
  smat_element * elements = Calloc(A->n, smat_element);
  for (int i = 0; i < A->n; ++i) {
    elements[i].i = A->i[i];
    elements[i].j = A->j[i];
    elements[i].x = A->x[i];
  }
  qsort(elements, A->n, sizeof(smat_element), smat_cmp_func_i_then_j);
  for (int i = 0; i < A->n; ++i) {
    A->i[i] = elements[i].i;
    A->j[i] = elements[i].j;
    A->x[i] = elements[i].x;
  }
  free(elements);
}

// Given an inla_cgeneric_smat_tp A, create a block diagonal matrix
// that contains A n times.
void block_diag_smat(inla_cgeneric_smat_tp * A, int n) {
  A->i = (int *) realloc(A->i, A->n * n * sizeof(int));
  A->j = (int *) realloc(A->j, A->n * n * sizeof(int));
  A->x = (double *) realloc(A->x, A->n * n * sizeof(double));
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < A->n; ++j) {
      A->i[j + A->n * i] = A->i[j] + A->nrow * i;
      A->j[j + A->n * i] = A->j[j] + A->ncol * i;
      A->x[j + A->n * i] = A->x[j];
    }
  }
  A->nrow *= n;
  A->ncol *= n;
  A->n *= n;
}

// Create a gsl_spmatrix that is a copy of an inla_cgeneric_smat_tp matrix
gsl_spmatrix * inla_smat_to_gsl_spmatrix(inla_cgeneric_smat_tp const * A) {
  gsl_spmatrix * tmp = gsl_spmatrix_alloc(A->nrow, A->ncol);
  for (int i = 0; i < A->n; ++i) {
    gsl_spmatrix_set(tmp, A->i[i], A->j[i], A->x[i]);
  }
  gsl_spmatrix * res = gsl_spmatrix_ccs(tmp);
  gsl_spmatrix_free(tmp);
  return res;
}

// Create a copy of a gsl_spmatrix A that is on the CSC format, where the
// new copy is on the COO format
gsl_spmatrix * csc_to_coo(gsl_spmatrix const * A) {
  gsl_spmatrix *res = gsl_spmatrix_alloc(A->size1, A->size2);
  int counter = 0, nums_in_col;
  for (int i = 0; i < A->size2; ++i) {
    if (i < (A->size2 - 1)) {
      nums_in_col = A->p[i + 1] - A->p[i];
    } else {
      nums_in_col = A->nz - A->p[i];
    }
    for (int j = 0; j < nums_in_col; ++j) {
      gsl_spmatrix_set(res, A->i[counter], i, A->data[counter]);
      counter += 1;
    }
  }
  return res;
}

// Create an inla_cgeneric_smat_tp matrix that is a copy of a gsl_spmatrix
inla_cgeneric_smat_tp gsl_spmatrix_to_inla_smat(gsl_spmatrix * A) {
  gsl_spmatrix * tmp;
  if (A->sptype == GSL_SPMATRIX_CSC) {
    tmp = csc_to_coo(A);
  } else {
    tmp = A;
  }
  assert(tmp->sptype == GSL_SPMATRIX_COO);
  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.nrow = tmp->size1;
  res.ncol = tmp->size2;
  res.n = tmp->nz;
  res.x = Calloc(res.n, double);
  res.i = Calloc(res.n, int);
  res.j = Calloc(res.n, int);
  for (int k = 0; k < tmp->nz; ++k) {
    res.x[k] = tmp->data[k];
    res.i[k] = tmp->i[k];
    res.j[k] = tmp->p[k];
  }
  if (A->sptype != GSL_SPMATRIX_COO) {
    gsl_spmatrix_free(tmp);
  }
  return res;
}

// Copy of the R functions INLA::theta2phiX(), where X = 0, 1, 2.
// Necessary for computing the precision matrix of the SPDE approximation.
gsl_vector * theta2phi(double log_rho,
		       double log_sigma,
		       inla_cgeneric_mat_tp const * B,
		       int exp_transform) {
  assert(B->ncol == 3);
  gsl_vector * res = gsl_vector_alloc(B->nrow);
  double tmp;
  for (int i = 0; i < B->nrow; ++i) {
    tmp = B->x[i * B->ncol] +
      B->x[i * B->ncol + 1] * log_rho +
      B->x[i * B->ncol + 2] * log_sigma;
    if (exp_transform) {
      tmp = exp(tmp);
    }
    gsl_vector_set(res, i, tmp);
  }
  return res;
}

// Given two gsl_spmatrices left and right, set the matrix left equal to the value
// of left + right.
void gsl_spmatrix_add_in_place(gsl_spmatrix * left, gsl_spmatrix const * right) {
  assert(left->size1 == right->size1);
  assert(left->size2 == right->size2);
  assert(left->sptype == right->sptype);
  gsl_spmatrix * tmp = gsl_spmatrix_alloc(left->size1, left->size2);
  tmp->sptype = left->sptype;
  gsl_spmatrix_memcpy(tmp, left);
  gsl_spmatrix_add(left, tmp, right);
  gsl_spmatrix_free(tmp);
}

// Compute the precision matrix of the SPDE approximation.
// This is a copy of the R function INLA::inla.spde2.precision()
inla_cgeneric_smat_tp spde_precision(double log_rho,
				     double log_sigma,
				     inla_cgeneric_mat_tp** B,
				     inla_cgeneric_smat_tp** M) {
  gsl_vector * phi0 = theta2phi(log_rho, log_sigma, B[0], 1);
  gsl_vector * phi1 = theta2phi(log_rho, log_sigma, B[1], 1);
  gsl_vector * phi12 = theta2phi(log_rho, log_sigma, B[2], 0);
  gsl_vector_mul(phi12, phi1);

  // Compute D1 * M0 * D1, where D1 is a diagonal matrix containing phi0,
  // seen in the R function INLA::inla.spde2.precision()
  gsl_spmatrix * tmp1 = inla_smat_to_gsl_spmatrix(M[0]);
  gsl_spmatrix_scale_columns(tmp1, phi1);
  gsl_spmatrix_scale_rows(tmp1, phi1);

  // Compute D12 * M1, where D12 is a diagonal matrix containing phi12,
  // seen in the R function INLA::inla.spde2.precision()
  gsl_spmatrix * tmp2 = inla_smat_to_gsl_spmatrix(M[1]);
  gsl_spmatrix_scale_rows(tmp2, phi12);

  // Compute the temporary precision matrix Q = D1 * M0 * D1 + D12 * M1
  gsl_spmatrix * Q = gsl_spmatrix_alloc(M[0]->nrow, M[0]->ncol);
  Q->sptype = GSL_SPMATRIX_CSC;
  gsl_spmatrix_add(Q, tmp1, tmp2);
  gsl_spmatrix_free(tmp1);

  // Transpose D12 * M1 to M1' * D12
  tmp1 = gsl_spmatrix_alloc(tmp2->size1, tmp2->size2);
  tmp1->sptype = tmp2->sptype;
  gsl_spmatrix_transpose_memcpy(tmp1, tmp2);
  gsl_spmatrix_free(tmp2);

  // Compute Q = Q + M1' * D12
  gsl_spmatrix_add_in_place(Q, tmp1);
  gsl_spmatrix_free(tmp1);

  // Compute Q = Q + M2
  tmp1 = inla_smat_to_gsl_spmatrix(M[2]);
  gsl_spmatrix_add_in_place(Q, tmp1);
  gsl_spmatrix_free(tmp1);

  // Compute Q = D0 * Q * D0
  gsl_spmatrix_scale_columns(Q, phi0);
  gsl_spmatrix_scale_rows(Q, phi0);

  // Avoid memory leaks
  gsl_vector_free(phi0);
  gsl_vector_free(phi1);
  gsl_vector_free(phi12);

  // Convert the precision matrix to inla_cgeneric_smat_tp format
  inla_cgeneric_smat_tp res = gsl_spmatrix_to_inla_smat(Q);
  gsl_spmatrix_free(Q);

  return res;
}
