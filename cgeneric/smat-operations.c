#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#include "cgeneric.h"
#include "smat-operations.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

inla_cgeneric_smat_tp smat_smat_add(inla_cgeneric_smat_tp const * A,
				    inla_cgeneric_smat_tp const * B) {
  // Assert that A and B have equal dimensions
  assert(A->ncol == B->ncol);
  assert(A->nrow == B->nrow);

  // Preallocate the result of A + B
  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.ncol = A->ncol;
  res.nrow = A->nrow;
  int max_alloc_size = int_min(A->n + B->n, res.ncol * res.nrow);
  res.x = Calloc(max_alloc_size, double);
  res.i = Calloc(max_alloc_size, int);
  res.j = Calloc(max_alloc_size, int);
  
  // Keep track of which rows of A and B contain nonzero elements,
  // and how many elements that are nonzero in each row
  int_loc A_row_locs = locate_nonzero_smat_rows(A);
  int_loc B_row_locs = locate_nonzero_smat_rows(B);

  // Preallocate structs for keeping track of all the nonzero elements whithin
  // a certain row of A, B and A + B
  nonzero_index A_nonzero_cols;
  nonzero_index B_nonzero_cols;
  nonzero_index A_B_sum;
  A_nonzero_cols.index = Calloc(A->ncol, int);
  A_nonzero_cols.val = Calloc(A->ncol, double);
  B_nonzero_cols.index = Calloc(B->ncol, int);
  B_nonzero_cols.val = Calloc(B->ncol, double);

  // Loop over all the rows of A and B, and compute A + B
  int count = 0;
  for (int i = 0; i < A->nrow; ++i) {
    if (A_row_locs.n[i] == 0 && B_row_locs.n[i] > 0) {
      // If there are only nonzero elements in the ith row of B,
      // Just set res[, i] = B[, i]
      for (int j = 0; j < B_row_locs.n[i]; ++j) {
	res.x[count] = B->x[B_row_locs.locs[i][j]];
	res.j[count] = B->j[B_row_locs.locs[i][j]];
	res.i[count] = i;
	count += 1;
      }
    } else if (A_row_locs.n[i] > 0 && B_row_locs.n[i] == 0) {
      // If there are only nonzero elements in the ith row of A,
      // Just set res[, i] = A[, i]
      for (int j = 0; j < A_row_locs.n[i]; ++j) {
	res.x[count] = A->x[A_row_locs.locs[i][j]];
	res.j[count] = A->j[A_row_locs.locs[i][j]];
	res.i[count] = i;
	count += 1;
      }
    } else if (A_row_locs.n[i] + B_row_locs.n[i] > 0) {
      // Everything is slightly more complicated if there are nonzero elements
      // in both A[, i] and B[, i].

      // First, loop over the nonzero elements of A and B, and add them into
      // the structs A_nonzero_cols and B_nonzero_cols
      for (int j = 0; j < A_row_locs.n[i]; j++) {
	A_nonzero_cols.index[j] = A->j[A_row_locs.locs[i][j]];
	A_nonzero_cols.val[j] = A->x[A_row_locs.locs[i][j]];
      }
      A_nonzero_cols.n = A_row_locs.n[i];
      for (int j = 0; j < B_row_locs.n[i]; j++) {
	B_nonzero_cols.index[j] = B->j[B_row_locs.locs[i][j]];
	B_nonzero_cols.val[j] = B->x[B_row_locs.locs[i][j]];
      }
      B_nonzero_cols.n = B_row_locs.n[i];

      // Then, compute the sum of A[, i] + B[, i]
      A_B_sum = sum_of_nonzero_indices(&A_nonzero_cols, &B_nonzero_cols);

      // Then, fill the result into res[, i]
      for (int j = 0; j < A_B_sum.n; j++) {
	res.x[count] = A_B_sum.val[j];
	res.j[count] = A_B_sum.index[j];
	res.i[count] = i;
	count += 1;
      }

      // Clean up after your self
      free_nonzero_index(A_B_sum);
    }
  }

  // Update the number of elements in res
  res.n = count;

  // Reallocate all the elements of res to free unused memory
  res.i = (int *) realloc(res.i, res.n * sizeof(int));
  res.j = (int *) realloc(res.j, res.n * sizeof(int));
  res.x = (double *) realloc(res.x, res.n * sizeof(double));

  // Clean up after your self
  free_nonzero_index(A_nonzero_cols);
  free_nonzero_index(B_nonzero_cols);

  return res;
}

void smat_smat_add_in_place(inla_cgeneric_smat_tp * A,
			    inla_cgeneric_smat_tp const * B) {
  inla_cgeneric_smat_tp res = smat_smat_add(A, B);
  free_smat(A);
  A->name = res.name;
  A->i = res.i;
  A->j = res.j;
  A->x = res.x;
  A->n = res.n;
}

void diag_smat_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D) {
  assert(D->dim == A->nrow);
  for (int i = 0; i < A->n; i++) {
    A->x[i] *= D->x[A->i[i]];
  }
}

void smat_diag_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D) {
  assert(D->dim == A->ncol);
  for (int i = 0; i < A->n; i++) {
    A->x[i] *= D->x[A->j[i]];
  }
}

void diag_smat_diag_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D) {
  assert(D->dim == A->ncol);
  assert(D->dim == A->nrow);
  for (int i = 0; i < A->n; i++) {
    A->x[i] *= (D->x[A->j[i]] * D->x[A->i[i]]);
  }
}


void transpose_smat(inla_cgeneric_smat_tp * A) {
  int * tmp_pnt = A->i;
  A->i = A->j;
  A->j = tmp_pnt;
  int tmp = A->nrow;
  A->nrow = A->ncol;
  A->ncol = tmp;
}

inla_cgeneric_smat_tp copy_smat(inla_cgeneric_smat_tp const * A) {
  inla_cgeneric_smat_tp res;
  res.nrow = A->nrow;
  res.ncol = A->ncol;
  res.n = A->n;
  res.x = Calloc(A->n, double);
  res.i = Calloc(A->n, int);
  res.j = Calloc(A->n, int);
  for (int k = 0; k < A->n; k++) {
    res.x[k] = A->x[k];
    res.i[k] = A->i[k];
    res.j[k] = A->j[k];
  }
  // We never use these names for anything, and they are just a hassle to deal with,
  // so we just set it equal to NULL
  res.name = NULL;
  return res;
}

inla_cgeneric_mat_tp copy_mat(const inla_cgeneric_mat_tp *A) {
  inla_cgeneric_mat_tp res;
  res.nrow = A->nrow;
  res.ncol = A->ncol;
  res.x = Calloc(A->nrow * A->ncol, double);
  for (int i = 0; i < A->nrow * A->ncol; i++) {
    res.x[i] = A->x[i];
  }
  // We never use these names for anything, and they are just a hassle to deal with,
  // so we just set it equal to NULL
  res.name = NULL;
  return res;
}

void free_smat(inla_cgeneric_smat_tp *A) {
  free(A->x);
  free(A->i);
  free(A->j);
  free(A->name);
}

void free_mat(inla_cgeneric_mat_tp *A) {
  free(A->x);
  free(A->name);
}

nonzero_index sum_of_nonzero_indices(nonzero_index const * x,
				     nonzero_index const * y) {
  // Allocate the result of x + y.
  // We assume that there are no duplicate elements inside each of the index-structs
  nonzero_index res;
  res.index = Calloc(x->n + y->n, int);
  res.val = Calloc(x->n + y->n, double);

  // Loop over the nonzero elements of x and add them to res
  int count = 0;
  for (int i = 0; i < x->n; ++i) {
    res.index[count] = x->index[i];
    res.val[count] = x->val[i];
    for (int j = 0; j < y->n; ++j) {
      // Loop over y and see if we have a nonzero element in the same
      // location as x->index[i]. If this is true, then add that element
      // to res.
      if (x->index[i] == y->index[j]) {
	res.val[count] += y->val[j];
	break;
      }
    }
    count += 1;
  }

  // Loop over the nonzero elements of y, check if they have already been added to res,
  // then add them if they have not already been added.
  int any_equal;
  for (int i = 0; i < y->n; i++) {
    any_equal = 0;
    for (int j = 0; j < x->n; j++) {
      if (y->index[i] == x->index[j]) {
	any_equal = 1;
	break;
      }
    }
    if (any_equal == 0) {
      res.index[count] = y->index[i];
      res.val[count] = y->val[i];
      count += 1;
    }
  }
  res.n = count;

  // Reallocate all the elements of res to free unused memory
  res.index = (int *) realloc(res.index, res.n * sizeof(int));
  res.val = (double *) realloc(res.val, res.n * sizeof(double));

  return res;
}

void free_nonzero_index(nonzero_index A) {
  free(A.index);
  free(A.val);
}

int_loc locate_nonzero_smat_rows(inla_cgeneric_smat_tp const * A) {
  return locate_ints_in_arr(A->i, A->nrow, A->n);
}

int_loc locate_ints_in_arr(int const * arr, int num_ints, int arr_length) {
  // Allocate the result
  int_loc res;
  res.locs = (int **) calloc(num_ints, sizeof(int *));
  res.n = Calloc(num_ints, int);
  res.num_ints = num_ints;

  // allocation_lengths is a list of the allocation lengths for each of the
  // num_ints arrays in res.locs
  int * allocation_lengths = Calloc(num_ints, int);
  for (int i = 0; i < num_ints; ++i) {
    res.locs[i] = Calloc(5, int);
    allocation_lengths[i] = 5;
    res.n[i] = 0;
  }

  // Loop over the array arr and record which int is found at arr[i]
  // for i = 0, ..., arr_length
  for (int i = 0; i < arr_length; ++i) {
    // Reallocate to increase the size of res.locs[arr[i]] if it becomes necessary
    if ((res.n[arr[i]] + 1) == allocation_lengths[arr[i]]) {
      allocation_lengths[arr[i]] *= 2;
      res.locs[arr[i]] = realloc(res.locs[arr[i]], allocation_lengths[arr[i]] * sizeof(int));
    }
    res.locs[arr[i]][res.n[arr[i]]] = i;
    res.n[arr[i]] += 1;
  }

  // Reallocate all arrays in res.locs to get the correct lengths and free
  // all unused memory
  for (int i = 0; i < res.num_ints; ++i) {
    res.locs[i] = realloc(res.locs[i], res.n[i] * sizeof(int));
  }

  // Clean up after yourself
  free(allocation_lengths);

  return res;
}

int int_min(int x, int y) {
  // Return x if y > x and y if y <= x
  return y > x ? x : y;
}
