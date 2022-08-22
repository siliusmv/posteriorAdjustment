
#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"
#include "smat-operations.h"
#include "spde-precision.h"


#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

double *inla_cgeneric_spde_model_with_b_func(inla_cgeneric_cmd_tp cmd,
					     double *theta,
					     inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Number of observations
  assert(!strcasecmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Which s0 do the different replications come from?
  assert(!strcasecmp(data->ints[2]->name, "s0_index"));
  int n_repl = data->ints[2]->len;
  int *s0_index = data->ints[2]->ints;

  // B0, B1 and B2, required for building the SPDE precision
  assert(!strcasecmp(data->mats[0]->name, "B0"));
  int n_spde = data->mats[0]->nrow;
  assert(n == (n_repl) * n_spde);
  assert(data->mats[0]->ncol == 3);
  assert(!strcasecmp(data->mats[1]->name, "B1"));
  assert(data->mats[1]->nrow == n_spde);
  assert(data->mats[1]->ncol == 3);
  assert(!strcasecmp(data->mats[2]->name, "B2"));
  assert(data->mats[2]->nrow == n_spde);
  assert(data->mats[2]->ncol == 3);

  // Distance to all conditioning sites from all mesh nodes
  assert(!strcasecmp(data->mats[3]->name, "dist_to_s0"));
  inla_cgeneric_mat_tp *dist_to_s0 = data->mats[3];

  // M0, M1 and M2, required for building the SPDE precision
  assert(!strcasecmp(data->smats[0]->name, "M0"));
  assert(data->smats[0]->nrow == n_spde);
  assert(data->smats[0]->ncol == n_spde);
  assert(!strcasecmp(data->smats[1]->name, "M1"));
  assert(data->smats[1]->nrow == n_spde);
  assert(data->smats[1]->ncol == n_spde);
  assert(!strcasecmp(data->smats[2]->name, "M2"));
  assert(data->smats[2]->nrow == n_spde);
  assert(data->smats[2]->ncol == n_spde);

  // Initial values and priors
  assert(!strcasecmp(data->doubles[0]->name, "init"));
  double *init = data->doubles[0]->doubles;
  int n_theta = data->doubles[0]->len;
  assert(!strcasecmp(data->doubles[1]->name, "rho_prior"));
  double *rho_prior = data->doubles[1]->doubles;
  assert(!strcasecmp(data->doubles[2]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[2]->doubles;
  assert(!strcasecmp(data->doubles[3]->name, "rho_b_prior"));
  double *rho_b_prior = data->doubles[3]->doubles;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, log_rho_b, rho_b;
  if (theta) {
    int count_theta = 0;
    if (rho_prior[0] == 0) {
      log_rho = rho_prior[1];
    } else {
      log_rho = theta[count_theta];
      count_theta += 1;
    }
    if (sigma_prior[0] == 0) {
      log_sigma = sigma_prior[1];
    } else {
      log_sigma = theta[count_theta];
      count_theta += 1;
    }
    if (rho_b_prior[0] == 0) {
      log_rho_b = rho_b_prior[1];
    } else {
      log_rho_b = theta[count_theta];
      count_theta += 1;
    }
    rho_b = exp(log_rho_b);
    assert(count_theta == n_theta);
  } else {
    log_rho = log_sigma = log_rho_b = rho_b = NAN;
  }
 
  // ==========================================================
  // This switch statement is the required method for implementing
  // cgeneric models in R-INLA
  // ==========================================================
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // return a vector of indices with format
      // c(n, M, ii, jj)
      // where ii<=jj and both ii and jj are non-decreasing
      // and M is the length of ii

      // Here we just compute the precision matrix, in order to find the graph structure
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);

      // Only keep the upper diagonal of the precision matrix
      upper_diag(&precision);

      // Sort the entries of the smat so they correspond to the R-INLA requirements
      sort_smat(&precision);

      // Replicate the precision matrix several times to create a block diagonal matrix
      block_diag_smat(&precision, n_repl);

      // Extract all the necessary i and j values from the block matrix
      ret = Calloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  /* i */
	ret[2 + precision.n + i] = precision.j[i];                    /* j */
      }

      // Free the precision matrix so we don't get memory leaks
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute the precision matrix
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);

      // Only keep the upper diagonal of the precision matrix
      upper_diag(&precision);

      // Sort the entries of the smat so they correspond to the R-INLA requirements
      sort_smat(&precision);

      // Replicate the precision matrix several times to create a block diagonal matrix
      block_diag_smat(&precision, n_repl);

      // Compute the values of b for each distance from a conditioning site
      inla_cgeneric_mat_tp b_vals = copy_mat(dist_to_s0);
      double tmp;
      for (int i = 0; i < b_vals.nrow * b_vals.ncol; i++) {
	// Ensure that we don't get numerical problems where b_vals.x[i] becomes INF or NAN
	tmp = dist_to_s0->x[i] / rho_b;
	if (tmp < 0.0000000001) {
	  tmp = 0.0000000001;
	}
	b_vals.x[i] = 1 / sqrt(1 - exp(-2 * tmp));
      }

      // Create a diagonal matrix containing the values of 1 / b, that can be used
      // for rescaling the precision matrix with the correct b values
      diag_mat_tp b_inv;
      b_inv.x = Calloc(precision.nrow, double);
      b_inv.dim = precision.nrow;
      for (int i = 0; i < n_repl; i++) {
	for (int j = 0; j < n_spde; j++) {
	  // use (s0_index[i] - 1) because R uses 1-indexing and C uses 0-indexing
	  b_inv.x[j + i * n_spde] = b_vals.x[j + (s0_index[i] - 1) * n_spde];
	}
      }
      // Multiply the diagonal matrix with the precision matrix: B_inv * Q * B_inv
      diag_smat_diag_mult(&precision, &b_inv);

      // Avoid memory leaks
      free(b_inv.x);
      free_mat(&b_vals);

      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij
      ret = Calloc(2 + precision.n, double);
      ret[0] = -1;                                   /* code for optimized output */
      ret[1] = precision.n;                          /* number of (i <= j) */
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.x[i];
      }
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0

      // Z_b has zero mean
      ret = Calloc(1, double);
      ret[0] = 0;
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters

      // The initial values depend on how many parameters that are fixed or estimated.
      ret = Calloc(n_theta + 1, double);
      ret[0] = n_theta;
      for (int i = 0; i < n_theta; i++) {
	ret[i + 1] = init[i];
      }
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

      // This is too difficult to compute for us
      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      // Only add terms for the parameters that are not fixed
      ret = Calloc(1, double);
      ret[0] = 0;
      if (rho_prior[0] != 0) {
	double lambda0 = -log(rho_prior[2]) * rho_prior[1];
	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
      }
      if (sigma_prior[0] != 0) {
	double lambda1 = -log(sigma_prior[2]) / sigma_prior[1];
	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
      }
      if (rho_b_prior[0] != 0) {
	double gaussian_const = -0.91893853320467;
	ret[0] += gaussian_const - log(rho_b_prior[2]) -
	  pow(log_rho_b - rho_b_prior[1], 2) / (2 * pow(rho_b_prior[2], 2));
      } 
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double *inla_cgeneric_iid_model_with_b_func(inla_cgeneric_cmd_tp cmd,
					    double *theta,
					    inla_cgeneric_data_tp * data) {
  // This is a special version of Z_b, where we assume that the process Z
  // is simply Gaussian white noise, meaning that we do not need to
  // provide an SPDE approximation, and that the precision matrix is diagonal.
  double *ret = NULL;

  // Number of observations
  assert(!strcasecmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  assert(!strcasecmp(data->doubles[0]->name, "dist_to_s0"));
  assert(data->doubles[0]->len == n);
  double *dist = data->doubles[0]->doubles;

  assert(!strcasecmp(data->doubles[1]->name, "init"));
  double *init = data->doubles[1]->doubles;
  int n_theta = data->doubles[1]->len;
  assert(!strcasecmp(data->doubles[2]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[2]->doubles;
  assert(!strcasecmp(data->doubles[3]->name, "rho_b_prior"));
  double *rho_b_prior = data->doubles[3]->doubles;

  double log_sigma, sigma, log_rho_b, rho_b;
  if (theta) {
    int count_theta = 0;
    if (sigma_prior[0] == 0) {
      log_sigma = sigma_prior[1];
    } else {
      log_sigma = theta[count_theta];
      count_theta += 1;
    }
    sigma = exp(log_sigma);
    if (rho_b_prior[0] == 0) {
      log_rho_b = rho_b_prior[1];
    } else {
      log_rho_b = theta[count_theta];
      count_theta += 1;
    }
    rho_b = exp(log_rho_b);
    assert(count_theta == n_theta);
  } else {
    sigma = log_sigma = log_rho_b = rho_b = NAN;
  }
 
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // return a vector of indices with format
      // c(n, M, ii, jj)
      // where ii<=jj and both ii and jj are non-decreasing
      // and M is the length of ii
      ret = Calloc(2 + 2 * n, double);
      ret[0] = n;
      ret[1] = n;
      for (int i = 0; i < n; i++) {
      	ret[2 + i] = i;			       /* i */
      	ret[2 + n + i] = i;		       /* j */
      }
      break;
    }

  case INLA_CGENERIC_Q:
    {
      ret = Calloc(2 + n, double);
      ret[0] = -1;			       /* code for optimized output */
      ret[1] = n;			       /* number of (i <= j) */
      double tmp;
      for (int i = 0; i < n; i++) {
	tmp = dist[i] / rho_b;
	if (tmp < 0.0000000001) {
	  tmp = 0.0000000001;
	}
	ret[2 + i] = 1 / (sigma * sigma * (1 - exp(-2 * tmp)));
      }
      break;
    }
  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0
      ret = Calloc(1, double);
      ret[0] = 0;
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters
      ret = Calloc(n_theta + 1, double);
      ret[0] = n_theta;
      for (int i = 0; i < n_theta; i++) {
	ret[i + 1] = init[i];
      }
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)
      ret = Calloc(1, double);
      ret[0] = 0;
      double gaussian_const = -0.91893853320467;
      if (sigma_prior[0] != 0) {
	ret[0] += gaussian_const - log(sigma_prior[2]) - pow(log_sigma - sigma_prior[1], 2) / (2 * pow(sigma_prior[2], 2));
      } 
      if (rho_b_prior[0] != 0) {
	ret[0] += gaussian_const - log(rho_b_prior[2]) - pow(log_rho_b - rho_b_prior[1], 2) / (2 * pow(rho_b_prior[2], 2));
      } 
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

