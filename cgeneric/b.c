
#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <omp.h>

#include "cgeneric.h"
#include "spde-funcs.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

double *inla_cgeneric_iid_model_with_b_func(inla_cgeneric_cmd_tp cmd,
					    double *theta,
					    inla_cgeneric_data_tp * data) {

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
  assert(!strcasecmp(data->doubles[3]->name, "rho_2_prior"));
  double *rho_2_prior = data->doubles[3]->doubles;

  double log_sigma, sigma, log_rho_2, rho_2;
  if (theta) {
    int count_theta = 0;
    if (sigma_prior[0] == 0) {
      log_sigma = sigma_prior[1];
    } else {
      log_sigma = theta[count_theta];
      count_theta += 1;
    }
    sigma = exp(log_sigma);
    if (rho_2_prior[0] == 0) {
      log_rho_2 = rho_2_prior[1];
    } else {
      log_rho_2 = theta[count_theta];
      count_theta += 1;
    }
    rho_2 = exp(log_rho_2);
    assert(count_theta == n_theta);
  } else {
    sigma = log_sigma = log_rho_2 = rho_2 = NAN;
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
	tmp = dist[i] / rho_2;
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
      if (rho_2_prior[0] != 0) {
	ret[0] += gaussian_const - log(rho_2_prior[2]) - pow(log_rho_2 - rho_2_prior[1], 2) / (2 * pow(rho_2_prior[2], 2));
      } 
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double *inla_cgeneric_spde_model_with_b_func(inla_cgeneric_cmd_tp cmd,
					     double *theta,
					     inla_cgeneric_data_tp * data) {

  double *ret = NULL;

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

  assert(!strcasecmp(data->doubles[0]->name, "init"));
  double *init = data->doubles[0]->doubles;
  int n_theta = data->doubles[0]->len;
  assert(!strcasecmp(data->doubles[1]->name, "rho_prior"));
  double *rho_prior = data->doubles[1]->doubles;
  assert(!strcasecmp(data->doubles[2]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[2]->doubles;
  assert(!strcasecmp(data->doubles[3]->name, "rho_2_prior"));
  double *rho_2_prior = data->doubles[3]->doubles;

  double log_rho, log_sigma, log_rho_2, rho_2;
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
    if (rho_2_prior[0] == 0) {
      log_rho_2 = rho_2_prior[1];
    } else {
      log_rho_2 = theta[count_theta];
      count_theta += 1;
    }
    rho_2 = exp(log_rho_2);
    assert(count_theta == n_theta);
  } else {
    log_rho = log_sigma = log_rho_2 = rho_2 = NAN;
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
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);
      block_diag_smat(&precision, n_repl);
      ret = Calloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  /* i */
	ret[2 + precision.n + i] = precision.j[i];                    /* j */
      }
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);
      block_diag_smat(&precision, n_repl);

      inla_cgeneric_mat_tp beta = copy_mat(dist_to_s0);
      double tmp;
      for (int i = 0; i < beta.nrow * beta.ncol; i++) {
	// Ensure that we don't get numerical problems where beta.x[i] becomes INF or NAN
	tmp = dist_to_s0->x[i] / rho_2;
	if (tmp < 0.0000000001) {
	  tmp = 0.0000000001;
	}
	beta.x[i] = 1 / sqrt(1 - exp(-2 * tmp));
      }

      diag_mat_tp b_inv;
      b_inv.x = Calloc(precision.nrow, double);
      b_inv.dim = precision.nrow;
      for (int i = 0; i < n_repl; i++) {
	for (int j = 0; j < n_spde; j++) {
	  // use (s0_index[i] - 1) because R uses 1-indexing and C uses 0-indexing
	  b_inv.x[j + i * n_spde] = beta.x[j + (s0_index[i] - 1) * n_spde];
	}
      }
      diag_smat_mult(&precision, &b_inv);
      smat_diag_mult(&precision, &b_inv);
      free(b_inv.x);
      free_mat(&beta);

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
      if (rho_prior[0] != 0) {
	double lambda0 = -log(rho_prior[2]) * rho_prior[1];
	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
      }
      if (sigma_prior[0] != 0) {
	double lambda1 = -log(sigma_prior[2]) / sigma_prior[1];
	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
      }
      if (rho_2_prior[0] != 0) {
	double gaussian_const = -0.91893853320467;
	ret[0] += gaussian_const - log(rho_2_prior[2]) - pow(log_rho_2 - rho_2_prior[1], 2) / (2 * pow(rho_2_prior[2], 2));
      } 
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double *inla_cgeneric_spde_model_with_replicates(inla_cgeneric_cmd_tp cmd,
						 double *theta,
						 inla_cgeneric_data_tp * data) {
  double *ret = NULL, log_rho, log_sigma;

  if (theta) {
    log_rho = theta[0];
    log_sigma = theta[1];
  } else {
    log_rho = log_sigma = NAN;
  }
  int n_theta = 2;

  // Number of observations
  assert(!strcasecmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Number of replicates
  assert(!strcasecmp(data->ints[2]->name, "n_repl"));
  int n_repl = data->ints[2]->ints[0];
  assert(n_repl > 0);

  // B0, B1 and B2, required for building the SPDE precision
  assert(!strcasecmp(data->mats[0]->name, "B0"));
  int n_spde = data->mats[0]->nrow;
  assert(data->mats[0]->ncol == 3);
  assert(!strcasecmp(data->mats[1]->name, "B1"));
  assert(data->mats[1]->nrow == n_spde);
  assert(data->mats[1]->ncol == 3);
  assert(!strcasecmp(data->mats[2]->name, "B2"));
  assert(data->mats[2]->nrow == n_spde);
  assert(data->mats[2]->ncol == 3);

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

  assert(!strcasecmp(data->doubles[0]->name, "init"));
  double *init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);
  assert(!strcasecmp(data->doubles[1]->name, "rho_prior"));
  double *rho_prior = data->doubles[1]->doubles;
  assert(!strcasecmp(data->doubles[2]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[2]->doubles;

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
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);
      block_diag_smat(&precision, n_repl);
      ret = Calloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  /* i */
	ret[2 + precision.n + i] = precision.j[i];                    /* j */
      }
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);
      block_diag_smat(&precision, n_repl);

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
      double lambda0, lambda1;
      ret = Calloc(1, double);
      lambda0 = -log(rho_prior[1]) * rho_prior[0];
      lambda1 = -log(sigma_prior[1]) / sigma_prior[0];
      ret[0] = log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
      ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double *inla_cgeneric_spde_model(inla_cgeneric_cmd_tp cmd,
				 double *theta,
				 inla_cgeneric_data_tp * data) {
  double *ret = NULL, log_rho, log_sigma;

  if (theta) {
    log_rho = theta[0];
    log_sigma = theta[1];
  } else {
    log_rho = log_sigma = NAN;
  }
  int n_theta = 2;

  // Number of observations
  assert(!strcasecmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // B0, B1 and B2, required for building the SPDE precision
  assert(!strcasecmp(data->mats[0]->name, "B0"));
  int n_spde = data->mats[0]->nrow;
  assert(data->mats[0]->ncol == 3);
  assert(!strcasecmp(data->mats[1]->name, "B1"));
  assert(data->mats[1]->nrow == n_spde);
  assert(data->mats[1]->ncol == 3);
  assert(!strcasecmp(data->mats[2]->name, "B2"));
  assert(data->mats[2]->nrow == n_spde);
  assert(data->mats[2]->ncol == 3);

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

  assert(!strcasecmp(data->doubles[0]->name, "init"));
  double *init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);
  assert(!strcasecmp(data->doubles[1]->name, "rho_prior"));
  double *rho_prior = data->doubles[1]->doubles;
  assert(!strcasecmp(data->doubles[2]->name, "sigma_prior"));
  double *sigma_prior = data->doubles[2]->doubles;

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
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);
      ret = Calloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  /* i */
	ret[2 + precision.n + i] = precision.j[i];                    /* j */
      }
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);
      upper_diag(&precision);
      sort_smat(&precision);

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
      double lambda0, lambda1;
      ret = Calloc(1, double);
      lambda0 = -log(rho_prior[1]) * rho_prior[0];
      lambda1 = -log(sigma_prior[1]) / sigma_prior[0];
      ret[0] = log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
      ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

