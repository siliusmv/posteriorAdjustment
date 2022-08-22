
#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))

double *inla_cgeneric_a_model(inla_cgeneric_cmd_tp cmd,
			      double *theta,
			      inla_cgeneric_data_tp * data) {
  double *ret = NULL;

  // This is an approximation for -0.5 * log(2π)
  double gaussian_const = -0.91893853320467;
  // The fixed precision used in the precision matrix of a
  double high_prec = exp(15);

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================
  assert(!strcasecmp(data->doubles[0]->name, "y0"));
  int n = data->doubles[0]->len;
  double *y0 = data->doubles[0]->doubles;

  assert(!strcasecmp(data->doubles[1]->name, "dist_to_s0"));
  assert(data->doubles[1]->len == n);
  double *dist = data->doubles[1]->doubles;

  // Initial values for the hyperparameters
  assert(!strcasecmp(data->doubles[2]->name, "init"));
  double *init = data->doubles[2]->doubles;
  double n_theta = data->doubles[2]->len;

  // PC priors for lambda and kappa
  assert(!strcasecmp(data->doubles[3]->name, "lambda_prior"));
  double *lambda_prior = data->doubles[3]->doubles;
  assert(!strcasecmp(data->doubles[4]->name, "kappa_prior"));
  double *kappa_prior = data->doubles[4]->doubles;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double lambda, kappa;
  if (theta) {
    int count_theta = 0;
    if (lambda_prior[0] == 0) {
      lambda = exp(lambda_prior[1]);
    } else {
      lambda = exp(theta[count_theta]);
      count_theta += 1;
    }
    if (kappa_prior[0] == 0) {
      kappa = exp(kappa_prior[1]);
    } else {
      kappa = exp(theta[count_theta]);
      count_theta += 1;
    }
    assert(count_theta == n_theta);
  } else {
    kappa = lambda = NAN;
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

      // The precision matrix is just a diagonal matrix with each element equal to high_prec
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
      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij

      // The precision matrix is just a diagonal matrix with each element equal to high_prec
      ret = Calloc(2 + n, double);
      ret[0] = -1;			       /* code for optimized output */
      ret[1] = n;			       /* number of (i <= j) */
      for (int i = 0; i < n; i++) {
	ret[2 + i] = high_prec;
      }
      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (n, mu)
      // if n==0 then mu is not needed as its taken to be mu[]==0
      ret = Calloc(1 + n, double);
      ret[0] = n;
      for (int i = 0; i < n; i++) {
	ret[1 + i] = y0[i] * exp(-pow(dist[i] / lambda, kappa));
      }
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
      ret = Calloc(1, double);
      ret[0] = n * (gaussian_const + 0.5 * log(high_prec));
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      // The prior is Gaussian for all parameters that are estimated.
      ret = Calloc(1, double);
      ret[0] = 0;
      if (lambda_prior[0] != 0) {
      	ret[0] += gaussian_const - log(lambda_prior[2]) - pow(log(lambda) - lambda_prior[1], 2) / (2 * pow(lambda_prior[2], 2));
      }
      if (kappa_prior[0] != 0) {
      	ret[0] += gaussian_const - log(kappa_prior[2]) - pow(log(kappa) - kappa_prior[1], 2) / (2 * pow(kappa_prior[2], 2));
      }
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return ret;
}
