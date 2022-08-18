// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cassert>

using namespace Rcpp;

// Some of the functions in this script are heavily inspired by
// https://gallery.rcpp.org/articles/dmvnorm_arma/

static double const log2pi = std::log(2.0 * M_PI);

// C++ version of the dtrmv BLAS function, that computes x * trimat
void inplace_tri_mat_mult(arma::vec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  for(unsigned j = n; j-- > 0;){
    double tmp(0.0);
    for(unsigned i = 0; i <= j; ++i) {
      tmp += trimat.at(i, j) * x[i];
    }
    x[j] = tmp;
  }
}

//' Compute the probability density function of a multivariate Gaussian
//' using the Rcpparmadillo package.
//' Here, x is a (d x n)-dimensional matrix, consisting of n d-dimensional random variables,
//' mean is a d-dimensional vector,
//' sigma is a (d x d)-dimensional covariance matrix,
//' and logd is a bool, stating whether we want the log-density or not.
// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat const &x,  
		       arma::vec const &mean,  
		       arma::mat const &sigma, 
		       bool const logd = false) { 
  arma::uword const n = x.n_cols;
  arma::uword const d = x.n_rows;
  arma::vec out(n);
  arma::mat const inverse_chol = arma::inv(trimatu(arma::chol(sigma)));
  double const lpdf_const =
    arma::sum(log(inverse_chol.diag())) - (double)d * 0.5 * log2pi;
    
  arma::vec z;
  for (arma::uword i = 0; i < n; ++i) {
    z = (x.col(i) - mean);
    inplace_tri_mat_mult(z, inverse_chol);
    out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
  }  

  if (logd) {
    return out;
  } else {
    return exp(out);
  }
}

// Perform fast forward substitution for solving the equation
// Lx = b for x, where L is lower triangular
void forward_substitution_fast(arma::vec &res,
			       arma::mat const &L,
			       arma::vec const &b) {
  arma::uword const n = b.size();
  for (arma::uword i = 0; i < n; ++i) {
    res.at(i) = b.at(i);
    for (arma::uword j = 0; j < i; ++j) {
      res.at(i) -= L.at(i, j) * res.at(j);
    }
    res.at(i) /= L.at(i, i);
  }
}

// Compute dmvnorm quickly when x is not a matrix, but a vector
double dmvnorm_double_fast(arma::vec const &x,  
			   arma::mat const &sigma) {
  arma::uword const d = x.n_rows;
  arma::mat const chol = arma::chol(sigma, "lower");
  assert(chol.n_cols == d);
  assert(chol.n_rows == d);
  double out =
    - arma::sum(log(chol.diag())) - (double)d * 0.5 * log2pi;
    
  arma::vec z(d);
  forward_substitution_fast(z, chol, x);
  out -= 0.5 * arma::dot(z, z);     

  return out;
}

// Compute Σ = A * B * Σ0 * B * A'
arma::mat sigma_func(arma::sp_mat const &A,
		     arma::sp_mat const &B,
		     arma::mat const &sigma,
		     double const nugget) {
  arma::mat tmp = B * sigma * B;
  arma::mat res = A * tmp * A.t();
  res.diag() += nugget;
  return res;
}

//' Compute the (log-)likelihood of the conditional extremes distribution when
//' a and b are given, in the sense that a has already been subtracted in x,
//' while b is assumed to not depend on y0 and given in a diagonal matrix form.
//' 
//' The input variables are:
//' x: an (n x d)-dimensional matrix of observations where a has already been subtracted.
//' A: The projection matrix used for building the SPDE approximation
//' B: A diagonal (m x m)-dimensional matrix containing the values of b at the m
//'   triangular mesh nodes. b is assumed to not depend on y0.
//' sigma0: The covariance matrix of the m Gaussian random variables in the triangular mesh.
//' nugget: The variance of the nugget effect.
//' logd: Boolean explaining if we should return the log-likelihood or the likelihood.
//' na_rm: Boolean explaining how we should treat NA values. If na_rm = false, then
//'   any column of x that returns an NA value will result in an NA value in the output.
//'   If na_rm = true, then we remove the NA values before computing the likelihood.
//'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
//'   the likelihood for a (d-3)-dimensional Gaussian random variable.
// [[Rcpp::export]]
arma::vec dconditional_arma(arma::mat const &x,  
			    arma::sp_mat const &A,
			    arma::sp_mat const &B,
			    arma::mat const &sigma0,
			    double const nugget,
			    bool const logd = true,
			    bool const na_rm = true) { 
  arma::mat sigma = sigma_func(A, B, sigma0, nugget);
  arma::uword const n = x.n_cols;
  arma::uword const d = x.n_rows;
  arma::vec out(n);
  arma::mat const chol = arma::chol(sigma, "lower");
  double const lpdf_const =
    - arma::sum(log(chol.diag())) - (double)d * 0.5 * log2pi;

  arma::uword count;
  arma::uvec good_index(d);
    
  arma::vec z(d);
  for (arma::uword i = 0; i < n; ++i) {
    z = x.col(i);

    if (na_rm) {
      count = 0;
      for (arma::uword j = 0; j < d; ++j) {
	if (!NumericVector::is_na(z(j))) {
	  good_index(count) = j;
	  ++count;
	}
      }
      if (count != d) {
	out(i) = dmvnorm_double_fast(z(good_index.head(count)),
				     sigma(good_index.head(count), good_index.head(count)));
      } else {
	forward_substitution_fast(z, chol, x.col(i));
	out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
      }
    } else {
      forward_substitution_fast(z, chol, x.col(i));
      out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
    }
  }
      
  if (logd) {
    return out;
  } else {
    return exp(out);
  }
}
