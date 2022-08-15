// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// These functions are modified versions of https://gallery.rcpp.org/articles/dmvnorm_arma/

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::vec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
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
  double const lpdf_const = arma::sum(log(inverse_chol.diag()))
    - (double)d * 0.5 * log2pi;
    
  arma::vec z;
  for (arma::uword i = 0; i < n; i++) {
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

double dconditional_no_na_double(arma::vec z,  
				 arma::mat const &sigma) {
  double const d = z.size();
  arma::mat const inverse_chol = arma::inv(trimatu(arma::chol(sigma)));
  double const lpdf_const = arma::sum(log(inverse_chol.diag()))
    - d * 0.5 * log2pi;

  inplace_tri_mat_mult(z, inverse_chol);
  double out = lpdf_const - 0.5 * arma::dot(z, z);     
      
  return out;
}

arma::mat sigma_func(const arma::sp_mat &A,
		     const arma::sp_mat &B,
		     const arma::mat &sigma) {
  arma::mat tmp = B * sigma * B;
  arma::mat res = A * tmp * A.t();
  return res;
}

// [[Rcpp::export]]
arma::vec dconditional_arma(const arma::mat &x,  
			    const arma::sp_mat &A,
			    const arma::sp_mat &B,
			    const arma::mat &sigma0,
			    const double nugget,
			    const bool logd = true,
			    const bool na_rm = true) { 
  arma::mat sigma = sigma_func(A, B, sigma0);
  arma::uword const n = x.n_cols;
  arma::uword const d = x.n_rows;
  arma::vec out(n);
  arma::mat sigma_with_nugget = sigma;
  sigma_with_nugget.diag() += nugget;
  arma::mat const inverse_chol = arma::inv(trimatu(arma::chol(sigma_with_nugget)));
  double const lpdf_const = arma::sum(log(inverse_chol.diag()))
    - (double)d * 0.5 * log2pi;

  arma::uword count;
  arma::uvec good_index(d);
    
  arma::vec z;
  for (arma::uword i = 0; i < n; i++) {
    z = x.col(i);

    if (na_rm) {
      count = 0;
      for (arma::uword j = 0; j < d; j++) {
	if (!NumericVector::is_na(z(j))) {
	  good_index(count) = j;
	  count += 1;
	}
      }
      if (count != d) {
	out(i) = dconditional_no_na_double(z(good_index.head(count)),
					   sigma_with_nugget(good_index.head(count), good_index.head(count)));
      } else {
	inplace_tri_mat_mult(z, inverse_chol);
	out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
      }
    } else {
      inplace_tri_mat_mult(z, inverse_chol);
      out(i) = lpdf_const - 0.5 * arma::dot(z, z);     
    }
  }
      
  if (logd) {
    return out;
  } else {
    return exp(out);
  }
}


