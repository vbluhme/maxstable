#include <Rcpp.h>
using namespace Rcpp;

// Function for computation of covariance matrix of fBs.

// [[Rcpp::export]]
NumericMatrix fBs_cov(NumericMatrix x, double alpha, double s) {
  int n = x.nrow();
  NumericMatrix cov(n, n);
  double normdiff;
  
  double scale = 1/(2 * pow(s, alpha));
  
  NumericVector normvec(n);
  for(int i = 0; i < n; i++) {
    normvec(i) = pow(pow(x(i,0), 2) + pow(x(i,1), 2), alpha/2);
  }
  
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      normdiff = pow(pow(x(i, 0) - x(j, 0), 2) + pow(x(i, 1) - x(j, 1),2), alpha/2);
      cov(i,j) = scale * (normvec(i) + normvec(j) - normdiff);
      cov(j,i) = cov(i,j);
    }
  }
  return(cov);
}

// [[Rcpp::export]]
NumericVector rPx_BR_sheet_cov(NumericMatrix x, int n, double alpha, double s) {
  int nrow = x.nrow();

  NumericVector out(nrow);
  for(int i = 0; i < nrow; i++) {
    out(i) = pow(pow(x(i,0) - x(n-1, 0), 2) + pow(x(i,1) - x(n-1, 1), 2), alpha / 2) / pow(s, alpha);
  }

  return(out);
}

/*
NumericVector sim_extremalfuns_BR(NumericMatrix x, double alpha, NumericMatrix L) {
  int N = x.nrow();
  NumericVector Z(N), Y(N);
  double Gamma;
  
  Gamma = rexp(1);
  Y = rPx_BR_rcpp(1, x, N);
  
  return(Z);
}
 */