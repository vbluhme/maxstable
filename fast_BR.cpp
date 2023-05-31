#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

mat rPx_BR(int n, mat& x, mat& L, mat& sig_n) {
  int nrow = x.n_rows;
  vec Z(nrow);
  std::generate( Z.begin(), Z.end(), ::norm_rand );
  vec W = L.t() * Z;
  
  vec Y = exp(W - W(n) - sig_n / 2);

  return Y;
}
  
int all_smaller(vec a, vec b) {
  int n = a.n_elem;
  int all_smaller = 1;
  for(int i = 0; i < n; i++) {
    if(a(i) >= b(i)) {all_smaller = 0;}
  }
  return all_smaller;
}
vec pmax(vec& a, vec& b) {
  int n = a.n_elem;
  vec out(n);
  for(int i = 0; i < n; i++) {
    out(i) = std::max(a(i), b(i));
  }
  return out;
}

// [[Rcpp::export]]
vec fast_BR(mat x, double alpha, double s, mat L) {
  int nx = x.n_rows;
  
  vec sig_n(nx);
  for(int i = 0; i < nx; i++) {
    sig_n(i) = pow(pow(x(i,0) - x(0, 0), 2) + pow(x(i,1) - x(0, 1), 2), alpha / 2) / pow(s, alpha);
  }
 
  double Gamma = rexp(1)(0);
  vec zY = rPx_BR(0, x, L, sig_n) / Gamma;

  vec Z = zY;

  if (nx == 1) return Z;

  for(int n = 1; n < nx; n++) {
    //Rcout << "\rn: " << n+1 << " / " << nx << " (rcpp)";
    for(int i = 0; i < nx; i++) {
      sig_n(i) = pow(pow(x(i,0) - x(n, 0), 2) + pow(x(i,1) - x(n, 1), 2), alpha / 2) / pow(s, alpha);
    }
    
    Gamma = rexp(1)(0);
    while(1/Gamma > Z(n)) {
      zY = rPx_BR(n, x, L, sig_n) / Gamma;

      if(all_smaller(zY.subvec(0, n-1), Z.subvec(0, n-1))) {
        Z = pmax(Z, zY);
      }
      
      Gamma = Gamma + rexp(1)(0);
    }
  }
  //Rcout << "\n";

  return Z;
}