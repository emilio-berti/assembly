#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute weighted Jaccard similarity
//'
//' @param A numeric matrix
//' @return The matrix with Jaccard similarity
// [[Rcpp::export]]
Rcpp::List wJaccard( arma::mat A ) {
  arma::mat sim;
  sim.copy_size(A);
  arma::vec Min = A.col(0);
  arma::vec Max = A.col(0);
  Min.zeros();
  Max.zeros();
  for (arma::uword i = 0; i < sim.n_rows; i++) {
    for (arma::uword j = 0; j < sim.n_cols; j++) {
      if (i < j) {
        for (arma::uword k = 0; k < sim.n_rows; k++) {
          Min(k) = A(k, i) < A(k, j) ? A(k, i) : A(k, j);
          Max(k) = A(k, i) > A(k, j) ? A(k, i) : A(k, j);
        }
        sim(i, j) = sum(Min) / sum(Max);
        sim(j, i) = sim(i, j);
        Min.zeros();
        Max.zeros();
      } else if (i == j) {
        sim(i, j) = 1;
      }
    }
  }
  return ( Rcpp::List::create(Rcpp::Named("Jaccard") = sim) );
}


