#include <RcppArmadillo.h>
#include <cmath>
//using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double cppJSD(const arma::vec& p, const arma::vec& q) {
  arma::vec m = 0.5 * (p + q);
  double jsd = 0.5 * (arma::accu(p % arma::log2(p / m)) + arma::accu(q % arma::log2(q / m)));
  return jsd;
}
// double cppJSD(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y) {
//   return (0.5 * (Rcpp::sum(x * Rcpp::log(x/((x+y)/2))) +
//                           Rcpp::sum(y * Rcpp::log(y/((x+y)/2)))));
// }

//' @name distJSD
//' @title Jensen-Shannon divergence
//' 
//' @description Computes Jensen-Shannon divergence, which is based on Kullback-Leibler divergence. deq(\text{JSD}(A||B)=\frac{1}{2} \left( \text{KL}(A||M)+\text{KL}(B||M) \right))
//' 
//' @param inMatrix a positive-valued numeric matrix.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix distJSD(const Rcpp::NumericMatrix& inMatrix) {
  size_t rows = inMatrix.nrow();
  Rcpp::NumericMatrix result(rows, rows);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < i; j++) {
      result(i,j) = cppJSD(inMatrix(i, Rcpp::_), inMatrix(j, Rcpp::_));
      result(j,i) = result(i,j);
      if (i == j) {
        result(i,j) = 0;
      }
    }
  }
  return result;
}



//' @name distJSD2
//' @title Jensen-Shannon divergence with two inputs
//' 
//' @description Computes Jensen-Shannon divergence, which is based on Kullback-Leibler divergence. deq(\text{JSD}(A||B)=\frac{1}{2} \left( \text{KL}(A||M)+\text{KL}(B||M) \right))
//' 
//' @param inMatrixTr a positive-valued numeric matrix of treated units.
//' @param inMatrixCo a positive-valued numeric matrix of controlled units.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix distJSD2(const Rcpp::NumericMatrix& inMatrixTr, const Rcpp::NumericMatrix& inMatrixCo) {
  size_t rows_tr = inMatrixTr.nrow();
  size_t rows_co = inMatrixCo.nrow();
  Rcpp::NumericMatrix result(rows_tr, rows_co);
  for (size_t i = 0; i < rows_tr; i++) {
    for (size_t j = 0; j < rows_co; j++) {
      result(i,j) = cppJSD(inMatrixTr(i, Rcpp::_), inMatrixCo(j, Rcpp::_));
    }
  }
  return result;
}
