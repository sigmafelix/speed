#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double cppJSD(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y) {
  return (0.5 * (Rcpp::sum(x * Rcpp::log(x/((x+y)/2))) +
                          Rcpp::sum(y * Rcpp::log(y/((x+y)/2)))));
}

//' @name distJSD
//' @title Jensen-Shannon divergence
//' 
//' @description Computes Jensen-Shannon divergence, which is based on Kullback-Leibler divergence. deq(\text{JSD}(A||B)=\frac{1}{2} \left( \text{KL}(A||M)+\text{KL}(B||M) \right))
//' 
//' @param inMatrix a positive-valued numeric matrix.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix distJSD(const Rcpp::NumericMatrix& inMatrix) {
  size_t cols = inMatrix.ncol();
  Rcpp::NumericMatrix result(cols, cols);
  for (size_t i = 0; i < cols; i++) {
    for (size_t j = 0; j < i; j++) {
      result(i,j) = cppJSD(inMatrix(Rcpp::_, i), inMatrix(Rcpp::_, j));
    }
  }
  return result;
}
