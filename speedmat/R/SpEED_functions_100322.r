### SpEED function code

if (!require(pacman)) { install.packages('pacman'); library(pacman) }
p_load(Rcpp, MatchIt)

### JSD-Matching part ####
library(Rcpp)

cppjsd = "
#include <Rcpp.h>

// [[Rcpp::export]]
double cppJSD(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y) {
  return (0.5 * (Rcpp::sum(x * Rcpp::log(x/((x+y)/2))) +
                          Rcpp::sum(y * Rcpp::log(y/((x+y)/2)))));
}

// [[Rcpp::export('distJSD')]]
Rcpp::NumericMatrix foo(const Rcpp::NumericMatrix& inMatrix) {
  size_t cols = inMatrix.ncol();
  Rcpp::NumericMatrix result(cols, cols);
  for (size_t i = 0; i < cols; i++) {
    for (size_t j = 0; j < i; j++) {
      result(i,j) = cppJSD(inMatrix(Rcpp::_, i), inMatrix(Rcpp::_, j));
    }
  }
  return result;
}
"
Rcpp::sourceCpp(code = cppjsd)





#' @description Returns Jensen-Shannon divergence and geographic distance matrix with an inflating factor
#' @param data data.frame
#' @param formula formula in y~x
#' @param outcome character(1) default value is "outcome."
#' @param treatment_flag character(1) designates the treat column name
#' @param coords character(2) designates the names of the columns with x- and y-dimension coordinates, respectively.
#' @param coords_factor numeric(1) how standardized coordinates will be exaggerated with respect to the factor. default value is 10. 
#' @author Insang Song
#' 
#'
matchingopt_jsd = function(data, formula, outcome = 'outcome', treatment_flag = 'treat', coords = c('X', 'Y'), coords_factor = 10) {
    
    scale_minmax = function(x) x / (max(x) - min(x))
    mat_mm = model.matrix(formula, data)[,-1]
    # formula_ext = update.formula(formula, as.formula(paste('.~.+', outcome)))
    formula_ext = formula
    mat_mf = model.frame(formula_ext, data)
    # 092422
    mat_mf = mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0 )]

    print(dim(mat_mf))
    mat_mmsc = mat_mm %>% 
        apply(2, function(x) abs(scale_minmax(x)) + scale_minmax(x) + 0.001)
        # apply(2, function(x)  as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001)
        #apply(2, function(x) if (length(unique(x)) > 2) as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001 else x)
    mat_mmsc[,coords] = mat_mmsc[,coords] * coords_factor
    mat_jsd = distJSD(t(mat_mmsc))#philentropy::JSD(mat_mmsc)
    mat_jsdt= t(mat_jsd)
    mat_jsd = mat_jsd + mat_jsdt
    gc()
    #mat_jsd_pm = pairmatch(mat_jsd)
    return(mat_jsd) #mat_jsd instead?

}

# Method 2: JSD * geodesic distance in km
#' @description Returns a SpEED matrix, which is the product of Jensen-Shannon divergence and geographic distance matrix
#' @param data data.frame
#' @param formula formula in y~x
#' @param outcome character(1) default value is "outcome."
#' @param treatment_flag character(1) designates the treat column name
#' @author Insang Song
#'
matchingopt_jsdd = function(data, formula, outcome = 'outcome', treatment_flag = 'treat') {
    
    scale_minmax = function(x) x / (max(x) - min(x))
    mat_mm = model.matrix(formula, data)[,-1]
    # formula_ext = update.formula(formula, as.formula(paste('.~.+', outcome)))
    formula_ext = formula
    mat_mf = model.frame(formula_ext, data)
    # 092422
    mat_mf = mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0 )]

    print(dim(mat_mf))
    mat_mmsc = mat_mm %>% 
        apply(2, function(x) abs(scale_minmax(x)) + scale_minmax(x) + 0.001)
        # apply(2, function(x)  as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001)
        #apply(2, function(x) if (length(unique(x)) > 2) as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001 else x)
    #  mat_mmsc[,coords] = mat_mmsc[,coords] * coords_factor
    mat_jsd = distJSD(t(mat_mmsc))#philentropy::JSD(mat_mmsc)
    mat_jsdt= t(mat_jsd)
    mat_jsd = mat_jsd + mat_jsdt

    geodist = sf::st_distance(data) / 1e3
    mat_jsdd = mat_jsd * geodist
    gc()
    #mat_jsd_pm = pairmatch(mat_jsd)
    return(mat_jsdd) #mat_jsd instead?

}

