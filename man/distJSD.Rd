% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{distJSD}
\alias{distJSD}
\title{Jensen-Shannon divergence}
\usage{
distJSD(inMatrix)
}
\arguments{
\item{inMatrix}{a positive-valued numeric matrix.}
}
\description{
Computes Jensen-Shannon divergence, which is based on Kullback-Leibler divergence. deq(\text{JSD}(A||B)=\frac{1}{2} \left( \text{KL}(A||M)+\text{KL}(B||M) \right))
}
