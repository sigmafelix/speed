% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{distJSD2}
\alias{distJSD2}
\title{Jensen-Shannon divergence with two inputs}
\usage{
distJSD2(inMatrixTr, inMatrixCo)
}
\arguments{
\item{inMatrixTr}{a positive-valued numeric matrix of treated units.}

\item{inMatrixCo}{a positive-valued numeric matrix of controlled units.}
}
\description{
Computes Jensen-Shannon divergence, which is based on Kullback-Leibler divergence. deq(\text{JSD}(A||B)=\frac{1}{2} \left( \text{KL}(A||M)+\text{KL}(B||M) \right))
}
