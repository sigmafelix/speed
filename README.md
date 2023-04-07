# `speedmat`
Spatially-Enhanced and Entropy-Derived Contiguity Matrix (SpEED)

# Installation
- Requirement: C++ compilers, Rcpp, RcppArmadillo, sf
- To install this package, please use the command below. It assumes that users have installed `remotes` package
```{r}
remotes::install_github('sigmafelix/speed/speedmat')
```

# Concept
- Calculating a Hadamard product of a Jensen-Shannon divergence matrix and a geographic distance matrix to control geographic confounding
- Emulates exact matching, with mixed variables
- Actual matching procedure is done with [https://github.com/kouskeimai/MatchIt](`MatchIt`)
- Inspired by [https://github.com/gpapadog/DAPSm](DAPSm)

# Future works
- ~~Improve the speed for scalable analyses by making the pairwise SpEED matrix more succinct (i.e., all records to treat-control alignments)~~ : added in 0.3.0
- ~~Enable calipers at the calculation steps for divergence, distance, and the product of both, respectively~~: added in 0.3.0
- TODO: Rcpp level parallelization to reduce execution time
- TODO: connecting/porting to Julia
