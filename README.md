[![codecov](https://codecov.io/gh/sigmafelix/speed/graph/badge.svg?token=EARXDXZPTY)](https://codecov.io/gh/sigmafelix/speed)

[![cov](https://sigmafelix.github.io/speed/badges/coverage.svg)](https://github.com/sigmafelix/speed/actions)


# `speed`
Spatially-Enhanced and Entropy-Derived Contiguity Matrix (SpEED)

# Installation
- Requirement: C++ compilers, Rcpp, RcppArmadillo, sf
- To install this package, please use the command below. It assumes that users have installed `remotes` package
```{r}
remotes::install_github('sigmafelix/speed')
```

# Concept
- Calculating a Hadamard product of a Jensen-Shannon divergence matrix and a geographic distance matrix to control geographic confounding
- Emulates exact matching, with mixed variables
- Actual matching procedure is done with [https://github.com/kouskeimai/MatchIt](`MatchIt`)
- Inspired by [https://github.com/gpapadog/DAPSm](DAPSm)

# Future works
- TODO: Rcpp level parallelization to reduce execution time
- TODO: connecting/porting to Julia
