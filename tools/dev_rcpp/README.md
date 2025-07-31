## Notes
- If you want to roll back to `Rcpp` version:

- First up, Rcpp should be linked in the `DESCRIPTION` file:

```
LinkingTo: Rcpp, RcppArmadillo
```

- Copy three files to the matching directories:
  - `RcppExports.R` to `R/`
  - `RcppExports.cpp` to `src/`
  - `distJSD.cpp` to `src/`