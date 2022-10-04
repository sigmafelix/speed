# `speed`
Spatially-Enhanced and Entropy-Derived Contiguity Matrix (SpEED)

```{r}
remote::install_github('sigmafelix/speed/speedmat')
```

# Concept
- Calculating a Hadamard product of a Jensen-Shannon divergence matrix and a geographic distance matrix to control geographic confounding
- Emulates exact matching, with mixed variables
- Actual matching procedure is done with [https://github.com/kouskeimai/MatchIt](`MatchIt`)
- Inspired by [https://github.com/gpapadog/DAPSm](DAPSm)

# Future works
- Improve the speed for scalable analyses by making the pairwise SpEED matrix more succinct (i.e., all records to treat-control alignments)
- Enable calipers at the calculation steps for divergence, distance, and the product of both, respectively
