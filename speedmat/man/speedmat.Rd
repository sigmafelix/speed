\name{speedmat}
\alias{speedmat}
\title{Spatially Enhanced and Entropy-Derived Contiguity Matrix (SpEED)}
\usage{
speedmat(sf,
         mode = "CDE", # "CE"/"DE"/"CDE"
         input_vars = c(), # Contiguity Distance Entropy
         bandwidth = NULL, # bandwidth: Currently only supports gaussian
         q_jsd = 0.05,     # quantile of j-s divergence
         kneigh = NULL,    # k-nearest neighbor for point data
         queen = TRUE,     # Queen's contiguity for "C" mode for polygon data
         cutoff_dist = NULL,
         cutoff_weightd = 0.001,
         cutoff_weightc = 0.001,
         sup_factor = NULL)
         }
\description{
\item{sf} \description{a `sf` object}
\item{mode} \description{Mode of generation}
}
\examples{
speed()
}
