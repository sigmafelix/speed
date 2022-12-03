\name{speedmat_old}
\alias{speedmat_old}
\title{Spatially Enhanced and Entropy-Derived Contiguity Matrix (SpEED)}
\usage{
speedmat_old(sf,
         mode = "CDE",
         input_vars = c(),
         bandwidth = NULL,
         q_jsd = 0.05,
         kneigh = NULL,
         queen = TRUE,
         cutoff_dist = NULL,
         cutoff_weightd = 0.001,
         cutoff_weightc = 0.001,
         sup_factor = NULL)
         }
\arguments{
\item{sf}{a \code{sf} object}
\item{mode}{Mode of generation}
\item{input_vars}{input variables to compute Jensen-Shannon divergence}
\item{bandwidth}{bandwidth distance for computing distance-based weights}
\item{q_jsd}{Quantile (centile) of Jensen-Shannon divergence. Less than `q_jsd` will be left for the subsequent analysis.}
\item{kneigh}{the number of nearest points when the k-nearest neighbor method is applied for point dataset}
\item{queen}{Queen's contiguity when polygon datasets will be used}
\item{cutoff_distance}{same as \code{bandwidth} otherwise set.}
\item{cutoff_weightd}{threshold for filtering excessively low distance-based weights}
\item{cutoff_weightc}{threshold for filtering weights if they are excessively low after the addition of two weight matrices (applied when \code{mode=="CDE"})}
\item{sup_factor}{suppressing factor for divergence-based weights. 1 if \code{NULL} is retained.}
}
\description{
}
\examples{
library(sf)
library(spatialreg)
library(dplyr)

download.file('https://geodacenter.github.io/data-and-lab/data/natregimes.zip', 'natregimes.zip')
unzip('natregimes.zip', files = 'natregimes.gpkg')
nats <- st_read('natregimes.gpkg')
nats <- st_transform(nats, 2163)

nats <- mutate_if(nats, is.numeric, list(~as.vector(scale(.))))
nats.form <- HR90 ~ RD90 + PS90 + UE90 + DV90 + MA90 + MFIL89
vvars <- as.character(nats.form)
vvars <- c(vvars[2], strsplit(vvars[3], " + ", TRUE)[[1]])
nats_speed <- speedmat(nats,
                    mode = "CDE",
                    input_vars = vvars,
                    bandwidth = 5e4,
                    q_jsd = 0.03,
                    sup_factor = 0.5)
nats_speed_lw <- mat2listw(nats_speed)
nats_plain_lw <- nb2listw(poly2nb(nats))

# may take some time (~10 min)
nats_se <- errorsarlm(formula = nats.form,
                      data = nats,
                      listw = nats_speed_lw)
nats_pl <- errorsarlm(formula = nats.form,
                      data = nats,
                      listw = nats_plain_lw)
# compare AICs
AIC(nats_se) # 2653.99
AIC(nats_pl) # 6848.12
}
