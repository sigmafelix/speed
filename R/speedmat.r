### SpEED function code
#' Legacy SpEED implementation
#' @description Initial implementation of SpEED. Users can choose contiguity-based, distance-based, and mixed spatial relation to get SpEED.
#' @param mode Computation mode.
#' @param input_vars character. Vector of column names to use.
#' @param cutoff_dist numeric(1). Cutoff distance. Should be positive.
#' @param cutoff_weightd numeric(1). Cutoff value of distance weight when \code{mode \%in\% c("CDE", "DE")} Below this value will be truncated to zero.
#' @param cutoff_weightc numeric(1). Cutoff value of distance weight when \code{mode \%in\% c("CDE", "CE")} Below this value will be truncated to zero.
#' @param bandwidth numeric(1). Bandwidth distance to cutoff nearby features beyond the given value. Polygon inputs will be represented by their centroids.
#' @param q_jsd numeric(1). Quantile of Jensen-Shannon divergence to apply cutoffs. Should be bounded [0, 1]
#' @param kneigh integer(1). Number of the nearest neighbors to select nearby features.
#' @param queen logical(1). Queen's (sharing a vertex) or Rook's (sharing a line segment) contiguity for polygon spatial data.
#' @param sup_factor numeric(1). Custom weight to compromise distance-based weights and Jensen-Shannon divergence. Should be bounded [0, 1].
#' @note This function is no longer in maintenance as of January 2024.
#' @author Insang Song (\email{geoissong@gmail.com})
#' @importFrom sf st_relate
#' @importFrom sf st_set_geometry
#' @importFrom sf st_is
#' @importFrom sf st_distance
#' @importFrom sf st_centroid
#' @importFrom nngeo st_nn
#' @importFrom philentropy JSD
#' @importFrom dplyr select
#' @importFrom dplyr mutate_at
#' @importFrom dplyr everything
#' @importFrom dplyr vars
#' @export
speedmat_legacy <-
  function(sf,
    mode = c("CDE", "CE", "DE"),
    input_vars = c(), # Contiguity Distance Entropy
    cutoff_dist = NULL,
    cutoff_weightd = 0.001,
    cutoff_weightc = 0.001,
    bandwidth = NULL, # bandwidth: Currently only supports gaussian
    q_jsd = 0.05,     # quantile of j-s divergence
    kneigh = NULL,    # k-nearest neighbor for point data
    queen = TRUE,     # Queen's contiguity for "C" mode for polygon data
    sup_factor = 0.5
  ) {
    mode <- match.arg(mode)
    if (is.null(bandwidth) && grepl("D", mode))
      stop("No bandwidth was entered")
    if (length(input_vars) == 0)
      stop("No input variable were specified")
    if (!is.null(bandwidth)) {
      cutoff_dist <- bandwidth
    }
    sfp <- sf[1, ]

    sf_ent <- sf |>
      dplyr::select(input_vars) |>
      sf::st_set_geometry(NULL) |>
      dplyr::mutate_at(
        .vars = dplyr::vars(input_vars),
        .funs = list(~as.vector(scale(.)))
      ) |>
      dplyr::mutate_at(
        .vars = dplyr::vars(dplyr::everything()),
        .funs = list(~. + abs(min(.)))
      ) |>
      as.matrix() |>
      philentropy::JSD()
    sf_ent <- sf_ent * (sf_ent <= quantile(sf_ent, q_jsd))
    sf_ent_ex <- exp(-1 * sf_ent)
    sf_ent_f <- (sf_ent_ex * (sf_ent_ex != 1))

    if (grepl("C", mode)){
      if (queen){
        pat <- "F***T****"
      } else {
        pat <- "F***1****"
      }

      if (any(sf::st_is(sfp, "POLYGON"),
              sf::st_is(sfp, "MULTIPOLYGON"),
              sf::st_is(sfp, "POLYGON Z"),
              sf::st_is(sfp, "MULTIPOLYGON Z"))) {
        sf_touch <-
          sf::st_relate(sf, pattern = pat) |>
          # poly2nb(sf, queen = queen) %>%
          as("matrix")
        sf_touch_ <- sf_touch
      } else {
        sf_touch <-
          nngeo::st_nn(
            sf, sf,
            k = kneigh + 1L,
            sparse = FALSE,
            returnDist = TRUE
          )
        sf_touch_ <- sf_touch
        sf_touch_nn <- sf_touch$nn
        diag(sf_touch_nn) <- FALSE
        sf_touch <- sf_touch_nn
      }

      # 122820
      sf_invdv <-
        ((1 - sup_factor) * sf_touch) +
        (sup_factor * sf_ent_f)
    }
    if (grepl("D", mode)) {
      if (any(sf::st_is(sfp, "POLYGON"),
              sf::st_is(sfp, "MULTIPOLYGON"),
              sf::st_is(sfp, "POLYGON Z"),
              sf::st_is(sfp, "MULTIPOLYGON Z"))) {
        sf_interd <- as(sf::st_distance(sf::st_centroid(sf)), "matrix")
      } else {
        sf_interd <- as(sf::st_distance(sf), "matrix")
      }

      sf_interd <- sf_interd * (sf_interd <= cutoff_dist)
      sf_invd <- exp(-1 * ((sf_interd ^ 2) / (bandwidth ^ 2)))
      sf_invd <- sf_invd * (sf_invd != 1)
      diag(sf_invd) <- 0

      if (grepl("^C", mode)) {
        if (
          !any(
            sf::st_is(sfp, "POLYGON"),
            sf::st_is(sfp, "MULTIPOLYGON"),
            sf::st_is(sfp, "POLYGON Z"),
            sf::st_is(sfp, "MULTIPOLYGON Z")
          )
        ) {
          sf_touch_ <- sf_touch
          sf_touch_d <- do.call(c, sf_touch_$dist)
          sf_touch_d <- sf_touch_d[sf_touch_d != 0]
          sf_touch_nn[sf_touch_nn == TRUE] <- sf_touch_d
          sf_touch <- sf_touch_nn
        }

        sf_invdv <-
          ((1 - sup_factor) * (sf_touch * sf_invd)) + (sup_factor * sf_ent_f)
        sf_invdv <-
          sf_invdv * (sf_invdv >= cutoff_weightc)
      } else {
        sf_invdv <-
          (
            (1 - sup_factor) *
            (sf_invd * (sf_invd >= cutoff_weightd))
          ) +
          (sup_factor * sf_ent_f)
        sf_invdv <- sf_invdv * (sf_invdv >= cutoff_weightc)
      }
    }

    # error check
    checksum <- apply(sf_invdv, 1, sum)
    if (any(is.nan(checksum) | checksum == 0)) {
      addr_check <- (is.nan(checksum) | checksum == 0)
      diag(sf_invdv)[addr_check] <- 1
    }
    speedmat <- sf_invdv / apply(sf_invdv, 1, sum)
    return(speedmat)
  }



scale_minmax <- function(x) {(x - min(x)) / (max(x) - min(x))}


#' @title SpEED matrix for matching analysis
#' @description Returns spatially enhanced and entropy-derived matrix (SpEED)
#'  for matching analysis
#' @param data sf object. Should include outcome, treatment, and coordinates (especially if \code{speed_mode == "coord"})
#' @param formula formula. in \code{y ~ x} form.
#' @param outcome character(1). Outcome variable name.
#' Default is \code{'outcome'}
#' @param treatment character(1). Treatment variable name.
#' Default is \code{'treatment'}
#' @param mode_speed character(1). One of \code{'product'} \code{(JSD\%*\%D)} or
#'  \code{'coord'} \code{(JSD(cbind(covariates, X, Y)))}
#' @param coords character(2). Names of the columns with x- and y-dimension
#'  coordinates. Only valid for \code{mode_speed=='coord'}
#' @param coords_factor numeric(1).
#'  Coordinate weights after standardization.
#'  Only valid when \code{mode_speed=='coord'}
#' @param caliper_s numeric(1). Caliper for distance
#' @param caliper_jsd numeric(1). Caliper for divergence value
#' @param scale logical(1). Min-max scaling for both distance and
#'  divergence, resulting in range [0,1]
#' @note Please note that SpEED is currently accepting binary treatment values only.
#' @author Insang Song (geoissong@gmail.com)
#' @export
#' @useDynLib speedmat
speedmat <- function(data,
                     formula,
                     outcome = "outcome",
                     treatment = "treatment",
                     mode_speed = c("product", "product2", "coord"),
                     coords = NULL,
                     coords_factor = NULL,
                     caliper_s = NULL,
                     caliper_jsd = NULL,
                     scale = FALSE) {
  mode_speed <- match.arg(mode_speed)

  if ((mode_speed == "coord")) {
    if (is.null(coords)) {
      stop("No input in the argument coords.
      Check if you set the correct value for mode_speed\n")
    }
    if (is.null(coords_factor)) {
      stop("No input in the argument coords_factor.
      Check if you set the correct value for mode_speed\n")
    }
  }

  speedmat_res <-
    switch(mode_speed,
           product = speedmat.jsdist(data,
                                     formula,
                                     outcome,
                                     treatment,
                                     caliper_s = caliper_s,
                                     caliper_jsd = caliper_jsd,
                                     scale = scale),
           product2 = speedmat.jsdist.m(data,
                                        formula,
                                        outcome,
                                        treatment,
                                        caliper_s = caliper_s,
                                        caliper_jsd = caliper_jsd,
                                        scale = scale),
           coord = speedmat.coord(data,
                                  formula,
                                  outcome,
                                  treatment,
                                  coords,
                                  coords_factor))
  return(speedmat_res)
}


#' @title SpEED matrix by accounting for the weighted coordinate variables
#' @description Returns the Jensen-Shannon divergence with weighted coordinates
#' @param data sf object. Should include outcome, treatment, and coordinates 
#' @param formula formula. in \code{y ~ x} form.
#' @param outcome character(1). Outcome variable name.
#' Default is \code{'outcome'}
#' @param treatment character(1). Treatment variable name.
#' Default is \code{'treatment'}
#' @param coords character(2). Names of the columns with x- and y-dimension
#' coordinates. Default is c("X", "Y")
#' @param coords_factor numeric(1). Coordinate weights after standardization.
#' @author Insang Song (geoissong@gmail.com)
#' @export
speedmat.coord <-
  function(data,
           formula,
           outcome = "outcome",
           treatment = "treatment",
           coords = c("X", "Y"),
           coords_factor = 10L) {

    mat_mm <- model.matrix(formula, data)[, -1]

    formula_ext <- formula
    mat_mf <- model.frame(formula_ext, data)
    # 092422
    mat_mf <- mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0)]

    print(dim(mat_mf))
    mat_mmsc <- mat_mm |>
      apply(2, function(x) scale_minmax(x) + 0.001)

    mat_mmsc[, coords] <- mat_mmsc[, coords] * coords_factor
    mat_jsd <- distJSD(t(mat_mmsc))
    mat_jsdt <- t(mat_jsd)
    mat_jsd <- mat_jsd + mat_jsdt
    gc()
    return(mat_jsd) #mat_jsd instead?

  }

#' @title SpEED matrix by multiplying Jensen-Shannon divergence and
#'  geodesic distance (in km)
#' @description Returns a SpEED matrix, which is the product of
#'  Jensen-Shannon divergence and geographic distance matrix
#' @param data sf object. Should include outcome, treatment, and
#'  coordinates
#' @param formula formula. in y ~ x form.
#' @param outcome character(1). Outcome variable name. default is 'outcome'
#' @param treatment character(1). Treatment variable name.
#' Default is 'treatment'
#' @param caliper_s numeric(1). Caliper for distance
#' @param caliper_jsd numeric(1). Caliper for divergence value
#' @param scale logical(1). Min-max scaling for both distance and
#'  divergence, resulting in range [0,1]
#' @param zero_adjust numeric(1). The minimum after min-max scaling.
#' This is for avoiding zeros in divergence computation.
#' @importFrom sf st_distance
#' @importFrom units drop_units
#' @author Insang Song (geoissong@gmail.com)
#' @export
speedmat.jsdist <- function(data,
                            formula,
                            outcome = "outcome",
                            treatment = "treatment",
                            caliper_s = NULL,
                            caliper_jsd = NULL,
                            scale = FALSE,
                            zero_adjust = 1) {

  mat_mm <- model.matrix(formula, data)[, -1]

  formula_ext <- formula
  mat_mf <- model.frame(formula_ext, data)
  # 092422
  mat_mf <- mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0)]

  mat_mf_dim <- dim(mat_mf)
  cat(sprintf("Valid data dimension: [%d, %d]\n",
              mat_mf_dim[1], mat_mf_dim[2]))
  mat_mmsc <-
    apply(mat_mm, 2,
          function(x) scale_minmax(x) + zero_adjust)

  mat_jsd <- distJSD(mat_mmsc)

  mat_geodist <- sf::st_distance(data) # meters
  mat_geodist <- units::drop_units(mat_geodist)
  if (scale) {
    mat_geodist_sc <- scale_minmax(mat_geodist) + zero_adjust
    mat_jsd <- scale_minmax(mat_jsd) + zero_adjust
  } else {
    mat_geodist_sc <- mat_geodist
  }
  if (!is.null(caliper_s)) {
    mat_geodist_sc[which(mat_geodist > caliper_s)] <- Inf
    # 1/min(mat_geodist)
  }
  if (!is.null(caliper_jsd)) {
    mat_jsd[which(mat_jsd - 0.001 > caliper_jsd)] <- Inf
    # 1/min(mat_jsd)
  }

  # Hadamard product
  mat_jsdd <- mat_jsd * mat_geodist_sc
  return(mat_jsdd) #mat_jsd instead?

}


#' @title SpEED matrix by multiplying Jensen-Shannon divergence and
#'  geodesic distance (in km) -- efficient version
#' @description Returns a SpEED matrix, which is the product of
#'  Jensen-Shannon divergence and geographic distance matrix
#' @param data sf object. Should include outcome, treatment, and
#'  coordinates
#' @param formula formula. in y ~ x form.
#' @param outcome character(1). Outcome variable name. default is 'outcome'
#' @param treatment character(1). Treatment variable name.
#'  Default is 'treatment'
#' @param caliper_s numeric(1). Caliper for distance
#' @param caliper_jsd numeric(1). Caliper for divergence value
#' @param scale logical(1). Min-max scaling for both distance and
#'  divergence, resulting in range [0,1]. Default is \code{FALSE}
#' @param zero_adjust numeric(1). The minimum after min-max scaling.
#' @author Insang Song (geoissong@gmail.com)
#' @export
speedmat.jsdist.m <- function(data,
                              formula,
                              outcome = "outcome",
                              treatment = "treatment",
                              caliper_s = NULL,
                              caliper_jsd = NULL,
                              scale = FALSE,
                              zero_adjust = 1) {

  if (length(unique(data[[treatment]])) < 2) {
    stop("The data appears to have less than two treatment values.
    Please check your data has relevant values in the treatment column
    (suggestion: 0/1).")
  }

  indx_tr <- which(data[[treatment]] > 0)
  indx_co <- which(data[[treatment]] == 0)

  mat_mm <- model.matrix(formula, data)[, -1]

  data_tr <- data[indx_tr, ]
  data_co <- data[indx_co, ]
  # formula_ext = update.formula(formula, as.formula(paste('.~.+", outcome)))
  formula_ext <- formula
  mat_mf <- model.frame(formula_ext, data)
  # 092422
  mat_mf <- mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0)]

  print(dim(mat_mf))
  mat_mmsc <- apply(mat_mm, 2, function(x) scale_minmax(x) + zero_adjust)
  mat_mm_tr <- mat_mmsc[indx_tr, ]
  mat_mm_co <- mat_mmsc[indx_co, ]

  mat_jsd <- distJSD2(mat_mm_tr, mat_mm_co)

  mat_geodist <- sf::st_distance(data_tr, data_co)
  mat_geodist <- units::drop_units(mat_geodist)

  if (scale) {
    mat_geodist_sc <- scale_minmax(mat_geodist) + zero_adjust
    mat_jsd <- scale_minmax(mat_jsd) + zero_adjust
  } else {
    mat_geodist_sc <- mat_geodist
  }
  if (!is.null(caliper_s)) {
    mat_geodist_sc[which(mat_geodist > caliper_s)] <- Inf
    # 1/min(mat_geodist)
  }
  if (!is.null(caliper_jsd)) {
    mat_jsd[which(mat_jsd - zero_adjust > caliper_jsd)] <- Inf
    # 1/min(mat_jsd)
  }
  # Hadamard product
  mat_jsdd <- mat_jsd * mat_geodist_sc
  return(mat_jsdd)

}


#' @title Plot the distribution of geographic distance and
#'  Jensen-Shannon divergence
#' @description Shows a multi-panel plot which displays histograms of
#'  geographic distance and JSD
#' @param mode One of "infunc" (internal use for other speed functions) or "independent" (other)
#' @param data sf object. Should include outcome, treatment, and
#'  coordinates
#' @param formula formula. in y ~ x form.
#' @param outcome character(1). Outcome variable name. Default is 'outcome'
#' @param treatment character(1). Treatment variable name.
#' Default is 'treatment'
#' @param zero_adjust numeric(1). A small value that supplements zero
#' for logarithm. Default is 0.001.
#' @author Insang Song (geoissong@gmail.com)
#' @importFrom sf st_distance
#' @importFrom units drop_units
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_histogram
#' @export
plot_dists <- function(mode = c("independent", "infunc"),
                       data = NULL,
                       formula = NULL,
                       outcome = NULL,
                       treatment = "treatment",
                       zero_adjust = 0.001) {
  if (length(unique(data[[treatment]])) < 2) {
    stop("The data appears to have less than two treatment values.
    Please check your data has relevant values in the treatment column
    (suggestion: 0/1).")
  }

  mode <- match.arg(mode)

  if (mode == "independent") {
    indx_tr <- which(data[[treatment]] > 0)
    indx_co <- which(data[[treatment]] == 0)

    mat_mm <- model.matrix(formula, data)[, -1]

    data_tr <- data[indx_tr, ]
    data_co <- data[indx_co, ]

    formula_ext <- formula
    mat_mf <- model.frame(formula_ext, data)
    # 092422
    mat_mf <- mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0)]

    cat(
      sprintf("Data has the valid dimensions of, [%d, %d]\n",
              dim(mat_mf)[1], dim(mat_mf)[2])
    )

    mat_mmsc <- apply(mat_mm, 2, function(x) scale_minmax(x) + zero_adjust)
    mat_mm_tr <- mat_mmsc[indx_tr, ]
    mat_mm_co <- mat_mmsc[indx_co, ]

    vec_jsd <- as.vector(distJSD2(mat_mm_tr, mat_mm_co))

    mat_geodist <- sf::st_distance(data_tr, data_co) 
    vec_geodist <- as.vector(units::drop_units(mat_geodist))
  } else {
    mat_geodist <- sf::st_distance(data_tr, data_co)
    vec_geodist <- as.vector(units::drop_units(mat_geodist))
  }

  data_dist <- data.frame(
    rbind(
      cbind(disttype = "Geographic distance", value = vec_geodist),
      cbind(disttype = "Jensen-Shannon divergence", value = vec_jsd)
    )
  )

  ggplot2::ggplot(
    data = data_dist,
    mapping = aes(x = value, color = disttype)
  ) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::facet_wrap(~disttype, scales = "free")
}
