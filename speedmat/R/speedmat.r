### SpEED function code
#' @export
speedmat_old <- function(sf,
                     mode = "CDE", # "CE"/"DE"/"CDE"
                     input_vars = c(), # Contiguity Distance Entropy
                     bandwidth = NULL, # bandwidth: Currently only supports gaussian
                     q_jsd = 0.05,     # quantile of j-s divergence
                     kneigh = NULL,    # k-nearest neighbor for point data
                     queen = TRUE,     # Queen's contiguity for "C" mode for polygon data
                     cutoff_dist = NULL,
                     cutoff_weightd = 0.001,
                     cutoff_weightc = 0.001,
                     sup_factor = 0.5){
    if (is.null(bandwidth) & grepl("D", mode))
        stop("No bandwidth was entered")
    if (length(input_vars) == 0)
        stop("No input variable were specified")
    if (!grepl("E", mode))
        stop("Please use the methods for standard spatial weight matrices for non-entropy matrices")
    if (!is.null(bandwidth)) {
        cutoff_dist <- bandwidth
    }
    sfp <- sf[1,]

    sf_ent <- sf %>%
        dplyr::select(input_vars) %>%
        st_set_geometry(NULL) %>%
        mutate_at(.vars = vars(input_vars),
                  .funs = list(~as.vector(scale(.)))) %>%
        mutate_at(.vars = vars(everything()),
                  .funs = list(~.+abs(min(.)))) %>%
        as.matrix %>%
        philentropy::JSD(.)
    sf_ent <- sf_ent * (sf_ent <= quantile(sf_ent, q_jsd))
    sf_ent_ex <- exp(-1 * sf_ent)
    sf_ent_f <- (sf_ent_ex * (sf_ent_ex != 1))

    if (grepl("C", mode)){
        if (queen){
            pat <- 'F***T****'
        } else {
            pat <- 'F***1****'
        }

        if (any(st_is(sfp, 'POLYGON'), st_is(sfp, 'MULTIPOLYGON'), st_is(sfp, "POLYGON Z"), st_is(sfp, "MULTIPOLYGON Z")))
            {sf_touch <- st_relate(sf, pattern = pat) %>% #poly2nb(sf, queen = queen) %>%
                as(., 'matrix')#nb2mat
             sf_touch_ <- sf_touch
            } else {
            sf_touch <- st_nn(sf, sf, k = kneigh + 1, sparse = FALSE, returnDist = TRUE) #%>%
            sf_touch_ <- sf_touch
            sf_touch_nn <- sf_touch$nn
            diag(sf_touch_nn) <- FALSE
            sf_touch <- sf_touch_nn
            }
        # 122820
        sf_invdv <- ((1-sup_factor)*sf_touch) + (sup_factor * sf_ent_f)
    }
    if (grepl("D", mode))
    {
        if (any(st_is(sfp, 'POLYGON'), st_is(sfp, 'MULTIPOLYGON'), st_is(sfp, "POLYGON Z"), st_is(sfp, "MULTIPOLYGON Z"))) {
            sf_interd <- as(st_distance(st_centroid(sf)), 'matrix')
        } else {
            sf_interd <- as(st_distance(sf), 'matrix')
        }
        sf_interd <- sf_interd * (sf_interd <= cutoff_dist)
        sf_invd <- exp(-1 * ((sf_interd^2)/(bandwidth^2)))
        sf_invd <- sf_invd * (sf_invd != 1)
        diag(sf_invd) <- 0

        if (grepl("^C", mode)) {
            if (!any(st_is(sfp, 'POLYGON'), st_is(sfp, 'MULTIPOLYGON'), st_is(sfp, "POLYGON Z"), st_is(sfp, "MULTIPOLYGON Z")))
            {
                sf_touch_ <- sf_touch
                sf_touch_d <- do.call(c, sf_touch_$dist)
                sf_touch_d <- sf_touch_d[sf_touch_d != 0]
                sf_touch_nn[sf_touch_nn == TRUE] <- sf_touch_d
                sf_touch <- sf_touch_nn
            }

            sf_invdv <- ((1-sup_factor) * (sf_touch * sf_invd)) + (sup_factor * sf_ent_f)
            sf_invdv <- sf_invdv * (sf_invdv >= cutoff_weightc)
        } else {
            sf_invdv <- ((1-sup_factor)* (sf_invd * (sf_invd >= cutoff_weightd))) + (sup_factor * sf_ent_f)
            sf_invdv <- sf_invdv * (sf_invdv >= cutoff_weightc)
        }
    }

    # error check
    checksum <- apply(sf_invdv, 1, sum)
    if (any(is.nan(checksum) | checksum == 0)){
        addr_check <- (is.nan(checksum) | checksum==0)
        diag(sf_invdv)[addr_check] <- 1
    }
    speedmat <- sf_invdv / apply(sf_invdv, 1, sum)
    return(speedmat)
}


#' @title SpEED matrix for matching analysis
#' @description Returns spatially enhanced and entropy-derived matrix (SpEED) for matching analysis
#' @param data sf object. Should include outcome, treatment, and coordinates (if \code{speed_mode == "coord"})
#' @param formula formula. in \code{y ~ x} form.
#' @param outcome character(1). Outcome variable name. default is \code{'outcome'}
#' @param treatment character(1). Treatment variable name. default is \code{'treatment'}
#' @param mode_speed character(1). One of \code{'product'} \code{(JSD\%*\%D)} or \code{'coord'} \code{(JSD(cbind(covariates, X, Y)))}
#' @param coords character(2). Names of the columns with x- and y-dimension coordinates. Only valid for \code{mode_speed=='coord'}
#' @param coords_factor numeric(1). Coordinate weights after standardization. Only valid for \code{mode_speed=='coord'}
#' @param caliper_s numeric(1). Caliper for distance
#' @param caliper_jsd numeric(1). Caliper for divergence value
#' @param scale logical(1). Min-max scaling for both distance and divergence, resulting in range [0,1]
#' @author Insang Song (sigmafelix@hotmail.com)
#' @export
#' @useDynLib speedmat
#' 
speedmat = function(data,
                    formula,
                    outcome = 'outcome',
                    treatment = 'treatment',
                    mode_speed = 'product',
                    coords = NULL,
                    coords_factor = NULL, 
                    caliper_s = NULL, 
                    caliper_jsd = NULL, 
                    scale = FALSE) {
    if (!mode_speed %in% c('product', 'coord')) {
        stop("The argument mode_speed should be one of 'product' or 'coord'")
    }
    if ((mode_speed == 'coord' && (is.null(coords) | is.null(coords_factor)))) {
        stop("No coords and coords_factor arguments. Check if you set the correct value for mode_speed")
    }
    
    speedmat_res = 
        switch(mode_speed,
              product = speedmat.jsdist(data, formula, outcome, treatment, caliper_s = caliper_s, caliper_jsd = caliper_jsd, scale = scale),
              coord = speedmat.coord(data, formula, outcome, treatment, coords, coords_factor))
    return(speedmat_res)
    }


#' @title SpEED matrix by accounting for the weighted coordinate variables
#' @description Returns the Jensen-Shannon divergence with weighted coordinates
#' @param data sf object. Should include outcome, treatment, and coordinates 
#' @param formula formula. in \code{y ~ x} form.
#' @param outcome character(1). Outcome variable name. default is \code{'outcome'}
#' @param treatment character(1). Treatment variable name. default is \code{'treatment'}
#' @param coords character(2). Names of the columns with x- and y-dimension coordinates. 
#' @param coords_factor numeric(1). Coordinate weights after standardization. 
#' @author Insang Song (sigmafelix@hotmail.com)
#' @export
#'
speedmat.coord = function(data, formula, outcome = 'outcome', treatment = 'treatment', coords = c('X', 'Y'), coords_factor = 10) {
    
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

#' @title SpEED matrix by multiplying Jensen-Shannon divergence and geodesic distance (in km)
#' @description Returns a SpEED matrix, which is the product of Jensen-Shannon divergence and geographic distance matrix
#' @param data sf object. Should include outcome, treatment, and coordinates 
#' @param formula formula. in \code{y ~ x} form.
#' @param outcome character(1). Outcome variable name. default is \code{'outcome'}
#' @param treatment character(1). Treatment variable name. default is \code{'treatment'}
#' @param caliper_s numeric(1). Caliper for distance
#' @param caliper_jsd numeric(1). Caliper for divergence value
#' @param scale logical(1). Min-max scaling for both distance and divergence, resulting in range [0,1]
#' @author Insang Song (sigmafelix@hotmail.com)
#' @export
#' 
speedmat.jsdist = function(data, formula, outcome = 'outcome', treatment = 'treatment', caliper_s = NULL, caliper_jsd = NULL, scale = FALSE) {
    
    scale_minmax = function(x) x / (max(x) - min(x))
    mat_mm = model.matrix(formula, data)[,-1]
    # formula_ext = update.formula(formula, as.formula(paste('.~.+', outcome)))
    formula_ext = formula
    mat_mf = model.frame(formula_ext, data)
    # 092422
    mat_mf = mat_mf[, sapply(mat_mf, function(x) length(unique(x)) != 0 )]

    print(dim(mat_mf))
    mat_mmsc = apply(mat_mm, 2, function(x) abs(scale_minmax(x)) + scale_minmax(x) + 0.001)
        # apply(2, function(x)  as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001)
        #apply(2, function(x) if (length(unique(x)) > 2) as.vector(scale(x)) + abs(min(as.vector(scale(x)))) + 0.001 else x)
    #  mat_mmsc[,coords] = mat_mmsc[,coords] * coords_factor
    mat_jsd = distJSD(t(mat_mmsc))#philentropy::JSD(mat_mmsc)
    mat_jsdt= t(mat_jsd)
    mat_jsd = mat_jsd + mat_jsdt

    mat_geodist = sf::st_distance(data) / 1e3
    if (!is.null(caliper_s)) {
        mat_geodist[which(mat_geodist > caliper_s)] = NA

    }
    if (!is.null(caliper_jsd)) {
        mat_jsd[which(mat_jsd > caliper_jsd)] = NA
    }
    if (scale) {
        mat_geodist = (mat_geodist - min(mat_geodist, na.rm = TRUE)) / max(mat_geodist, na.rm = TRUE)
        mat_jsd = (mat_jsd - min(mat_jsd, na.rm = TRUE)) / max(mat_jsd, na.rm = TRUE)
        
    }
    gc()
    #mat_jsd_pm = pairmatch(mat_jsd)
    # Hadamard product
    mat_jsdd = mat_jsd * mat_geodist
    return(mat_jsdd) #mat_jsd instead?

}

