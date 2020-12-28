## Insang Song (isong@uoregon.edu)
## Spatially Enhanced and Entropy-Derived MATrix (SpEED-MAT)
speedmat <- function(sf,
                     mode = "CDE", # "CE"/"DE"/"CDE"
                     input_vars = c(), # Contiguity Distance Entropy
                     bandwidth = NULL, # bandwidth: Currently only supports gaussian
                     q_jsd = 0.05,     # quantile of j-s divergence
                     kneigh = NULL,    # k-nearest neighbor for point data
                     queen = TRUE,     # Queen's contiguity for "C" mode for polygon data
                     cutoff_dist = NULL,
                     cutoff_weightd = 0.001,
                     cutoff_weightc = 0.001,
                     sup_factor = NULL){
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
    if (!is.null(sup_factor)){
        sf_ent_f <- sf_ent_f * sup_factor
    }

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
                #knn2nb %>%
                #nb2mat
        }
        sf_invdv <- sf_touch + sf_ent_f
    }
    if (grepl("D", mode))
    {
        if (any(st_is(sfp, 'POLYGON'), st_is(sfp, 'MULTIPOLYGON'), st_is(sfp, "POLYGON Z"), st_is(sfp, "MULTIPOLYGON Z"))) {
            sf_interd <- as(st_distance(st_centroid(sf)), 'matrix')
        } else {
            sf_interd <- as(st_distance(sf), 'matrix')
        }
        #if (!is.null(cutoff_dist)){
        sf_interd <- sf_interd * (sf_interd <= cutoff_dist)
        #}
        sf_invd <- exp(-1 * ((sf_interd^2)/(bandwidth^2)))
        sf_invd <- sf_invd * (sf_invd != 1)
        #
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

            sf_invdv <- (sf_touch * sf_invd) + sf_ent_f
            sf_invdv <- sf_invdv * (sf_invdv >= cutoff_weightc)
        } else {
            sf_invdv <- (sf_invd * (sf_invd >= cutoff_weightd)) + sf_ent_f
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
