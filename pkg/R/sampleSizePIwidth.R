sampleSizePIwidth <- function(zo = NA,
                              w,
                              conf.level = 0.95,
                              designPrior = "predictive",
                              d = 0) { 
    
    # vectorize function in all arguments
    cV <- mapply(FUN = function(zo, w, conf.level, designPrior, d) {
        # sanity checks
        if (!(designPrior %in% c("conditional", "predictive", "EB")))
            stop('designPrior must be either "conditional", "predictive", or "EB"')
        if (!is.numeric(d) || d < 0)
            stop("d must be numeric and cannot be negative")
        if (!is.numeric(conf.level) || (conf.level < 0 || conf.level > 1)) 
            stop("conf.level must be numeric and in [0, 1]")
        if (!is.numeric(w) || w <= 0) 
            stop("w must be numeric and positive")
        if (designPrior == "EB") {
            if (is.na(zo) == TRUE) 
                stop("For EB designPrior 'zo' needs to be specified")
        }
        
        # compute relative sample size for specified relative width
        z <- qnorm(p = (1 + conf.level)/2)
        if (designPrior == "conditional") {
            c <- 1/w^2
        } else {
            if (designPrior == "predictive") s <- 1
            if (designPrior == "EB") s <- pmax(1 - (1 + d)/zo^2, 0)
            c <- 1/(w^2 - s*(1 + d) - d)
        }
        
        # c negative: impossible for PI at this conf.level to have desired width
        if (c <= 0) c <- NA
        
        return(c)
        
        }, zo, w, conf.level, designPrior, d)
    
    return(cV)
}
