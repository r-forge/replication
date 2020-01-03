sampleSizePIwidth <- function(po = NULL,
                              to = p2t(po, alternative = alternative),
                              w,
                              conf.level = 0.95,
                              alternative = "two.sided",
                              designPrior = "predictive",
                              d = 0) { 
    # sanity checks
    if (!(designPrior %in% c("conditional", "predictive", "EB")))
        stop('designPrior must be either "conditional", "predictive", or "EB"')
    if (min(d, na.rm = TRUE) < 0)
        stop("d cannot be negative")
    if (conf.level <= 0 | conf.level >= 1) 
        stop("conf.level must be in (0, 1)")
    if (any(w <= 0)) 
        stop("w must be positive")
    if (designPrior == "EB") {
        if (length(to) == 0) 
            stop("For EB designPrior 'to' needs to be specified")
    }
    
    # compute relative sample size for specified relative width
    z <- qnorm(p = (1 + conf.level)/2)
    if (designPrior == "conditional") {
        c <- 1/w^2
    } else {
        if (designPrior == "predictive") s <- 1
        if (designPrior == "EB") s <- pmax(1 - (1 + d)/to^2, 0)
        c <- 1/(w^2 - s*(1 + d) - d)
    }
    
    # c negative: impossible for PI at this conf.level to have desired width
    c <- ifelse(c <= 0, NA, c) 
    return(c)
}