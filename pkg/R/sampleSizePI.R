sampleSizePI <- function(po = NULL,
                         to = p2t(po, alternative = alternative),
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
    
    # compute relative sample size for PI at this conf.level not to include 0
    z <- qnorm(p = (1 + conf.level)/2)
    if (designPrior == "conditional") {
        c <- z^2/to^2
    } else {
        if (designPrior == "predictive") s <- 1
        if (designPrior == "EB") s <- pmax(1 - (1 + d)/to^2, 0)
        c <- 1/(s^2*to^2/z^2 - s*(1 + d) - d)
    }
    
    # c negative: impossible for PI at this conf.level not to include 0
    c <- ifelse(c <= 0, NA, c) 
    return(c)
}