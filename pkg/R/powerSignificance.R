powerSignificance <- function(zo,
                              c = 1, 
                              level = 0.025,
                              designPrior = "conditional",
                              alternative = "one.sided",
                              d = 0,
                              shrinkage = 0){
    # sanity checks
    if (!(designPrior %in% c("conditional", "predictive", "EB")))
        stop('designPrior must be either "conditional", "predictive", or "EB"')
    if (min(c, na.rm = TRUE) < 0)
        stop("c must be larger than 0")
    if (min(d, na.rm = TRUE) < 0)
        stop("d cannot be negative")
    if ((min(shrinkage, na.rm = TRUE) < 0 || max(shrinkage, na.rm = TRUE) > 1)) 
        stop("shrinkage must be in [0, 1]")
    
    # determine direction of alternative and critical value of zr
    v <- p2z(p = level, alternative = alternative) 
    lowertail <- ifelse(alternative == "less", TRUE, FALSE)
    if (alternative %in% c("one.sided", "two.sided")) zo  <- abs(zo)
    
    # shrinkage is the shrinkage factor; s is 1-shrinkage factor
    s <- 1 - shrinkage
    
    # determine parameters of predictive distribution of tr
    if(designPrior == "conditional"){
        mu <- s*zo*sqrt(c)
        sigma <- 1
    }
    if(designPrior == "predictive"){
        mu <- s*zo*sqrt(c)
        sigma <- sqrt(c + 1 + 2*d*c)
    }
    if (designPrior == "EB"){
        s <- pmax(1 - (1 + d)/zo^2, 0)
        mu <- s*zo*sqrt(c)
        sigma <- sqrt(s*c*(1 + d) + 1 + d*c)
    }
    
    # compute replication probability
    pSig <- pnorm(q = v, mean = mu, sd = sigma, lower.tail = lowertail)
    # if (alternative == "two.sided") pSig + pnorm(q = -v, mean = mu, sd = sigma)
    return(pSig)
}
