sampleSizeReplicationSuccess <- function(zo,
                                         power, 
                                         level = thresholdSceptical(level = 0.025, 
                                                                    alternative = alternative, 
                                                                    type = "golden"),
                                         designPrior = "conditional",
                                         alternative = "one.sided"){
    
    # target function for uniroot
    target <- function(zo, c, p, level, designPrior, alternative){
        zr2 <- zr2.quantile(zo = zo, c = c, p = p, designPrior = designPrior)
        pC <- pSceptical(zo = zo, zr = sqrt(zr2), c = c, alternative = alternative)
        return(pC - level)
    }
    mylower <- 0
    myupper <- 1000
    
    # vectorize function in all arguments
    cV <- mapply(FUN = function(zo, power, level, designPrior, alternative) {
        # sanity checks
        if (!(designPrior %in% c("conditional", "predictive")))
            stop('designPrior must be either "conditional" or "predictive"')
        if(!is.numeric(power) || (power <= 0 || power >= 1))
            stop("power must be numeric and in (0,1)!")
        if(!is.numeric(level) || (level <= 0 || level >= 1))
            stop("level must be numeric and in (0,1)!")
        
        target.l <- target(c = mylower, 
                           zo = zo, 
                           p = 1 - power, 
                           level = level,
                           designPrior = designPrior,
                           alternative = alternative)
        target.u <- target(c = myupper, 
                           zo = zo, 
                           p = 1 - power, 
                           level = level,
                           designPrior = designPrior,
                           alternative = alternative)
        if (sign(target.l) == sign(target.u)) {
            if(sign(target.u) > 0)
                c <- Inf
            else 
                c <- NA
        }
        else {
            c <- uniroot(f = target, 
                         lower = mylower, 
                         upper = myupper, 
                         zo = zo, 
                         p = 1 - power, 
                         level = level,
                         designPrior = designPrior,
                         alternative = alternative)$root
        }
        return(c)
    }, zo, power, level, designPrior, alternative)
    
    return(cV)
}
    
