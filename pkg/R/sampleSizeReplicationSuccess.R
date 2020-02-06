sampleSizeReplicationSuccess <- function(zo,
                                         power, 
                                         level = thresholdSceptical(0.025, alternative = "one.sided", type="controlled"),
                                         designPrior = "conditional",
                                         alternative = "one.sided"){
    if(min(power) < 0 || max(power) > 1)
        stop("The power must lie in (0,1)!")
    
    # target function for uniroot
    target <- function(zo, 
                       c, 
                       p, 
                       level, 
                       designPrior, 
                       alternative){
        zr2 <- zr2.quantile(zo = zo, c = c, p = p, designPrior = designPrior)
        pC <- pSceptical(zo = zo, zr = sqrt(zr2), c = c, alternative = alternative)
        return(pC - level)
    }
    
    mylower <- 0
    myupper <- 1000
    c <- numeric()
    for(i in seq_len(length(zo))){
        target.l <- target(c = mylower, 
                           zo = zo[i], 
                           p = 1 - power, 
                           level = level,
                           designPrior = designPrior,
                           alternative = alternative)
        target.u <- target(c = myupper, 
                           zo = zo[i], 
                           p = 1 - power, 
                           level = level,
                           designPrior = designPrior,
                           alternative = alternative)
        if(sign(target.l) == sign(target.u)){
            if(sign(target.u) > 0)
                c[i] <- Inf
            else c[i] <- NA
        }
        else{
            c[i] <- uniroot(f = target, 
                            lower = mylower, 
                            upper = myupper, 
                            zo = zo[i], 
                            p = 1 - power, 
                            level = level,
                            designPrior = designPrior,
                            alternative = alternative)$root
        }
    }
    return(c)
}
    
