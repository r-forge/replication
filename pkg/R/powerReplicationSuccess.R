# function that returns z_r^2 quantile for given z_o, c
zr2.quantile <- function(zo, 
                         c, 
                         p, 
                         designPrior){
    
    if(designPrior == "predictive"){
        lambda <- zo^2/(1 + 1/c)
        factor <- c + 1 
        # factor <- c + 1 + 2*d
    }
    if(designPrior == "conditional"){
        lambda <- c*zo^2
        factor <- 1
    }
    # if (designPrior == "EB") {
    #     d <- 0 # heterogeneity, implement later
    #     s <- pmax(1 - (1 + d)/zo^2, 0)
    #     lambda <- zo^2*s^2*c/(s*c*(1 + d) + 1 + d*c)
    #     factor <- s*c*(1 + d) + 1 + d*c
    # }
    if(lambda < 100)
        res <- qchisq(p = p, df = 1, ncp = lambda)
    else
        res <- qnorm(p = p, mean = sqrt(lambda), sd = 1)^2
    return(factor*res)
}

powerReplicationSuccess <- function(zo,
                                    c = 1, 
                                    level = 0.025,
                                    designPrior = "conditional",
                                    alternative = "one.sided",
                                    type = "golden"){
    
    targetPower <- function(power, zo, c, level, designPrior, alternative, type){
        zr2 <- zr2.quantile(zo = zo, c = c, p = 1 - power, 
                            designPrior = designPrior)
        pC <- pSceptical(zo = zo, zr = sqrt(zr2), c = c,
                         alternative = alternative, type = type)
        return(pC - level)
        # pC <- pSceptical(zo = zo, zr = sqrt(zr2), c = c, 
        #                  alternative = alternative)
        # return(pC - levelSceptical(level = level, alternative = alternative,
        #                                type = type))
    }
    
    # vectorize function in all arguments
    resV <- mapply(FUN = function(zo, c, level, designPrior, alternative, type) {
        
        # sanity checks
        if (is.na(zo))
            return(NA)
        if (!(designPrior %in% c("conditional", "predictive")))
            stop('designPrior must be either "conditional" or "predictive"')
        if (!is.numeric(c) || c < 0)
            stop("c must be numeric and larger than 0")
        if (!is.numeric(level) || (level <= 0 || level >= 1))
            stop("level must be numeric and in (0, 1)!")
        
        # check if original study was not significant, then power is zero
        zo <- abs(zo)
        p <- z2p(z = zo, alternative = alternative)
        eps <- 1e-5
        mylower <- eps
        myupper <- 1 - eps
        
        if (p > levelSceptical(level = level, 
                                   alternative = alternative, 
                                   type = type)) res <- 0
        else {
            target.l <- targetPower(power = mylower, 
                                    zo = zo, 
                                    c = c, 
                                    level = level,
                                    designPrior = designPrior,
                                    alternative = alternative,
                                    type = type)
            target.u <- targetPower(power = myupper, 
                                    zo = zo, 
                                    c = c, 
                                    level = level,
                                    designPrior = designPrior,
                                    alternative = alternative,
                                    type = type)
            if (sign(target.l) == sign(target.u)) {
                if ((sign(target.l) >= 0) & (sign(target.u) >= 0))
                    res <- 0
                if ((sign(target.l) < 0) & (sign(target.u) < 0))
                    res <- 1
            }
            else {
                res <- uniroot(f = targetPower, 
                               lower = mylower, 
                               upper = myupper,
                               zo = zo, 
                               c = c, 
                               level = level,
                               designPrior = designPrior,
                               alternative = alternative,
                               type = type)$root
            }
        }
        return(res)
    }, zo, c, level, designPrior, alternative, type)
    
    return(resV)
}
