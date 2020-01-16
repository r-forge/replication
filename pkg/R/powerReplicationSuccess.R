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
                                    alternative = "one.sided"){
    
    targetPower <- function(power, 
                            zo, 
                            c, 
                            level,
                            designPrior, 
                            alternative){
        zr2 <- zr2.quantile(zo = zo, c = c, p = 1 - power, 
                            designPrior = designPrior)
        pC <- pSceptical(zo = zo, zr = sqrt(zr2), c = c, 
                         alternative = alternative)
        return(pC - level)
    }
    
    # check if original study was not significant, then power is zero
    zo <- abs(zo)
    p <- z2p(z = zo, alternative = alternative)
    if((length(c) == 1) & (length(p) > 1))
        c <- rep(c, length(p))
    eps <- 1e-5
    mylower <- eps
    myupper <- 1 - eps
    res <- numeric()
    for(i in seq_len(length(zo))){
        if(p[i] > level) res[i] <- 0
        else{
            # now the more interesting case
            target.l <- targetPower(power = mylower, 
                                    zo = zo[i], 
                                    c = c[i], 
                                    level = level,
                                    designPrior = designPrior,
                                    alternative = alternative)
            target.u <- targetPower(power = myupper, 
                                    zo = zo[i], 
                                    c = c[i], 
                                    level = level,
                                    designPrior = designPrior,
                                    alternative = alternative)
            if(sign(target.l) == sign(target.u)){
                if((sign(target.l) >= 0) & (sign(target.u) >= 0))
                    res[i] <- 0
                if((sign(target.l) < 0) & (sign(target.u) < 0))
                    res[i] <- 1
            }
            else{
                res[i] <- uniroot(f = targetPower, 
                                  lower = mylower, 
                                  upper = myupper,
                                  zo = zo[i], 
                                  c = c[i], 
                                  level = level,
                                  designPrior = designPrior,
                                  alternative = alternative)$root
            }
        }
    }
    return(res)
}
