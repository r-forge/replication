# Target function for calculating required sample size using uniroot
ClassicalTarget <- function(c, to, level = 0.05, power, alternative = alternative,
                            d = 0, shrinkage = 1, designPrior){
    term <- powerSignificance(to = to, c = c, level = level,
                              designPrior = designPrior,
                              alternative = alternative,
                              d = d,
                              shrinkage = shrinkage)
    return(term - power)
}

sampleSizeSignificance <- function(po = NULL,
                                   to = p2t(po, alternative = alternative),
                                   power,
                                   level = 0.05,
                                   designPrior = "conditional",
                                   alternative = "two.sided",
                                   d = 0,
                                   shrinkage = 1){
    # sanity checks
    if (!(designPrior %in% c("conditional", "predictive", "EB")))
        stop('designPrior must be either "conditional", "predictive", or "EB"')
    if (power <= 0 | power >= 1) 
        stop("power must be in (0, 1)")
    if (level <= 0 | level >= 1) 
        stop("level must be in (0, 1)")
    if (min(d, na.rm = TRUE) < 0)
        stop("d cannot be negative")
    if ((min(shrinkage, na.rm = TRUE) < 0 || max(shrinkage, na.rm = TRUE) > 1)) 
        stop("shrinkage must be in [0, 1]")
    
    c <- numeric()
    for(i in seq_len(length(to))){
        
        # for conditional designPrior use analytical solution
        if(designPrior == "conditional"){
            u <- qnorm(p = power)
            v <- p2t(level, alternative = alternative)
            c <- (u + v)^2*(1/(shrinkage*to))^2
        }
        
        # for predictive and EB designPrior use uniroot
        if(designPrior %in% c("predictive", "EB")){
            
            # compute upper bound of power
            if (designPrior == "predictive") s <- shrinkage
            if (designPrior == "EB") s <- pmax(1 - (1 + d)/to[i]^2, 0)
            power.limit <- pnorm(sqrt(1/(s*(1 + d) + d))*s*abs(to[i]))
            if (power > power.limit) {
                power.limit.r <- floor(power.limit * 1000)/1000
                stop(paste("power too large, power should not exceed",
                           power.limit.r))
            }
            
            # check whether desired power can be achieved for max c = 100
            n.l <- 0
            n.u <- 100
            target.l <- ClassicalTarget(c = n.l, 
                                        to = to[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        d = d,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            target.u <- ClassicalTarget(c = n.u, 
                                        to = to[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        d = d,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            if (sign(target.l) == sign(target.u)) 
                c[i] <- NA
            
            # determine c to achieve desired power
            else c[i] <- uniroot(ClassicalTarget, 
                                 lower = n.l, 
                                 upper = n.u,
                                 to = to[i],
                                 level = level,
                                 power = power, 
                                 alternative = alternative,
                                 d = d,
                                 shrinkage = shrinkage,
                                 designPrior = designPrior)$root
        }
    }
    return(c)
}
