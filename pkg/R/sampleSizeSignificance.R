sampleSizeSignificance <- function(zo,
                                   power,
                                   level = 0.025,
                                   designPrior = "conditional",
                                   alternative = "one.sided",
                                   d = 0,
                                   shrinkage = 0){
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
    
    # s is 1 - shrinkage
    s <- 1 - shrinkage
    
    # for conditional designPrior use analytical solution
    if(designPrior == "conditional"){
        u <- qnorm(p = power)
        v <- p2z(p = level, alternative = alternative)
        c <- (u + v)^2*(1/(s*zo))^2
    }
    
    # Target function for calculating required sample size using uniroot
    ClassicalTarget <- function(c, 
                                zo, 
                                level, 
                                power, 
                                alternative,
                                d, 
                                shrinkage, 
                                designPrior){
        term <- powerSignificance(zo = zo, 
                                  c = c, 
                                  level = level,
                                  designPrior = designPrior,
                                  alternative = alternative,
                                  d = d,
                                  shrinkage = shrinkage)
        return(term - power)
    }
    
    # for predictive and EB designPrior use uniroot
    if(designPrior %in% c("predictive", "EB")){
        c <- numeric()
        for(i in seq_len(length(zo))){

            # compute upper bound of power
            if (designPrior == "EB") s <- pmax(1 - (1 + d)/zo[i]^2, 0)
            power.limit <- pnorm(sqrt(1/(s*(1 + d) + d))*s*abs(zo[i]))
            if (power > power.limit) {
                power.limit.r <- floor(power.limit * 1000)/1000
                # warning(paste("power too large, power should not exceed",
                #                power.limit.r,
                #               "for a zo of",
                #               zo,
                #               "\n"))
                c[i] <- NaN
            }
            
            # check whether desired power can be achieved for max c = 1000
            n.l <- 0
            n.u <- 1000
            target.l <- ClassicalTarget(c = n.l, 
                                        zo = zo[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        d = d,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            target.u <- ClassicalTarget(c = n.u, 
                                        zo = zo[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        d = d,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            if (sign(target.l) == sign(target.u)) 
                c[i] <- NaN
            
            # determine c to achieve desired power
            else c[i] <- uniroot(f = ClassicalTarget, 
                                 lower = n.l, 
                                 upper = n.u,
                                 zo = zo[i],
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
