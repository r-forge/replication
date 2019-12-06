
###############################################
## Target function for calculating required
## sample size using powerClassical
###############################################

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
    c <- numeric()
    for(i in seq_len(length(to))){
        if(designPrior == "conditional"){
            u <- qnorm(p = power)
            v <- p2t(level, alternative = alternative)
            c <- (u + v)^2*(1/(shrinkage*to))^2
        }
        if(designPrior %in% c("predictive", "EB")){
            if (designPrior == "predictive") s <- shrinkage
            if (designPrior == "EB") s <- pmax(1 - (1 + d)/to[i]^2, 0)
            power.limit <- pnorm(sqrt(1/(s*(1 + d) + d))*s*abs(to[i]))
            if (power > power.limit) {
                power.limit.r <- floor(power.limit * 1000)/1000
                stop(paste("power too large, power should not exceed",
                           power.limit.r))
            }
            n.l <- 0
            n.u <- 100
            target.l <- ClassicalTarget(c = n.l, 
                                        to = to[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            target.u <- ClassicalTarget(c = n.u, to = to[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative,
                                        shrinkage = shrinkage,
                                        designPrior = designPrior)
            if (sign(target.l) == sign(target.u)) 
                c[i] <- NA
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
