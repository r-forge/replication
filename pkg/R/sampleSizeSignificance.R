
###############################################
## Target function for calculating required
## sample size using powerClassical
###############################################

ClassicalTarget <- function(c, to, level=0.05, power, alternative=alternative,
                            shrinkage=FALSE, d=0){
    term <- powerSignificance(to=to, c=c, level=level,
                              designPrior="predictive",
                              alternative=alternative,
                              shrinkage=shrinkage,
                              d=d)
    return(term - power)
}

sampleSizeSignificance <- function(po=NULL,
                                   to=p2t(po, alternative=alternative),
                                   power,
                                   level = 0.05,
                                   designPrior="conditional",
                                   alternative="two.sided",
                                   shrinkage=FALSE,
                                   d=0){
    c <- numeric()
    for(i in seq_len(length(to))){
        if(designPrior=="conditional"){
            u <- qnorm(p = power)
            v <- p2t(level, alternative=alternative)
            c <- (u+v)^2*(1/to)^2
        }
        if(designPrior=="predictive"){
            power.limit <- pnorm(abs(to[i]))
            if (power > power.limit) {
                power.limit.r <- floor(power.limit * 1000)/1000
                stop(paste("power too large, power should not exceed",
                           power.limit.r))
            }
            nstart <- sampleSizeSignificance(to = to[i], 
                                             power = power,
                                             level = level,
                                             alternative = alternative)
            n.l <- 0
            n.u <- 100
            target.l <- ClassicalTarget(c = n.l, 
                                        to = to[i],
                                        level = level,
                                        power = power,
                                        alternative=alternative)
            target.u <- ClassicalTarget(c=n.u, to=to[i],
                                        level = level,
                                        power = power,
                                        alternative = alternative)
            if (sign(target.l) == sign(target.u)) 
                c[i] <- NA
            else c[i] <- uniroot(ClassicalTarget, lower = n.l, upper = n.u,
                                 to=to[i],
                                 level = level,
                                 power = power, 
                                 alternative = alternative,
                                 shrinkage = shrinkage,
                                 d = d)$root
        }
    }
    return(c)
}
