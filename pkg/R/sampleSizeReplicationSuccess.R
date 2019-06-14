
target <- function(to, c, p, level=0.05, designPrior, alternative=alternative){
    tr2 <- tr2.quantile(to=to, c=c, p=p, designPrior=designPrior)
    pC <- pSceptical(to=to, tr=sqrt(tr2), c=c, alternative=alternative)
    return(pC-level)
}

sampleSizeReplicationSuccess <- function(po=NULL,
                                         to=p2t(po, alternative=alternative),
                                         power, level=0.05,
                                         designPrior="conditional",
                                         alternative="two.sided"){
    if(min(power) < 0 || max(power) >1)
        stop("The power must lie in (0,1)!")
    mylower <- 0
    myupper <- 1000
    c <- numeric()
    for(i in seq_len(length(to))){
        target.l <- target(mylower, to=to[i], p=1-power, level=level,
                           designPrior=designPrior,
                           alternative=alternative)
        target.u <- target(myupper, to=to[i], p=1-power, level=level,
                           designPrior=designPrior,
                           alternative=alternative)
        if(sign(target.l)==sign(target.u)){
            if(sign(target.u)>0)
                c[i] <- Inf
            else c[i] <- NA
        }
        else{
            c[i] <- uniroot(target, lower=mylower, upper=myupper, 
                            to=to[i], p=1-power, level=level,
                            designPrior=designPrior,
                            alternative=alternative)$root
        }
    }
    return(c)
}
    
