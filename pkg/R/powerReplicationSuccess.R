
tr2.quantile <- function(to, c, p, designPrior){
    if(designPrior=="predictive"){
        lambda <- to^2/(1+1/c)
        factor <- c+1
    }
    if(designPrior=="conditional"){
        lambda <- c*to^2
        factor <- 1
    }
    # if (designPrior == "EB") {
    #     d <- 0 # heterogeneity, implement later
    #     s <- pmax(1 - (1 + d)/to^2, 0)
    #     lambda <- to^2*s^2*c/(s*c*(1 + d) + 1 + d*c)
    #     factor <- s*c*(1 + d) + 1 + d*c
    # }
    if(lambda < 100)
        res <- qchisq(p, df=1, ncp=lambda)
    else
        res <- qnorm(p, mean=sqrt(lambda), sd=1)^2
    return(factor*res)
}

targetPower <- function(power, to, c, level=0.05,
                        designPrior, alternative=alternative){
    tr2 <- tr2.quantile(to=to, c=c, p=1-power, designPrior=designPrior)
    pC <- pSceptical(to=to, tr=sqrt(tr2), c=c, alternative=alternative)
    return(pC-level)
}

powerReplicationSuccess <- function(po=NULL,
                                    to=p2t(po, alternative=alternative),
                                    c=1, level=0.05,
                                    designPrior="conditional",
                                    alternative="two.sided"){
    ## check if original study was not significant
    ## then power is zero
    to <- abs(to)
    p <- t2p(to, alternative=alternative)
    if((length(c)==1)&(length(p)>1))
        c <- rep(c, length(p))
    eps <- 1e-5
    mylower <- eps
    myupper <- 1-eps
    res <- numeric()
    for(i in seq_len(length(to))){
        if(p[i] > level) res[i] <- 0
        else{
            ## now the more interesting case
            target.l <- targetPower(mylower, to=to[i], c=c[i], 
                                    level=level,
                                    designPrior=designPrior,
                                    alternative=alternative)
            target.u <- targetPower(myupper, to=to[i], c=c[i], 
                                    level=level,
                                    designPrior=designPrior,
                                    alternative=alternative)
            if(sign(target.l)==sign(target.u)){
                if((sign(target.l)>=0) & (sign(target.u)>=0))
                    res[i] <- 0
                if((sign(target.l)<0) & (sign(target.u)<0))
                    res[i] <- 1
            }
            else{
                res[i] <- uniroot(targetPower, lower=mylower, upper=myupper,
                                  to=to[i], c=c[i], level=level,
                                  designPrior=designPrior,
                                  alternative=alternative)$root
            }
        }
    }
    return(res)
}
