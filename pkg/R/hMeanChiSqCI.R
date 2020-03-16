## confidence interval based on harmonic mean chi-squared test

hMeanChiSqCI <- function(thetahat, se, w=rep(1, length(thetahat)),
                         alternative="two.sided", level=0.95){
    stopifnot(min(w)>0)
    stopifnot(min(se)>0)
    if (!(alternative %in% c("two.sided", "greater", "less")))
        stop('alternative must be either "two.sided", "greater" or "less"')
    stopifnot((level>0)|(level<1))
    ## target function to compute the limits of the CI
    target <- function(limit, thetahat, se, w=w, alternative=alternative, alpha){
        res <- hMeanChiSqMu(thetahat, se, w=w, mu=limit, alternative=alternative)-alpha
        return(res)
    }
    mini <- which.min(thetahat)
    maxi <- which.max(thetahat)
    mint <- thetahat[mini]
    maxt <- thetahat[maxi]
    minse <- se[mini]
    maxse <- se[maxi]
    alpha <- 1-level
    z <- -qnorm(alpha)
    eps <- 1e-6
    if(alternative=="two.sided"){
        lower <- uniroot(target, thetahat=thetahat, se=se, w=w, alternative=alternative,
                         alpha=alpha, lower=mint-5*z*minse, upper=mint-eps*minse)$root
        upper <- uniroot(target, thetahat=thetahat, se=se, w=w, alternative=alternative,
                         alpha=alpha, lower=maxt+eps*maxse, upper=maxt+5*z*maxse)$root
    }
    if(alternative=="greater"){
        lower <- uniroot(target, thetahat=thetahat, se=se, w=w, alternative=alternative,
                         alpha=alpha, lower=mint-5*z*minse, upper=mint-eps*minse)$root
        upper <- Inf
    }
    if(alternative=="less"){
        lower <- -Inf
        upper <- uniroot(target, thetahat=thetahat, se=se, w=w, alternative=alternative,
                         alpha=alpha, lower=maxt+eps*maxse, upper=maxt+5*z*maxse)$root
    }
    return(c(lower, upper))
}
