## harmonic mean chi-squared test based on estimates and standard errors

hMeanChiSqMu <- function(thetahat, se, w=rep(1, length(thetahat)), mu=0, alternative="greater"){
    stopifnot(min(w)>0)
    stopifnot(min(se)>0)
    if (!(alternative %in% c("greater", "less", "two.sided")))
        stop('alternative must be either "greater", "less" of "two.sided"')
    require(ReplicationSuccess)
    n <- length(thetahat)
    z <- (thetahat-mu)/se
    zH2 <- sum(sqrt(w))^2/sum(w/z^2)
    res <- pchisq(zH2, df=1, lower.tail=FALSE)
    check.greater <- (min(thetahat-mu)>0)
    check.less <- (max(thetahat-mu)<0)
    
    if(alternative=="greater")
        res <- ifelse(check.greater, res/(2^n), NaN)
    if(alternative=="less")
        res <- ifelse(check.less, res/(2^n), NaN)
    if(alternative=="two.sided")
        res <- ifelse((check.greater|check.less), res/(2^(n-1)), NaN)
    return(res)
}



