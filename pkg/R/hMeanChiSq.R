
## harmonic mean chi-squared test based on z-values

hMeanChiSq <- function(z, w=rep(1, length(z)), alternative="greater"){
    stopifnot(min(w)>0)
    if (!(alternative %in% c("greater", "less", "two.sided")))
        stop('alternative must be either "greater", "less" of "two.sided"')
    require(ReplicationSuccess)
    n <- length(z)
    zH2 <- sum(sqrt(w))^2/sum(w/z^2)
    res <- pchisq(zH2, df=1, lower.tail=FALSE)
    check.greater <- (min(z)>0)
    check.less <- (max(z)<0)
    if(alternative=="greater")
        res <- ifelse(check.greater, res/(2^n), NaN)
    if(alternative=="less")
        res <- ifelse(check.less, res/(2^n), NaN)
    if(alternative=="two.sided")
        res <- ifelse((check.greater|check.less), res/(2^(n-1)), NaN)
    return(res)
}

