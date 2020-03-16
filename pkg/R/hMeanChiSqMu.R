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
    break.p <- res/(2^n)
    if(alternative=="greater")
        res <- ifelse(check.greater, break.p, paste(">", format(break.p, scientific = FALSE)))
    if(alternative=="less")
        res <- ifelse(check.less, break.p, paste(">", format(break.p, scientific = FALSE)))
    if(alternative=="two.sided")
        res <- ifelse((check.greater|check.less), 2*break.p, paste(">", format(2*break.p, scientific = FALSE)))
    return(res)
}



