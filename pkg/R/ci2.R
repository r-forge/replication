
ci2se <- function(lower, upper, conf.level=0.95, ratio=FALSE){
    stopifnot(length(lower)==length(upper))
    stopifnot(sum(lower >= upper)==0)
    level <- 1-conf.level
    q <- qnorm(1-level/2, lower.tail=TRUE)
    if(ratio==TRUE){
        stopifnot(sum(lower <= 0)==0)
        lower <- log(lower)
        upper <- log(upper)
    }
    se <- (upper-lower)/(2*q)
    return(se)
}

ci2estimate <- function(lower, upper, ratio=FALSE, antilog=FALSE){
    stopifnot(length(lower)==length(upper))
    stopifnot(sum(lower >= upper)==0)
    if(ratio==TRUE){
        stopifnot(sum(lower <= 0)==0)
        lower <- log(lower)
        upper <- log(upper)
    }
    res <- (lower+upper)/2
    if((ratio==TRUE)&(antilog==TRUE))
        res <- exp(res)
    return(res)
}


ci2t <- function(lower, upper, conf.level=0.95, ratio=FALSE){
    stopifnot(length(lower)==length(upper))
    stopifnot(sum(lower >= upper)==0)
    estimate <- ci2estimate(lower, upper, ratio=ratio)
    se <- ci2se(lower, upper, conf.level=conf.level, ratio=ratio)
    t <- estimate/se
    return(t)
}

ci2p <- function(lower, upper, conf.level=0.95, ratio=FALSE,
                    alternative="two.sided"){
    t <- ci2t(lower, upper, conf.level=conf.level, ratio=ratio)
    p <- t2p(t, alternative=alternative)
    return(p)
}
