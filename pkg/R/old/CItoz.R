
CI.to.SE <- function(lower, upper, conf.level=0.95, ratio=FALSE){
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

CI.to.estimate <- function(lower, upper, ratio=FALSE, antilog=FALSE){
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


CI.to.z <- function(lower, upper, conf.level=0.95, ratio=FALSE){
    stopifnot(length(lower)==length(upper))
    stopifnot(sum(lower >= upper)==0)
    estimate <- CI.to.estimate(lower, upper, ratio=ratio)
    se <- CI.to.SE(lower, upper, conf.level=conf.level, ratio=ratio)
    z <- estimate/se
    return(z)
}

CI.to.p <- function(lower, upper, conf.level=0.95, ratio=FALSE,
                    alternative="two.sided"){
    z <- CI.to.z(lower, upper, conf.level=conf.level, ratio=ratio)
    p <- z.to.p(z, alternative=alternative)
    return(p)
}
