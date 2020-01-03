predictionInterval <- function(po = NULL,
                               to = p2t(po, alternative = alternative),
                               c = 1,
                               conf.level = 0.95,
                               designPrior = "predictive",
                               alternative = "two.sided",
                               d = 0){
    # sanity checks
    if (!(designPrior %in% c("conditional", "predictive", "EB")))
        stop('designPrior must be either "conditional", "predictive", or "EB"')
    if (min(c, na.rm = TRUE) < 0)
        stop("c must be larger than 0")
    if (min(d, na.rm = TRUE) < 0)
        stop("d cannot be negative")
    
    # determine parameters of predictive distribution of tr
    if(designPrior == "conditional"){
        mu <- to*sqrt(c)
        sigma <- 1
    }
    if(designPrior == "predictive"){
        mu <- to*sqrt(c)
        sigma <- sqrt(c + 1 + 2*d*c)
    }
    if (designPrior == "EB"){
        s <- pmax(1 - (1 + d)/to^2, 0)
        mu <- s*to*sqrt(c)
        sigma <- sqrt(s*c*(1 + d) + 1 + d*c)
    }
    
    # compute prediction interval
    lower <- qnorm(p = (1 - conf.level)/2, mean = mu, sd = sigma)
    upper <- qnorm(p = (1 + conf.level)/2, mean = mu, sd = sigma)
    result <- data.frame(lower = lower, mean = mu, upper = upper)
    return(result)
}
