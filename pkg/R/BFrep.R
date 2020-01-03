BFrep <- function(po = NULL, 
                  to = p2t(po, alternative = alternative),
                  pr = NULL,
                  tr = p2t(pr, alternative = alternative),
                  c = 1,
                  alternative = "two.sided"){
    # sanity checks
    if (min(c, na.rm = TRUE) < 0)
        stop("c must be larger than 0")

    # compute BF for H0: theta = 0 vs. H1: theta ~ N(hat(theta)_o, sigma^2_o)
    bf <- dnorm(x = tr, mean = 0, sd = 1)/
          dnorm(x = tr, mean = to*sqrt(c), sd = sqrt(c + 1))
    return(bf)
}