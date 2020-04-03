## Function to compute density of truncated normal
dnormtrunc <- function(x, mean = 0, sd = 1, a = -Inf, b = Inf, log = FALSE) {
    k <- diff(pnorm(q = c(a, b), mean = mean, sd = sd))
    I <- as.numeric(x >= a & x <= b)
    dens <- dnorm(x = x, mean = mean, sd = sd)*I/k
    if (log == TRUE) dens <- log(dens)
    return(dens)
}

## BF comparing H0: theta = 0 vs. H1: theta ~ N(hat(theta)_o, sigma_o^2) 
## for data from replication study hat(theta)_r ~ N(theta, sigma^2_r)
BFrep <- function(zo, zr, c, truncate = FALSE) {
    if (!is.numeric(c) || c < 0) 
        stop("c must be numeric and larger than zero")
    
    if (truncate == FALSE) {
        logbf <- dnorm(x = z_r, mean = 0, sd = 1, log = TRUE) -
            dnorm(x = z_r, mean = z_o*sqrt(c), sd = sqrt(1 + c), log = TRUE)
    } else {
        logbf <- dnorm(x = z_r, mean = 0, sd = 1, log = TRUE) -
            dnormtrunc(x = z_r, mean = z_o*sqrt(c), sd = sqrt(1 + c), log = TRUE,
                       a = ifelse(z_o >= 0, 0, -Inf), 
                       b = ifelse(z_o >= 0, Inf, 0))
    }
    return(exp(logbf))
}