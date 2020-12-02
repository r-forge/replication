## compute project power of sceptical p-value for specified 
## level of replication success, variance ratio, alternative hypothesis,
## type-I error and power of original study
## TODO need to implement alternative "less" and "greater"

PPpSceptical <- function(level, c, alpha, power, alternative = "one.sided") {  
  
  ## vectorize function in all arguments
  ppV <- mapply(FUN = function(level, c, alpha, power, alternative) {
    ## sanity checks
    if (!(alternative %in% c("one.sided", "two.sided")))
      stop('alternative must be either "one.sided" or "two.sided"')
    if (!is.numeric(c) || c < 0)
      stop("c must be numeric and larger than 0")
    if (!is.numeric(level) || (level <= 0 || level >= 1))
      stop("level must be numeric and in (0, 1)!")
    if (!is.numeric(alpha) || (alpha <= 0 || alpha >= 1))
      stop("alpha must be numeric and in (0, 1)!")
    if (!is.numeric(power) || (power <= 0 || power >= 1))
      stop("power must be numeric and in (0, 1)!")
    
    ## compute normal quantile corresponding to level
    zas <- p2z(p = level, alternative = alternative)
    
    ## compute mean based on alpha and power 
    mu <- p2z(p = alpha, alternative = alternative) + stats::qnorm(p = power)
    
    ## compute project power with numerical integration
    if (alternative == "two.sided") {
      ## define function to integrate over zo
      intFun <- function(zo) {
        ## compute minimal zr to achieve replication success given zo and level
        K <- zo^2/zas^2
        zrmin <- zas*sqrt(1 + c/(K - 1))
        
        ## compute integrand
        ifelse(sign(zo) == 1,
           ## on positive side of plane (zo, zr > 0): P(|zr| >= zrmin)*dnorm(zo)
           (stats::pnorm(q = zrmin, mean = sqrt(c)*mu, lower.tail = FALSE) +
            stats::pnorm(q = -zrmin, mean = sqrt(c)*mu, lower.tail = TRUE))*
           stats::dnorm(x = zo, mean = mu),
           
           ## on negative side of plane (zo, zr < 0): P(zr <= -zrmin)*dnorm(zo)
           (stats::pnorm(q = zrmin, mean = sqrt(c)*mu, lower.tail = FALSE) +
            stats::pnorm(q = -zrmin, mean = sqrt(c)*mu, lower.tail = TRUE))*
           stats::dnorm(x = zo, mean = mu)
        )
      } 
    }
    
    if (alternative == "one.sided") {
      # define function to integrate over zo
      intFun <- function(zo) {
        ## compute minimal zr to achieve replication success given zo and level
        K <- zo^2/zas^2
        zrmin <- zas*sqrt(1 + c/(K - 1))
        
        ## compute integrand
        ifelse(sign(zo) == 1,
          ## on positive side of plane (zo, zr > 0): P(zr >= zrmin)*dnorm(zo)
          stats::pnorm(q = zrmin, mean = sqrt(c)*mu, lower.tail = FALSE)*
          stats::dnorm(x = zo, mean = mu),
        
          ## on negative side of plane (zo, zr < 0): P(zr <= -zrmin)*dnorm(zo)
          ## (will be very small for large mu)
          stats::pnorm(q = -zrmin, mean = sqrt(c)*mu, lower.tail = TRUE)*
          stats::dnorm(x = zo, mean = mu)
          )
      }
    }
    
    if (alternative %in% c("one.sided", "two.sided")) {
      ## integrate zo, zr over region where replication succcess possible
      pp <- stats::integrate(f = intFun, lower = zas, upper = Inf)$value +
            stats::integrate(f = intFun, lower = -Inf, upper = -zas)$value
      return(pp)
    }

    
    return(pp)
  }, level, c, alpha, power, alternative)
  
  return(ppV)
}
