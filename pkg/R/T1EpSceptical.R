## compute type-I error of sceptical p-value for specified 
## level of replication success, variance ratio, and 
## alternative hypothesis
## TODO need to implement alternative "less" and "greater"

T1EpSceptical <- function(level, c, alternative = "one.sided") {  
  
  ## vectorize function in all arguments
  t1errV <- mapply(FUN = function(level, c, alternative) {
    ## sanity checks
    if (!(alternative %in% c("one.sided", "two.sided")))
      stop('alternative must be either "one.sided" or "two.sided"')
    if (!is.numeric(c) || c < 0)
      stop("c must be numeric and larger than 0")
    if (!is.numeric(level) || (level <= 0 || level >= 1))
      stop("level must be numeric and in (0, 1)!")
    
    ## compute normal quantile corresponding to level
    zas <- p2z(p = level, alternative = alternative)
    
    if (alternative == "two.sided") {
      ## if c = 1 compute analytically
      if (c == 1) {
        t1err <- 2*(1 - stats::pnorm(q = 2*zas))
        return(t1err)
      } else {  ## if c != 1 use numerical integration
        
        ## define function to integrate over zo from zas to Infty
        intFun <- function(zo) {
          K <- zo^2/zas^2
          ## compute minimal zr to achieve replication success given zo and level
          zrmin <- zas*sqrt(1 + c/(K - 1))
          ## return integrand: P(|zr| >= zrmin)*dnorm(zo)
          2*(1 - stats::pnorm(q = zrmin))*stats::dnorm(x = zo)
        } 
      }
    }
    
    if (alternative == "one.sided") {
      ## if c = 1 compute analytically
      if (c == 1) {
        t1err <- 1 - stats::pnorm(q = 2*zas)
        return(t1err)
      } else { ## if c != 1 use numerical integration
        
        # define function to integrate over zo from zas to Infty
        intFun <- function(zo) {
          K <- zo^2/zas^2
          ## compute minimal zr to achieve replication success given zo and level
          zrmin <- zas*sqrt(1 + c/(K - 1))
          ## compute integrand: P(zr >= zrmin)*dnorm(zo)
          (1 - stats::pnorm(q = zrmin))*stats::dnorm(x = zo)
        }
      }
    }
    
    if (alternative %in% c("one.sided", "two.sided")) {
      ## the integral is symmetric around zero for "one.sided" and "two.sided" 
      ## so we can multiply the integral from zas to Infty by 2
      t1err <- 2*stats::integrate(f = intFun, lower = zas, upper = Inf)$value
      return(t1err)
    }
    
    return(t1err)
  }, level, c, alternative)
  
  return(t1errV)
}
