## compute type-I error of sceptical p-value for specified 
## level of replication success, variance ratio, and 
## alternative hypothesis
## TODO need to implement alternative "less" and "greater"
## LH 8.2.2021: I have added type="golden" as default
## LH 8.2.2021: I have added alternative "greater" and "less"

T1EpSceptical <- function(level, c, alternative = "one.sided", type = "golden") {  
    
    ## vectorize function in all arguments
    t1errV <- mapply(FUN = function(level, c, alternative, type) {
        ## sanity checks
        if (!(alternative %in% c("one.sided", "two.sided", "greater", "less")))
            stop('alternative must be either "one.sided", "two.sided", "greater" or "less"')
        ## if (!(alternative %in% c("one.sided", "two.sided")))
        ##     stop('alternative must be either "one.sided" or "two.sided"')
        if (!is.numeric(c) || c < 0)
            stop("c must be numeric and larger than 0")
        if (!is.numeric(level) || (level <= 0 || level >= 1))
            stop("level must be numeric and in (0, 1)!")
        if (!(type %in% c("nominal", "liberal", "controlled", "golden")))
            stop('type must be either "nominal", "liberal", "controlled", or "golden"')
        
        ## compute normal quantile corresponding to level and type
        if (alternative %in% c("one.sided", "two.sided")) {
            alphas <- levelSceptical(level = level, 
                                     alternative = alternative, 
                                     type = type)
            zas <- p2z(alphas, alternative = alternative)
        }
        ## LH: not very elegant, levelSceptical needs to be updated 
        if (alternative %in% c("greater", "less")) {
            alphas <- levelSceptical(level = level, 
                                     alternative = "one.sided", 
                                     type = type)
            zas <- p2z(alphas, alternative = "one.sided")
        }
        
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
        if (alternative %in% c("greater", "less")) {
            ## if c = 1 compute analytically
            if (c == 1) {
                t1err <- (1 - stats::pnorm(q = 2*zas))/2
                return(t1err)
            } else { ## if c != 1 use numerical integration
                
                ## define function to integrate over zo from zas to Infty
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
        if (alternative %in% c("greater", "less")) {
            t1err <- stats::integrate(f = intFun, lower = zas, upper = Inf)$value
            return(t1err)
        }
        
        return(t1err)
    }, level, c, alternative, type)
    
    return(t1errV)
}
