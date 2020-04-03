pSceptical <- function(zo,
                       zr,
                       c, 
                       alternative = "one.sided"){
    
    ## vectorize function in all arguments
    resV <- mapply(FUN = function(zo, zr, c, alternative) {
        ## sanity checks
        if (!(alternative %in% c("one.sided", "two.sided")))
            stop('alternative must be either "one.sided" or "two.sided"')
        if (!is.numeric(c) || c < 0)
            stop("c must be numeric and larger than 0")
        
        z <- zSceptical(zo, zr, c)
        res <- z2p(z = z)
        if(alternative == "one.sided") {
            if(sign(zo) == sign(zr)) 
                res <- res/2
            else 
                res <- 1 - res/2
        }
        return(res)
    }, zo, zr, c, alternative)
    
    return(resV)
}
