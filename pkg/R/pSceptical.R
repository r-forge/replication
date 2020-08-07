pSceptical <- function(zo,
                       zr,
                       c, 
                       alternative = "one.sided",
                       type="nominal"){
    
    ## vectorize function in all arguments
    resV <- mapply(FUN = function(zo, zr, c, alternative, type) {
        ## sanity checks
        if (!(alternative %in% c("one.sided", "two.sided")))
            stop('alternative must be either "one.sided" or "two.sided"')
        if (!is.numeric(c) || c < 0)
            stop("c must be numeric and larger than 0")
        if (!(type %in% c("nominal", "liberal", "controlled", "golden")))
            stop('type must be either "nominal", "liberal", "controlled", or "golden"')
        
        z <- zSceptical(zo, zr, c)
        if(type=="nominal")
            result <- z
        if(type=="liberal"){
            result <- z*sqrt(2)
        }
        if(type=="controlled"){
            result <- p2z(sqrt((1-pnorm(2*z))/2), alternative = "greater")
        }
        if(type=="golden"){
            ## golden ratio 
            phi <- (sqrt(5)+1)/2  
            result <- z*sqrt(phi)
        }
        res <- z2p(z = result)
        if(alternative == "one.sided") {
            if(sign(zo) == sign(zr)) 
                res <- res/2
            else 
                res <- 1 - res/2
        }
        return(res)
    }, zo, zr, c, alternative, type)
    
    return(resV)
}
