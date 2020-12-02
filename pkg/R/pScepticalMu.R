## Computes the sceptical p-value when the null hypothesis is H0: theta = mu. 

pScepticalMu <- function(thetao, 
                         thetar, 
                         se_thetao, 
                         se_thetar, 
                         mu, 
                         alternative = "two.sided",
                         type = "golden"){
  resV <- mapply(FUN = function(thetao, thetar, se_thetao,  se_thetar, mu, 
                                alternative, type) {
    ## sanity checks
    if (!is.numeric(se_thetao) || se_thetao < 0)
      stop("se_thetao must be numeric and larger than 0")
    if (!is.numeric(se_thetar) || se_thetar < 0)
      stop("se_thetar must be numeric and larger than 0")
    
    ## compute z-statistics and c
    zo <- (thetao - mu)/se_thetao
    zr <- (thetar - mu)/se_thetar
    c <- se_thetao^2/se_thetar^2
    
    ## compute sceptical p-value
    pS <- pSceptical(zo = zo, zr = zr, c = c, alternative = alternative, 
                     type = type)
    return(pS)
  }, thetao, thetar, se_thetao, se_thetar, mu, alternative, type)
  return(resV)
}

