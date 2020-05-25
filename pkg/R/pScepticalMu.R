## Computes the sceptical p-value when the null hypothesis is H0: theta = mu. 

pScepticalMu = function(thetao, 
                        thetar, 
                        se_thetao, 
                        se_thetar, 
                        mu, 
                        alternative = "two.sided"){
  resV <- mapply(FUN = function(thetao, thetar, se_thetao,  se_thetar, mu, alternative) {
    if (!(alternative %in% c("one.sided", "two.sided"))) 
      stop("alternative must be either \"one.sided\" or \"two.sided\"")
    zo = (thetao-mu)/se_thetao
    zr = (thetar-mu)/se_thetar
    pS = pSceptical(zo, zr, se_thetao^2/se_thetar^2, alternative = alternative)
    if (alternative == "one.sided") {
      if (sign(zo) == sign(zr)) 
        res <- res/2
      else res <- 1 - res/2
    }
    return(pS)
  }, thetao, thetar, se_thetao, se_thetar, mu, alternative)
  return(resV)
}

