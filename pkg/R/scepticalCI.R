## Compute the sceptical confidence interval

scepticalCI <- function(thetao, thetar, se_thetao, se_thetar, alphas){
  fun <- function(mu){
    pS <- pScepticalMu(thetao, thetar, se_thetao, se_thetar, mu, 
                       alternative = "two.sided")
    return(pS - alphas)
  }
  rootSolve::uniroot.all(fun, lower = -5, upper = 5)
  
}
