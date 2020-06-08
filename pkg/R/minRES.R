minRES <- function(zo, c, level = thresholdSceptical(level = 0.025, 
                                                   alternative = "one.sided", type = "golden")){
  mRV <- mapply(FUN = function(zo, c, level){
    if (!is.numeric(level) || (level <= 0 || level >= 1)) 
      stop("level must be numeric and in (0,1)!")
    zalphas <- p2z(level, alternative = "one.sided")
    K <- zo^2/zalphas^2
    if (zalphas > zo) {
      warning(paste("Replication success is not achievable at this level as", 
                    zo, " < ", round(p2z(level, alternative = "one.sided"), 
                                     3)))
      s <- NA
    } else {
    s <- ifelse(c<Inf, sqrt(1+c/(K-1))/(sqrt(K*c)), 1/sqrt(K*(K-1)))
    }
    return(s)
  }, zo, c, level)
  return(mRV)
}
