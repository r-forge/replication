sampleSizeRES <- function(zo, s = 1, 
                          level = thresholdSceptical(level = 0.025, 
                                                     alternative = "one.sided", 
                                                     type = "golden")){
  sSRV <- mapply(FUN = function(zo, s, level) {
    if (!is.numeric(s)) 
      stop("s must be numeric")
    if (!is.numeric(level) || (level <= 0 || level >= 1)) 
      stop("level must be numeric and in (0,1)!")
    zalphas <- p2z(level, alternative="one.sided")
    K <- zo^2/zalphas^2
    denom <- s^2*K - 1/(K-1)
    if (zalphas > zo){
      warning(paste("Replication success is not achievable at this level as", 
                    zo, " < ",  round(p2z(level, alternative = "one.sided"), 3)))
      c <- NA
    } else { 
    c <- ifelse(denom > 0, 1/denom, NA)  # denom < 0: smin > s_\infty
    }
    return(c)
    
    }, zo, s, level)
  return(sSRV)
}



