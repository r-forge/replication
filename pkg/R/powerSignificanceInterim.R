powerSignificanceInterim <- function(to, ti, c = 1, f = 1 / 2,
                                     level = 0.05,
                                     designPrior = "conditional",
                                     analysisPrior = "flat",
                                     alternative = "two.sided",
                                     shrinkage = 1) 
{
  if (!(designPrior %in% c("conditional", "predictive", "flat"))) 
    stop("designPrior must be either \"conditional\", \"predictive\", or \"flat\"")
  if (!(analysisPrior %in% c("flat", "original"))) 
    stop("analysisPrior must be either \"flat\" or \"original\"")
  if (min(c, na.rm = TRUE) < 0) 
    stop("c must be larger than 0")
  if ((min(f, na.rm = TRUE) < 0 || max(f, na.rm = TRUE) > 
       1)) 
    stop("f must be in [0, 1]")
  if ((min(shrinkage, na.rm = TRUE) < 0 || max(shrinkage, na.rm = TRUE) > 
       1)) 
    stop("shrinkage must be in [0, 1]")
  v <- p2t(level, alternative = alternative)
  to <- shrinkage * to
  if (designPrior == "conditional")
    if (analysisPrior == "original")
      return(NA) ## For now, we are not interested in the case where the design prior is conditional and the analysis prior normal.
  pSig <-
    pnorm(to * sqrt(c * (1 / f - 1)) + ti / (sqrt(1 / f - 1)) - sqrt(1 / (1 -
                                                                            f)) * v)
  if (designPrior == "predictive") {
    if (analysisPrior == "original") {
      term1 <- sqrt(1 + (c * (1 - f)) / (f * (c + 1)))
      term2 <- sqrt(f / (c * (1 - f))) * to
      term3 <- sqrt(f / (1 - f)) * ti
      term4 <- sqrt((f * (c + 1)) / (c * (1 - f))) * v
      
      pSig = pnorm(term1 * (term2 + term3) - term4)
    }
    else if (analysisPrior == "flat") {
      term1 <- sqrt(((1 - f) * c) / ((c + 1) * (f + c))) * to
      term2 <- sqrt((f + c) / ((1 - f) * (c + 1))) * ti
      term3 <- sqrt((f * (c + 1)) / ((f + c) * (1 - f))) * v
      pSig <- pnorm(term1 + term2 - term3)
    }
  }
  if (designPrior == "flat")
    if (analysisPrior == "flat") {
      pSig <- pnorm((ti - sqrt(f) * v) / sqrt(1 - f))
    } else if (analysisPrior == "original")
      return(NA)
  if (designPrior == "conditional") {
    if (analysisPrior == "original")
      return(NA)
    pSig <-
      pnorm(ti * sqrt(c * (1 / f - 1)) + ti / (sqrt(1 / f - 1))
            - sqrt(1 / (1 - f)) * v)
  }
  
  return(pSig)
}
