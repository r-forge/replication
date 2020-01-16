
p2z <- function(p, 
                alternative = "two.sided"){
  
  if(any(!is.na(p)) && (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >1))
    stop("All elements of p must lie in (0,1]!")
  if (!(alternative %in% c("less", "greater", "two.sided", "one.sided")))
    stop('alternative must be either "less", "greater", "two.sided", or "one.sided"')
  
  if(alternative == "two.sided")
    z <- qnorm(p = p/2, lower.tail = FALSE)
  if (alternative == "less")
    z <- qnorm(p = p, lower.tail = TRUE)
  if (alternative == "greater" | alternative == "one.sided")
    z <- qnorm(p = p, lower.tail = FALSE)
  
  return(z)
}
