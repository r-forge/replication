
p2t <- function(p, alternative="two.sided"){
  if(any(!is.na(p)) && (min(p, na.rm=TRUE) <=0 || max(p, na.rm=TRUE) >1))
    stop("All elements of p must lie in (0,1]!")
  if(alternative=="two.sided")
    t <- qnorm(p/2, lower.tail=FALSE)
  # if(alternative=="one.sided")
  #   t <- qnorm(p, lower.tail=FALSE)
  if (alternative == "less")
    t <- qnorm(p, lower.tail = TRUE)
    if (alternative == "greater" | alternative == "one.sided")
    t <- qnorm(p = p, lower.tail = FALSE)
  return(t)
}
