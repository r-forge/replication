
p.to.t <- function(p, alternative="two.sided"){
  if(min(p) <=0 || max(p) >1)
    stop("All elements of p must lie in (0,1]!")
  if(alternative=="two.sided")
    t <- qnorm(p/2, lower.tail=FALSE)
  if(alternative=="one.sided")
    t <- qnorm(p, lower.tail=FALSE)
  return(t)
}
