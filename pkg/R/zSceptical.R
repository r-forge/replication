zSceptical <- function(zo, 
                       zr, 
                       c){
  ## arithmetic mean
  arit.mean <- function(x, y)
    return((x + y)/2)
  ## harmonic mean
  harm.mean <- function(x, y)
    return(2/(1/x + 1/y))
  
  z2H <- harm.mean(zo^2, zr^2)
  z2A <- arit.mean(zo^2, zr^2)
  if(length(c) == 1){
    if(c == 1)
      z2 <- z2H/2
    else
      z2 <- (sqrt(z2A*(z2A + (c - 1)*z2H)) - z2A)/(c - 1)
  }
  if(length(c) > 1)
    z2 <- ifelse(c == 1, z2H/2, 
                 (sqrt(z2A*(z2A + (c - 1)*z2H)) - z2A)/(c - 1))
  return(sqrt(z2))
}