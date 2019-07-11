## Replace this function with ReplicationSuccess::p2t
# p.to.z <- function(p)
#   return(qnorm(1-p/2))

## Replace functions with ReplicationSuccess::powerSignificance
# StandardPower <- function(n, delta, sd, sig.level=0.05, type="two.sample", strict=FALSE){
#   stopifnot(type=="two.sample" | type=="one.sample")
#   if(type=="two.sample")
#     sd <- sqrt(2)*sd
#   v <- qnorm(p = sig.level/2)
#   if(strict==FALSE)
#     return(pnorm((delta*sqrt(n))/sd +v))
#   else
#     # power function (3.1) in Matthews (p. 32) (in Spiegelhalter notation)
#     return(pnorm((delta*sqrt(n))/sd +v) + 1 - pnorm((delta*sqrt(n))/sd -v))
# }
# 
## Replace functions with ReplicationSuccess::powerSignificance
# ClassicalPower <- function(n0, n, delta, sd, sig.level=0.05, type="two.sample"){
#   stopifnot(type=="two.sample" | type=="one.sample")
#   if(type=="two.sample")
#     sd <- sqrt(2)*sd
#   term1 <- sqrt(n0/(n0+n))
#   v <- qnorm(p = sig.level/2)
#   term2 <- delta*sqrt(n)/sd+v
#   return(pnorm(term1*term2))
# }


########################
## Plots
########################
pval_label <- expression(paste(italic(p), "-Wert der Originalstudie (zweiseitig)"))

powerPlot <- function(p.min = 0.00001, p.max = 0.05){
  ## Function to draw box for power curves depending on selected max p-value
  p <- exp(seq(from = log(p.min), to = log(p.max), length.out = 250))
  par(las = 1)
  plot(0, 0, col = "white", xlim = c(0, p.max), ylim = c(0,100), ylab = "Power (in %)", 
       xlab = pval_label, axes = FALSE, cex.lab = 1.2) #two-sided p-value in original study
  
  # x-axis
  if(p.max <= 0.03)
    axis(1, at = seq(from = 0, to = 1, by = 0.005), 
         labels = as.character(seq(from = 0, to = 1, by = 0.005)), cex.axis = 1.2)
  if(p.max > 0.03 && p.max <= 0.1)
    axis(1, at = seq(from = 0, to = 1, by = 0.01), 
         labels = as.character(seq(from = 0, to = 1, by = 0.01)), cex.axis = 1.2)
  if(p.max > 0.1)
    axis(1, at = seq(from = 0, to = 1, by = 0.02), 
         labels = as.character(seq(from = 0, to = 1, by = 0.02)), cex.axis = 1.2)
  
  # y-axis
  axis(2, at = seq(from = 0, to = 100, by = 10), 
       labels = as.character(seq(from = 0, to = 100, by = 10)), cex.axis = 1.2)
  box()
}

## Replace with RepicationSuccess::powerSignificance 
# Power <- function(p, c, alpha, prior="point"){
#   n0 <- 50
#   n <- c*50
#   myt <- p2t(p)
#   if(prior=="point"){
#     power <- StandardPower(n=n, delta=myt/sqrt(n0), sd=1, sig.level=alpha, type="one.sample")
#   }
#   if(prior=="normal"){
#     power <- ClassicalPower(n0=n0, n=n, delta=myt/sqrt(n0), sd=1, sig.level=alpha, type="one.sample")
#   }
#   return(power)
# }

segm.power <- function(p, c, alpha, prior = "conditional", col = 2){
  ## Function to draw segment lines on existing plot
  power <- 100*powerSignificance(po = p, c = c, level = alpha, designPrior = prior)
  segments(x0 = p, y0 = 0, x1 = p, y1 = power, col = col, lty = 2, lwd = 1.5)
  segments(x0 = 0, y0 = power, x1 = p, y1 = power, col = col, lty = 2, lwd = 1.5)
}

## Replace with one function using ReplicationSuccess::powerSignificance
# standardPowerCurve <- function(p, c, alpha, col=3){
#   n0 <- 50
#   n <- c*n0
#   sd <- 1
#   myt <- p2t(p)
#   mydelta <- myt*sd/sqrt(n0)
#   standardPower <- numeric()
#   for(i in 1:length(mydelta)){
#     standardPower[i] <- StandardPower(n=n, delta=mydelta[i], sd=1, sig.level=alpha, type="one.sample")
#   }
#   lines(p, 100*standardPower, col=4, lwd=2)
# }
# 
# classicalPowerCurve <- function(p, c, alpha, col=5){
#   n0 <- 50
#   n <- c*n0
#   sd <- 1
#   myt <- p2t(p)
#   mydelta <- myt*sd/sqrt(n0)
#   classicalPower <- numeric()
#   for(i in 1:length(mydelta)){
#     classicalPower[i] <- ClassicalPower(n0=n0, n=n, delta=mydelta[i], sd=1, sig.level=alpha, type="one.sample")
#   }
#   lines(p, 100*classicalPower, col=3, lwd=2)
# }

PowerCurve <- function(p, c, alpha, prior = "conditional", ...){
  ## Function to draw power curve on existing plot
  power <- powerSignificance(po = p, c = c, level = alpha, designPrior = prior)
  lines(p, 100*power, ...)
}
