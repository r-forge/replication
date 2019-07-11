## Replace with ReplicationSuccess::powerSignificace
# pSignificant <- function(p=NULL, t=p2t(p), c=1, sig.level=0.05){
#   pSig <- pnorm(p2t(sig.level), mean=sqrt(c)*t, sd=sqrt(1+c), lower.tail=FALSE)
#   return(pSig)
# }

## Replace with ReplicationSuccess::sampleSizeSignificance
## classical sample size calculation
# sampleSizeZtest <- function(delta, sd, sig.level = 0.05, power, type="two.sample"){
#   stopifnot(type=="two.sample" | type=="one.sample")
#   if(type=="two.sample")
#     sd <- sqrt(2)*sd
#   u <- qnorm(p = power)
#   v <- qnorm(p = 1-sig.level/2)
#   n <- (u+v)^2*sd^2/delta^2
#   return(n)
# }

########################
## Plots
########################
sampleSizePlot <- function(p.min = 0.00001, p.max = 0.05, c.max){
  ## Function to draw box for sample size curves depending on selected max p-value and c
  p <- exp(seq(from = log(p.min), to = log(p.max), length.out = 250))
  par(las = 1)
  plot(0, 0, col = "white", xlim = c(0, p.max), ylim = c(0, c.max), 
       ylab = "Relative Stichprobengrösse", #Relative sample size
       xlab = expression(paste(italic(p), "-Wert der Originalstudie (zweiseitig)")), 
       axes = FALSE, cex.lab = 1.2) #two-sided p-value in original study
  
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
  if(c.max < 2)
    axis(2, at = seq(from = 0, to = 2, by = 0.2), 
         labels = as.character(seq(from = 0, to = 2, by = 0.2)), cex.axis = 1.2)
  if(c.max >= 2 && c.max < 6)
    axis(2, at = seq(from = 0, to = 6, by = 0.5), 
         labels = as.character(seq(from = 0, to = 6, by = 0.5)), cex.axis = 1.2)
  if(c.max >= 6 && c.max < 16)
     axis(2, at = seq(from = 0, to = c.max + 1, by = 1), 
          labels = as.character(seq(from = 0, to = c.max + 1, by = 1)), cex.axis = 1.2)
  if(c.max >= 16 && c.max < 30)
    axis(2, at = seq(from = 0, to = c.max + 1, by = 2), 
         labels = as.character(seq(from = 0, to = c.max + 1, by = 2)), cex.axis = 1.2)
  if(c.max >= 30 && c.max < 60)
    axis(2, at = seq(from = 0, to = 200, by = 5),
         labels = as.character(seq(from = 0, to = 200, by = 5)), cex.axis = 1.2)
  if(c.max >= 60)
    axis(2, at = seq(from = 0, to = 200, by = 10),
         labels = as.character(seq(from = 0, to = 200, by = 10)), cex.axis = 1.2)
  box()
}

## Replace with ReplicationSuccess::sampleSizeSignificance
# classicalTarget <- function(c, p, alpha=0.05, power=0.8){
#   term <- pSignificant(p=p, c=c, sig.level=alpha)
#   return(term - power)
# }

## Replace with ReplicationSuccess::sampleSizeSignificance
# relSampleSize <- function(p, alpha = 0.05, power = 0.8, prior = "point"){
#   n0 <- 50
#   sd <- 1
#   myt <- p2t(p)
#   mydelta <- myt*sd/sqrt(n0)
#   if(prior=="point"){
#     SampleSize <- sampleSizeZtest(delta=mydelta, sd=sd, sig.level=alpha, 
#                                   power=power, type="one.sample")
#     c <- SampleSize/n0
#   }
#   if(prior=="normal"){
#     
#     c.l <- 0.001
#     c.u <- 1000
#     # c.l <- 0.01
#     # c.u <- 150
#     target.l <- classicalTarget(c=c.l, p=p, alpha=alpha, power=power)
#     target.u <- classicalTarget(c=c.u, p=p, alpha=alpha, power=power)
#     if(sign(target.l)==sign(target.u))
#       c <- NA
#     else
#       c <- uniroot(classicalTarget, lower=c.l, upper=c.u, p=p, alpha=alpha, power=power)$root
#     # return(c)
#     
#     # c <- uniroot(classicalTarget, lower=0.001, upper=10^6, p=p, alpha=alpha, power=power)$root
#     # SampleSize <- ClassicalSampleSize(n0=n0, delta=mydelta, sd=sd, sig.level=alpha, 
#     #                                   power=power, type="one.sample")
#   }
#   return(c)
# }


segm.rss <- function(p, alpha = 0.05, power = 0.8, prior = "conditional", col = 2){
  ## Function to draw segment lines on existing plot
  c <- sampleSizeSignificance(po = p, power = power, level = alpha, designPrior = prior)
  segments(x0 = p, y0 = 0, x1 = p, y1 = c, col = col, lty = 2, lwd = 1.5)
  segments(x0 = 0, y0 = c, x1 = p, y1 = c, col = col, lty = 2, lwd = 1.5)
}

# SamSizeCurve <- function(p, alpha=0.05, power=0.8, prior="point", col){
#   n0 <- 50
#   sd <- 1
#   myt <- p2t(p)
#   mydelta <- myt*sd/sqrt(n0)
#   # myt <- p2t(p)
#   # n.o <- ((myt*sd)/delta)^2
#   c <- numeric()
#   # if(prior=="point"){
#       for(i in 1:length(mydelta)){
#         c[i] <- relSampleSize(p=p[i], alpha=alpha, power=power, prior=prior)
#       }
#   lines(p, c, col=col, lwd=2)
# }

SamSizeCurve <- function(p, alpha = 0.05, power = 0.8, prior = "conditional", ...){
  ## Function to draw relative sample size curve on existing plot
  c <- sampleSizeSignificance(po = p, power = power, level = alpha, designPrior = prior)
  lines(p, c, ...)
}


