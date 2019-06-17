
p.to.z <- function(p)
  return(qnorm(1-p/2))

StandardPower <- function(n, delta, sd, sig.level=0.05, type="two.sample", strict=FALSE){
  stopifnot(type=="two.sample" | type=="one.sample")
  if(type=="two.sample")
    sd <- sqrt(2)*sd
  v <- qnorm(p = sig.level/2)
  if(strict==FALSE)
    return(pnorm((delta*sqrt(n))/sd +v))
  else
    # power function (3.1) in Matthews (p. 32) (in Spiegelhalter notation)
    return(pnorm((delta*sqrt(n))/sd +v) + 1 - pnorm((delta*sqrt(n))/sd -v))
}

ClassicalPower <- function(n0, n, delta, sd, sig.level=0.05, type="two.sample"){
  stopifnot(type=="two.sample" | type=="one.sample")
  if(type=="two.sample")
    sd <- sqrt(2)*sd
  term1 <- sqrt(n0/(n0+n))
  v <- qnorm(p = sig.level/2)
  term2 <- delta*sqrt(n)/sd+v
  return(pnorm(term1*term2))
}


########################
## Plots
########################

powerPlot <- function(p.min=0.00001, p.max=0.05){
  p <- exp(seq(log(p.min), log(p.max), length.out=250))
  par(las=1)
  # par(las=1, pty="s")
  plot(0,0, col="white", xlim=c(0, p.max), ylim=c(0,100), ylab="Power (in %)", xlab="two-sided p-value in original study", axes=FALSE, cex.lab=1.2)
       # main="Significance")
  if(p.max <= 0.03)
    axis(1, at=seq(from=0, to=1, by=0.005), labels=as.character(seq(from=0, to=1, by=0.005)), cex.axis=1.2)
  if(p.max > 0.03 && p.max <= 0.1)
    axis(1, at=seq(from=0, to=1, by=0.01), labels=as.character(seq(from=0, to=1, by=0.01)), cex.axis=1.2)
  if(p.max > 0.1)
    axis(1, at=seq(from=0, to=1, by=0.02), labels=as.character(seq(from=0, to=1, by=0.02)), cex.axis=1.2)
  # if(p.max < 0.3)
  #   axis(1, at=seq(from=0, to=0.3, by=0.005), labels=as.character(seq(from=0, to=0.3, by=0.005)), cex.axis=1.2)
  # if(p.max >= 0.3 & p.max <= 0.1)
  #   axis(1, at=seq(from=0, to=1, by=0.01), labels=as.character(seq(from=0, to=1, by=0.01)), cex.axis=1.2)
  # if(p.max > 0.1)
  #   axis(1, at=seq(from=0, to=1, by=0.02), labels=as.character(seq(from=0, to=1, by=0.02)), cex.axis=1.2)
  # if(p.max > 0.2)
  # axis(1, at=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), label=as.character(c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05)))
  axis(2, at=seq(from=0, to=100, by=10), labels=as.character(seq(from=0, to=100, by=10)), cex.axis=1.2)
  box()
}

Power <- function(p, c, alpha, prior="point"){
  n0 <- 50
  n <- c*50
  myt <- p.to.z(p)
  if(prior=="point"){
    power <- StandardPower(n=n, delta=myt/sqrt(n0), sd=1, sig.level=alpha, type="one.sample")
  }
  if(prior=="normal"){
    power <- ClassicalPower(n0=n0, n=n, delta=myt/sqrt(n0), sd=1, sig.level=alpha, type="one.sample")
  }
  return(power)
}

segm.power <- function(p, c, alpha, prior="point", col=2){
  power <- 100*Power(p=p, c=c, alpha=alpha, prior=prior)
  segments(x0=p, y0=0, x1=p, y1=power, col=col, lty=2, lwd=1.5)
  segments(x0=0, y0=power, x1=p, y1=power, col=col, lty=2, lwd=1.5)
}

standardPowerCurve <- function(p, c, alpha, col=3){
  n0 <- 50
  n <- c*n0
  sd <- 1
  myt <- p.to.z(p)
  mydelta <- myt*sd/sqrt(n0)
  standardPower <- numeric()
  for(i in 1:length(mydelta)){
    standardPower[i] <- StandardPower(n=n, delta=mydelta[i], sd=1, sig.level=alpha, type="one.sample")
  }
  lines(p, 100*standardPower, col=col, lwd=2)
}

classicalPowerCurve <- function(p, c, alpha, col=5){
  n0 <- 50
  n <- c*n0
  sd <- 1
  myt <- p.to.z(p)
  mydelta <- myt*sd/sqrt(n0)
  classicalPower <- numeric()
  for(i in 1:length(mydelta)){
    classicalPower[i] <- ClassicalPower(n0=n0, n=n, delta=mydelta[i], sd=1, sig.level=alpha, type="one.sample")
  }
  lines(p, 100*classicalPower, col=3, lwd=2)
}
