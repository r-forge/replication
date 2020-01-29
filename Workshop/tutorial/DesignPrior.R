
sampleSizeZtest = function(delta, sd, sig.level = 0.05, power){
  u <- qnorm(p = power)
  v <- qnorm(p = 1 - sig.level/2)
  n <- 2*(u + v)^2*sd^2/delta^2
  return(n)
}

grid <- seq(0.075, 0.225, 0.001)
deltahat <- 0.15
se <- 0.025
ss <- sampleSizeZtest(delta = grid, sd = 0.4, sig.level = 0.01, power = 0.95)
ss1 <- sampleSizeZtest(delta = deltahat, sd = 0.4, sig.level = 0.01, power = 0.95)

set.seed(12345)
nsim <- 100000
deltahat2 <- rnorm(nsim, mean=deltahat, sd=se)
ss2 <- sampleSizeZtest(delta = deltahat2, sd = 0.4, sig.level = 0.01, power = 0.95)

library(gplots)
par(las=1)

plot(grid, ss, type="l", xlab="effect", ylab="sample size", ylim=c(0, max(ss)))
lines(rep(deltahat, 2), c(0, ss1), lty=2)
points(deltahat, ss1, pch=19, cex=1, col=2) 
points(deltahat, 0, pch=19, cex=1, col=1) 
lines(c(0, deltahat), rep(ss1, 2), lty=2, col=2)
axis(2, at=ss1, label=round(ss1), col=2, col.axis=2)
nshow <- 20
for(i in 1:nshow){
    lines(rep(deltahat2[i], 2), c(0, ss2[i]), lty=2, col="grey")
    points(deltahat2[i], ss2[i], pch=19, cex=1, col=2) 
    points(deltahat2[i], 0, pch=19, cex=1, col=1) 
    lines(c(0, deltahat2[i]), rep(ss2[i], 2), lty=2, col="pink")
}
plotCI(x = deltahat, y=0, li=quantile(deltahat2, 0.025), ui = quantile(deltahat2, 0.975),
       col=1, lwd=2, add=TRUE, err="x")

plotCI(x = min(grid), y=mean(ss2), li=quantile(ss2, 0.025), ui = quantile(ss2, 0.975),
       col=2, lwd=2, add=TRUE, err="y")
lines(rep(quantile(deltahat2, 0.025), 2), c(0, quantile(ss2, 0.975)), lty=2)
lines(rep(quantile(deltahat2, 0.975), 2), c(0, quantile(ss2, 0.025)), lty=2)
lines(c(0, quantile(deltahat2, 0.975)), rep(quantile(ss2, 0.025), 2), lty=2, col=2)
lines(c(0, quantile(deltahat2, 0.025)), rep(quantile(ss2, 0.975), 2), lty=2, col=2)
axis(2, at=mean(ss2), label=round(mean(ss2)), col=2, col.axis=2)

