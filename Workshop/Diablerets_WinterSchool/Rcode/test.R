
source("plotLevel.R")

par(mfrow=c(2,1))
plotLevel(method="RS")
plotLevel(method="BFs")


plotLevel(method="RS", c=1)
plotLevel(method="BFs", c=1)

plotLevel(method="signif", c=1)
plotLevel(method="BFs", c=1)
plotLevel(method="BFs", c=1)

plotLevel(method="signif", c=100, level=0.05)
plotLevel(method="BFs", c=100, level=0.05)

par(mfrow=c(2,2))
plotLevel(method="signif", c=1, level=0.025)
plotLevel(method="RS", c=1, level=0.025, type="nominal")
plotLevel(method="BFs", c=1, level=0.025)
