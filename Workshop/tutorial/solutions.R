# ======================================================================
# Solutions to exercises from tutorial on R-package ReplicationSuccess #
# ======================================================================


# Loading package and data
# ----------------------------------------------------------------------
library("ReplicationSuccess")
data("RProjects")


# Some data manipulations
# ----------------------------------------------------------------------
# computing z-values and variance ratio c
RProjects$zo <- RProjects$fiso/RProjects$se_fiso
RProjects$zr <- RProjects$fisr/RProjects$se_fisr
RProjects$c <- RProjects$se_fiso^2/RProjects$se_fisr^2

# computing 1-sided p-values (all original estimates were positive)
RProjects$po1 <- z2p(RProjects$zo, alternative = "greater")
RProjects$pr1 <- z2p(RProjects$zr, alternative = "greater")


# Exercise 1.1 (Success based on statistical significance)
# ----------------------------------------------------------------------
## compute for all projects
significant_O <- RProjects$po1 < 0.025
significant_R <- RProjects$pr1 < 0.025
RProjects$successSignif <- significant_O & significant_R
allDF <- data.frame(project = "all", 
                    success = mean(RProjects$successSignif)*100)

## compute for each project
successSignif <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             success = mean(data_project$successSignif)*100)
})
(signSuccessDF <- rbind(do.call("rbind", successSignif), allDF))

## plots
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$successSignif
  col_success <- color <- ifelse(success == FALSE, "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(mean(success)*100, 0), 
                  "% (", sum(success), "/", length(success), ")")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_success, pch = 20)
  legend("topleft", 
         legend = c(expression(paste("both ", italic(p), "-values < 0.025")), 
                    expression(paste("not both ", italic(p), "-values < 0.025"))),
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


# Exercise 1.2 (Detecting incompatibility with Q-test)
# ----------------------------------------------------------------------
## compute for all projects
RProjects$pQ <- Qtest(thetao = RProjects$fiso, thetar = RProjects$fisr,
                      seo = RProjects$se_fiso, ser = RProjects$se_fisr)
RProjects$Qincompatible <- RProjects$pQ <= 0.05
allQDF <- data.frame(project = "all", 
                     incomp = mean(RProjects$Qincompatible)*100)

## compute for each project
Qproject <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             incomp = mean(data_project$Qincompatible)*100)
})
(QDF <- rbind(do.call("rbind", Qproject), allQDF))

## plot
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  incompatible <- data_project$pQ < 0.05
  PropIncomp <- mean(incompatible)
  color <- ifelse(incompatible == FALSE, "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(PropIncomp*100, 0), '% incompatible')
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = color, pch = 20)
  legend("topleft", 
         c(expression(italic(p)[Q] < 0.05), expression(italic(p)[Q] >= 0.05)),
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}

# Exercise 1.3 (sceptical p-value)
# ----------------------------------------------------------------------
## compute for all projects
RProjects$ps <- with(RProjects, pSceptical(zo = zo, zr = zr, c = c, 
                                           alternative = "one.sided",
                                           type = "golden"))
RProjects$pSsuccess <- RProjects$ps < 0.025
allpsDF <- data.frame(project = "all", 
                      success = mean(RProjects$pSsuccess)*100)

## compute for each project
psProjects <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             success = mean(data_project$pSsuccess)*100)
})
(psDF <- rbind(do.call("rbind", psProjects), allpsDF))


## plot
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$pSsuccess
  col_success <- ifelse(success == FALSE, "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(mean(success)*100, 0), "% (",
                  sum(success), "/",  length(success), ")")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_success, pch = 20)
  legend("topleft", c(expression(italic(p)[s] < 0.025),
                      expression(italic(p)[s] >= 0.025)),
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


# Exercise 1.4 (Looking closer at studies where discrepancies between methods)
# ----------------------------------------------------------------------
## plot
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  discrep <- data_project$pSsuccess != data_project$successSignif
  col_discord <- ifelse(discrep == TRUE,
                        ifelse(data_project$pSsuccess == TRUE,
                               "#8B0000B3", "#00008AB3"), "#B2B2B299")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = p, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_discord, pch = 20)
  legend("topleft", legend = c(expression(paste(italic(p)[s], " only")),
                               "significance only"),
         pch = 20, col = c("#8B0000B3", "#00008AB3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}

## looking closer at sample size (ratio), effect estimates, p-values
discrepantDF <- subset(RProjects, pSsuccess != successSignif)
discrepantDF[, c("study", "project", "no", "nr", "c", "ro", "rr", "po1", "pr1", 
                 "ps", "pSsuccess", "successSignif")]


# Exercise 2.1
# ----------------------------------------------------------------------
par(las = 1)

# first plot
pval.or <- c(0.0001, 0.001, 0.005, 0.01, 0.025)
pow.cond <- powerSignificance(zo = p2z(pval.or, alternative = "one.sided"), 
                              c = 1, designPrior = "conditional")
pow.pred <- powerSignificance(zo = p2z(pval.or, alternative = "one.sided"),
                              c = 1, designPrior = "predictive")
plot(pval.or, pow.cond*100, type = "p", col = "red", ylim = c(0, 100), 
     xlab = "One-sided original p-value", ylab = "Power (in %)", pch = 20)
points(pval.or, pow.pred*100, col = "blue", pch = 20)
legend("topright", c("Conditional", "Predictive"), col = c("red", "blue"), 
       pch = 20, bty = "n")

# second plot
po.plot <- seq(0.000001, 0.025, by = 0.0001)
pow.cond1 <- powerSignificance(zo = p2z(po.plot, alternative = "one.sided"),
                               c = 1, designPrior = "conditional")
pow.pred1 <- powerSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                               c = 1, designPrior = "predictive")
plot(po.plot, pow.cond1*100, col = "red", ylim = c(0, 100), type = "l", 
     xlab = "One-sided original p-value", ylab = "Power (in %)")
lines(po.plot, pow.pred1*100, col = "blue")
abline(h = 50, lty = 3)
axis(2, at = 50, label = "50", col = "gray40")
legend("topright",  c("Conditional", "Predictive"),  col = c("red", "blue"), 
       lty = 1, bty = "n")

# third plot
pow.cond <- powerSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                              c = 1, designPrior = "conditional")
pow.condshrink <- powerSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                                    c = 1, designPrior = "conditional", 
                                    shrinkage = 0.75)
plot(po.plot, pow.cond*100, col = "red", ylim = c(0, 100), type = "l", 
     xlab = "One-sided original p-value", ylab = "Power (in %)")
lines(po.plot, pow.condshrink*100, col = "red", lty = 2)
abline(h = 50, lty = 3)
axis(2, at = 50, label = "50", col = "gray40")
legend("topright", 
       c("Conditional without shrinkage", "Conditional with 25% shrinkage"), 
       col = c("red", "red"), lty = c(1, 2), bty = "n")


# Exercise 2.2
# ----------------------------------------------------------------------
# first plot
ss.cond <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                                  power = 0.8, designPrior = "conditional")
ss.pred <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"),
                                  power = 0.8, designPrior = "predictive" )
plot(po.plot, ss.cond, type = "l", ylim = c(0, 4), col = "red", 
     xlab = "One-sided original p-value", ylab = "Relative sample size")
lines(po.plot, ss.pred, col = "blue")
legend("topleft", c("Conditional", "Predictive"), col = c("red", "blue"), 
       lty = 1, bty = "n")

# second plot
ss.cond <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"),
                                  power = 0.8, designPrior = "conditional")
ss.condshrink <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                                        power = 0.8, designPrior = "conditional", 
                                        shrinkage = 0.75)
plot(po.plot, ss.cond, type = "l", ylim = c(0, 4), col = "red", 
     xlab = "One-sided original p-value", ylab = "Relative sample size")
lines(po.plot, ss.condshrink, col = "red", lty = 2)
legend("topleft", 
       c("Conditional without shrinkage", "Conditional with 25% shrinkage"), 
       col = c("red", "red"), lty = c(1, 2), bty = "n")


# Exercise 2.3
# ----------------------------------------------------------------------
library("lattice")
library("ggplot2")
eco <- subset(RProjects, project == "Experimental Economics")
pow_c1 <- powerSignificance(zo = p2z(eco$po), c = eco$c, level = 0.025, 
                            alternative = "one.sided", designPrior = "conditional")

pow_p1 <- powerSignificance(zo = p2z(eco$po), c = eco$c, level = 0.025, 
                            alternative = "one.sided", designPrior = "predictive")

mat1 <- matrix(c(rep(eco$study, times = 2), pow_c1*100, pow_p1*100,
                 rep(c("Conditional","Predictive"), times = c(nrow(eco), nrow(eco)))), 
               ncol = 3)

colnames(mat1) <- c("ID", "power", "group")
mat1 <- as.data.frame(mat1)
mat1$ID <- as.factor(mat1$ID)
mat1$group <- factor(mat1$group, levels = c("Conditional", "Predictive"), 
                     order = TRUE)
mat1$power <- as.numeric(as.character(mat1$power))

trellis.par.set(
  list(
    plot.symbol = list(col = 1, pch = 20, cex = 0.7),
    box.rectangle = list(col = alpha("black", 0.5)),
    box.umbrella = list(lty = 1, col = alpha("black", 0.5)),
    strip.background = list(col = "white")
  )
)
panel_bw <- function(x, y, groups, subscripts, ...) {
  panel.bwplot(x = x, y = y, ...)
  tapply(1:length(y), groups[subscripts], function(i) {
    llines(
      x = 1:nlevels(x),
      y = y[i][order(x[i])],
      col = rgb(.2, .2, .2, .2)
    )
    lpoints(
      x = 1:nlevels(x),
      y = y[i][order(x[i])],
      col = rgb(.2, .2, .2, .2),
      pch = 16,
      cex = 1.2
    )
  })
  panel.abline(h = 50, col="red", lty=2)
}
bwplot(power~group, data = mat1, groups = ID, panel = panel_bw, 
       xlab = list(""), ylab = list("Power (in %)"),
       between = list(x = 1), scales = list(x = "free", y = "free", rot = 0),
       pch = "|")


# Exercise 2.4
# ----------------------------------------------------------------------
philo <- subset(RProjects, project == "Experimental Philosophy")
ss_cond1 <- sampleSizeSignificance(zo = p2z(philo$po), power = 0.95,
                                   designPrior = "conditional")
plot(philo$c, ss_cond1,  xlim = c(0.01, 50), ylim = c(0.01, 50), 
     xlab = "Reported relative sample size", 
     ylab = "Calculated relative sample size", 
     pch = 20, log = "xy", xaxt = "n", yaxt = "n")
abline(a = 0, b = 1, col = "lightpink")
axis(1, las = 1, at = c(0.01, 0.1, 1, 10, 50), 
     labels = c("1/100", "1/10", "1", "10", "50"))
axis(2, las = 1, at = c(0.01, 0.1, 1, 10, 50), 
     labels = c("1/100", "1/10", "1", "10", "50"))


# Exercise 3.1
# ----------------------------------------------------------------------
pow.cond2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "conditional", 
                                     level  = 0.065, alternative = "one.sided")

pow.pred2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "predictive", 
                                     level = 0.065, alternative = "one.sided")

plot(po.plot, pow.cond2*100, type = "l", col = "red", ylim = c(0, 100),
     xlab = "One-sided original p-value", ylab = "Power (in %)", lty = 1)
lines(po.plot, pow.pred2*100, col = "blue")
lines(po.plot, pow.cond1*100, col = "red", lty = 2)
lines(po.plot, pow.pred1*100, col = "blue", lty = 2)
abline(h = 50, lty = 3)
legend("topright", c("Conditional", "Predictive"), col = c("red", "blue"), 
       lty = 1, bty = "n")
legend("bottomleft", c("Replication success", "Significance"), 
       lty = c(1, 2), bty = "n")


# Exercise 3.2
# ----------------------------------------------------------------------
pow_c2 <- powerReplicationSuccess(zo = p2z(eco$po), c = eco$c, level = 0.065, 
                                  alternative = "one.sided", designPrior = "conditional")
pow_p2 <- powerReplicationSuccess(zo = p2z(eco$po), c = eco$c, level = 0.065, 
                                  alternative = "one.sided", designPrior = "predictive")
mat2 <- matrix(c(rep(eco$study, times = 2),
                 pow_c2*100, pow_p2*100, 
                 rep(c("Conditional", "Predictive"), times=c(nrow(eco),nrow(eco)))),
               ncol = 3)
colnames(mat2) <- c("ID", "power", "group")
mat2 <- as.data.frame(mat2)
mat2$ID <- as.factor(mat2$ID)
mat2$group <- factor(mat2$group, levels = c("Conditional", "Predictive"), order = TRUE)
mat2$power <- as.numeric(as.character(mat2$power))
bwplot(power~group, data = mat2, groups = ID, panel = panel_bw, 
       xlab = list(""), ylab = list("Power (in %)"),
       between = list(x = 1), scales = list(x = "free", y = "free", rot = 0),
       pch = "|")


# Exercise 3.3
# ----------------------------------------------------------------------
ss_c2 <- sampleSizeReplicationSuccess(zo = p2z(philo$po), power = 0.95, 
                                      level = 0.065, 
                                      alternative = "one.sided",
                                      designPrior = "conditional")
plot(philo$c, ss_c2, xlim = c(0.01, 50), ylim = c(0.01, 50), 
     xlab = "Reported relative sample size", 
     ylab = "Calculated relative sample size", 
     pch = 20, log = "xy", xaxt = "n", yaxt = "n")
abline(a = 0, b = 1, col = "lightpink")
axis(1, las = 1, at = c(0.01, 0.1, 1, 10, 50), 
     labels = c("1/100", "1/10", "1", "10", "50"))
axis(2, las = 1, at = c(0.01, 0.1, 1, 10, 50), 
     labels = c("1/100", "1/10", "1", "10", "50"))