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


# Exercise 1.1
# ----------------------------------------------------------------------
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  significantO <- data_project$po1 < 0.025
  significantR <- data_project$pr1 < 0.025
  success <- significantO & significantR 
  
  cat(paste0(p, ": \n"))
  cat(paste0(round(mean(significantO)*100, 1), 
             "% original studies significant (", 
             sum(significantO), "/", length(significantO), ")\n"))
  cat(paste0(round(mean(significantR)*100, 1), 
             "% replications significant (", 
             sum(significantR), "/", length(significantR), ")\n"))
  cat(paste0(round(mean(success)*100, 1), 
             "% both significant, same direction (",
             sum(success), "/", length(success), ") \n \n"))
}

# plots
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$po1 < 0.025 & data_project$pr1 < 0.025 
  col_success <- color <- ifelse(success == FALSE, "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(mean(success)*100, 0), 
                  "% (", sum(success), "/", length(success), ")")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_success, pch = 20)
  legend("topleft", c("both significant", "not both significant"), 
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


# Exercise 1.2
# ----------------------------------------------------------------------
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  PI <- predictionInterval(thetao = data_project$fiso, 
                           seo = data_project$se_fiso, 
                           ser = data_project$se_fisr)
  
  PI <- tanh(PI) # transforming back to correlation scale
  
  within <- (data_project$rr < PI$upper) & (data_project$rr > PI$lower)
  coverage <- mean(within)
  color <- ifelse(within == TRUE, "#333333B3", "#8B0000B3")
  studynr <- seq(1, nrow(data_project))
  plot(data_project$rr, studynr, col = color, pch = 20, 
       xlim = c(-0.5, 1), xlab = expression(italic(r)[r]), 
       main = paste0(p, ": ", round(coverage*100, 0), "% coverage"), 
       yaxt = "n", ylab = "")
  arrows(PI$lower, studynr, PI$upper, studynr, length = 0.02, 
         angle = 90, code = 3, col = color)
  abline(v = 0, lty = 3)
}


# Exercise 1.3
# ----------------------------------------------------------------------
RProjects$ps <- with(RProjects, pSceptical(zo = zo, zr = zr, c = c, 
                                           alternative = "one.sided"))
threshNominal <- thresholdSceptical(level = 0.025, alternative = "one.sided")
threshControlled <- thresholdSceptical(level = 0.025, type = "controlled", 
                                       alternative = "one.sided")

for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  cat(paste0(p, ": \n"))
  success_scept_nominal <- (data_project$ps < threshNominal)
  success_scept_controlled <- (data_project$ps < threshControlled)
  cat(paste0(round(mean(success_scept_controlled)*100, 1), 
             "% smaller than controlled threshold 0.065 (",
             sum(success_scept_controlled), "/", 
             length(success_scept_controlled), ") \n"))
  cat(paste0(round(mean(success_scept_nominal)*100, 1), 
             "% smaller than nominal threshold 0.025 (",
             sum(success_scept_nominal), "/", 
             length(success_scept_nominal), ") \n"))
  cat("\n")
}

# plots
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success_scept_controlled <- (data_project$ps < threshControlled)
  col_success <- ifelse(success_scept_controlled == FALSE, 
                        "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(mean(success_scept_controlled)*100, 0), 
                  "% (", sum(success_scept_controlled), "/", 
                  length(success_scept_controlled), ")")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_success, pch = 20)
  legend("topleft", c("p-sceptical > 0.065", "p-sceptical < 0.065"), 
         pch = 20, col = c("#333333B3", "#8B0000B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


# Exercise 1.4
# ----------------------------------------------------------------------
options(scipen = 5)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  cat(paste0(p, ": \n"))

  success_scept_controlled <- (data_project$ps < threshControlled)
  success_tradit <- data_project$po1 < 0.025 & data_project$pr1 < 0.025 
  
  if(any(success_tradit != success_scept_controlled) == TRUE){
    discrepant <- which(success_scept_controlled != success_tradit)
    discrepant_df <- data_project[discrepant,
                                  c("ro", "rr", "c", "po1", "pr1", "ps")]
    # print effect estimates, 1-sided p-values, and c of discrepant studies
    print(signif(discrepant_df, 2), row.names = FALSE)
  }
  cat("\n")
}

# plot
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success_scept_controlled <- (data_project$ps < threshControlled)
  success_tradit <- data_project$po1 < 0.025 & data_project$pr1 < 0.025 
  discrepant <- success_scept_controlled != success_tradit
  col_discrepant <- ifelse(discrepant == TRUE, 
                        ifelse(success_scept_controlled == TRUE, 
                               "#8B0000B3", "#00008AB3"), "#B2B2B299")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = p, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_discrepant, pch = 20)
  legend("topleft", c("only p-sceptical successful", 
                      "only significance successful"), 
         pch = 20, col = c("#8B0000B3", "#00008AB3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


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