# ======================================================================
# Solutions to exercises from tutorial on R-package ReplicationSuccess #
# C. Micheloud, S.Pawel 
# charlotte.micheloud@uzh.ch
# ======================================================================


# Loading package and data
# ----------------------------------------------------------------------
library("ReplicationSuccess")
data("RProjects")


# Some data manipulations
# ----------------------------------------------------------------------
## computing z-values and variance ratio c
RProjects$zo <- RProjects$fiso/RProjects$se_fiso
RProjects$zr <- RProjects$fisr/RProjects$se_fisr
RProjects$c <- RProjects$se_fiso^2/RProjects$se_fisr^2

## computing 1-sided p-values (all original estimates coded in positive direction)
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


# Exercise 1.5 (calculating d and dmin for discrepant studies)
# ----------------------------------------------------------------------

discrepantDF$d <- discrepantDF$fisr/discrepantDF$fiso
discrepantDF$dminSign <- effectSizeSignificance(zo = discrepantDF$zo,
                                                c = discrepantDF$se_fiso^2/discrepantDF$se_fisr^2, level = 0.025)
discrepantDF$dminRS <- effectSizeReplicationSuccess(zo = discrepantDF$zo,
                                                    c = discrepantDF$se_fiso^2/discrepantDF$se_fisr^2, level = 0.025,
                                                    type = "golden")

discrepantDF[, c("d", "dminSign", "dminRS")]

# Exercise 2.1
# ----------------------------------------------------------------------
par(mfrow = c(1, 1), las = 1)

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
                                    shrinkage = 0.25)
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
                                        shrinkage = 0.25)
plot(po.plot, ss.cond, type = "l", ylim = c(0, 4), col = "red", 
     xlab = "One-sided original p-value", ylab = "Relative sample size")
lines(po.plot, ss.condshrink, col = "red", lty = 2)
legend("topleft", 
       c("Conditional without shrinkage", "Conditional with 25% shrinkage"), 
       col = c("red", "red"), lty = c(1, 2), bty = "n")



# Exercise 3.1
# ----------------------------------------------------------------------
pow.cond2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "conditional", 
                                     level  = 0.025, alternative = "one.sided")

pow.pred2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "predictive", 
                                     level = 0.025, alternative = "one.sided")

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


# Exercise 4.1
# ----------------------------------------------------------------------


for(i in unique(protzko2020$experiment)){
  dat = subset(protzko2020, experiment == i)
  ma <- metagen(TE = smd, seTE = se, studlab = lab,
                exclude = (type == "original" | type == "self-replication"),
                data = dat, sm = "SMD")
  print(ma)
  forest(ma)
}

# Exercise 4.2
# ----------------------------------------------------------------------

par(las = 1, mar = c(10, 5, 4, 2) + 0.1, 
    mfrow = c(1,1))
ma_mat <-   matrix(data = NA,
                   nrow = length(unique(protzko2020$experiment)),
                   ncol = 4)

colnames(ma_mat) <-  c("I2", "zr", "sefixed", "experiment")

for(i in 1:length(unique(protzko2020$experiment))){
  ex <- unique(protzko2020$experiment)[i]
  dat = protzko2020[protzko2020$experiment == ex, ]
  ma <- metagen(TE = smd, seTE = se, data = dat, sm = "SMD",
                exclude = (type == "original" | type == "self-replication"))
  ma_mat[i, 1] <-  ma$I2*100
  ma_mat[i, 2] <- ma$zval.fixed
  ma_mat[i, 3] <- ma$seTE.fixed
  ma_mat[i, 4] <-  ex
  
}

matplot(ma_mat[, c(1, 4)],
        pch = 19,
        xaxt = "n",
        xlab = "",
        ylab = expression(paste("Higgins ", I^2, " (%)")),
        cex.axis = 1.5,
        cex.lab = 1.5,
        ylim = c(0, 100),
        cex = 2.5,
        col = alpha("blue", 0.6))

axis(side = 1, at = 1:16, labels = unique(protzko2020$experiment),
     las = 3,
     cex.axis = 1.5)



# Exercise 4.3
# ----------------------------------------------------------------------


zr_ma <- as.numeric(ma_mat[, 2])
pr_ma <- z2p(zr_ma, alternative = "one.sided")
ori_subset <- subset(protzko2020, type == "original")
zo_ma <- ori_subset$smd/ori_subset$se
po_ma <- z2p(zo_ma, alternative = "one.sided")
c_ma = ori_subset$se^2/as.numeric(ma_mat[, 3])^2
pS_ma <- pSceptical(zo = zo_ma, zr = zr_ma, c = c_ma )


matrix(c(unique(protzko2020$experiment),
                        po_ma,
                        pr_ma,
                        pS_ma),
                      ncol = 4)
