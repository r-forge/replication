# ======================================================================
# Solutions to exercises from tutorial on R-package ReplicationSuccess #
# C. Micheloud, L. Held
# charlotte.micheloud@uzh.ch
# ======================================================================


# install required packages 
# ----------------------------------------------------------------------

# Some packages
installpkg <- function(pkg) {
  if (!require(pkg, char = TRUE, quietly = TRUE))
    install.packages(pkg); require(pkg, char = TRUE)
}
pkgs <- sapply(c("meta", "devtools", "ggplot2"), installpkg)

# ReplicationSuccess package 
devtools::install_github(repo = "SamCH93/ReplicationSuccess",
                         ref = "T1Econtrol", 
                         build_vignettes = TRUE)
vignette("ReplicationSuccess")



# Loading package and data
# ----------------------------------------------------------------------
library("meta")
library("ggplot2")
library("ReplicationSuccess")
data("RProjects")


# Some data manipulations
# ----------------------------------------------------------------------
## computing z-values and variance ratio c
RProjects$zo <- RProjects$fiso/RProjects$se_fiso
RProjects$zr <- RProjects$fisr/RProjects$se_fisr
RProjects$c <- RProjects$se_fiso^2/RProjects$se_fisr^2



##############################################################################
##########################  EXERCISE SESSION 1  ##############################
##############################################################################


# Exercise 1.1 (Original versus replication effect estimate)
# ----------------------------------------------------------------------

par(mfrow = c(1,1))
with(RProjects, 
     plot(x = fiso, y = fisr, 
          pch = 19, 
          xlim = c(-0.5, 2), 
          ylim = c(-0.5, 2)
          ))
abline(a = 0, b = 1, col = "grey")


# Exercise 1.2 (Success based on statistical significance)
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


# Exercise 1.3 (Detecting incompatibility with Q-test)
# ----------------------------------------------------------------------
 
## PART A

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
         c(expression(italic(p)[Q] <= 0.05), expression(italic(p)[Q] > 0.05)),
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}

## PART B

incompSuccess <- sum(RProjects$pQ <= 0.05 & RProjects$successSignif)
compSuccess <-  sum(RProjects$pQ > 0.05 & RProjects$successSignif)
incompNoSuccess <- sum(RProjects$pQ <= 0.05 & !RProjects$successSignif)
compNoSuccess <- sum(RProjects$pQ > 0.05 & !RProjects$successSignif)

matrix(c(compSuccess, compNoSuccess, incompSuccess, incompNoSuccess), byrow = T,
       nrow = 2)


# Exercise 1.4 (replication success based on the meta-analysis)
# ----------------------------------------------------------------------

## PART A

pval_meta <- c()
for(i in 1:nrow(RProjects)){
  pval_meta[length(pval_meta) + 1] <- with(RProjects[i, ],
                                           metagen(TE = c(fiso, fisr),
                                                   seTE = c(se_fiso, se_fisr))$pval.fixed)
}

RProjects$pMeta <- pval_meta
RProjects$pMsuccess <- RProjects$pMeta < 0.05

allMetaDF <- data.frame(project = "all", 
                     success = mean(RProjects$pMsuccess)*100)

Metaproject <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             success = mean(data_project$pMsuccess )*100)
})
(MetaDF <- rbind(do.call("rbind", Metaproject), allMetaDF))

## plot
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$pMeta < 0.05
  PropSuccess <- mean(success)
  color <- ifelse(success == FALSE, "#333333B3", "#8B0000B3")
  title <- paste0(p, ": ", round(PropSuccess*100, 0), '% success')
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = title, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = color, pch = 20)
  legend("topleft", 
         c(expression(italic(p)[M] <= 0.05), expression(italic(p)[M] > 0.05)),
         pch = 20, col = c("#8B0000B3", "#333333B3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}


## PART B 

ind.example <- which(RProjects$pMsuccess & RProjects$Qincompatible & !RProjects$successSignif)

RProjects[ind.example,]


##############################################################################
##########################  EXERCISE SESSION 2  ##############################
##############################################################################

# Exercise 2.1 (replication success with the sceptical p-value)
# ----------------------------------------------------------------------


## PARTS A and B

## compute for all projects
RProjects$psG <- with(RProjects, pSceptical(zo = zo, zr = zr, c = c, 
                                           alternative = "one.sided",
                                           type = "golden"))
RProjects$psC <- with(RProjects, pSceptical(zo = zo, zr = zr, c = c, 
                                            alternative = "one.sided",
                                            type = "controlled"))
RProjects$pSsuccessG <- RProjects$psG < 0.025
RProjects$pSsuccessC <- RProjects$psC < 0.025



### golden
allpsGDF <- data.frame(project = "all", 
                      success = mean(RProjects$pSsuccessG)*100)

## compute for each project
psGProjects <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             success = mean(data_project$pSsuccessG)*100)
})
(psGDF <- rbind(do.call("rbind", psGProjects), allpsGDF))


## plot
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$pSsuccessG
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


### controlled
allpsCDF <- data.frame(project = "all", 
                       success = mean(RProjects$pSsuccessC)*100)

## compute for each project
psCProjects <- lapply(X = unique(RProjects$project), FUN = function(p) {
  data_project <- subset(RProjects, project == p)
  data.frame(project = p, 
             success = mean(data_project$pSsuccessC)*100)
})
(psGDF <- rbind(do.call("rbind", psCProjects), allpsCDF))


## plot
par(mfrow = c(2, 2), las = 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  success <- data_project$pSsuccessC
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


## PART C 
# Looking closer at studies where discrepancies between methods

par(mfrow = c(1, 1), 
    las = 1)
with(RProjects, 
     plot(psG, psC, 
     pch = 19, 
     xlim = c(0.0001, 1), 
     ylim = c(0.0001, 1), 
     log = "xy"))
abline(a = 0, b = 1, col = "grey")
abline(v = 0.025, col = "red", lty = 2)
abline(h = 0.025, col = "red", lty = 2)

mean(RProjects$pSsuccessG == RProjects$pSsuccessC)


## plot
par(mfrow = c(2, 2))
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  discrep <- data_project$pSsuccessG != data_project$pSsuccessC
  col_discord <- ifelse(discrep == TRUE,
                        ifelse(data_project$pSsuccessG == TRUE,
                               "#8B0000B3", "#00008AB3"), "#B2B2B299")
  plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), cex = 2.5,
       xlim = c(-0.5, 1), main = p, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]), col = col_discord, pch = 20)
  legend("topleft", legend = c(expression(paste("golden ", italic(p)[s], " only")),
                               expression(paste("controlled ", italic(p)[s], " only"))),
         pch = 20, col = c("#8B0000B3", "#00008AB3"), bty = "n")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}

## looking closer at sample size (ratio), effect estimates, p-values
discrepantDFpS <- subset(RProjects, pSsuccessG != pSsuccessC)
discrepantDFpS[, c("study", "project", "no", "nr", "c", "ro", "rr", "po1", "pr1", 
                 "psG", "psC"),]



# Exercise 2.2 (comparison two-trials rule vs sceptical p-value)
# ----------------------------------------------------------------------

## PART A: pS controlled vs two-trials rule

discrepantDF_ttrcontrolled <- subset(RProjects, successSignif != pSsuccessC)
discrepantDF_ttrcontrolled[, c("study", "project", "no", "nr", "c", "ro", 
                               "rr", "po1", "pr1", "psC"),]

## PART B: 

discrepantDF_ttrgolden <- subset(RProjects, successSignif != pSsuccessG)
discrepantDFpS[, c("study", "project", "no", "nr", "c", "ro", "rr", "po1", "pr1", 
                   "psG"),]


# Exercise 2.3 (p-value function and 95% sceptical confidence interval)
# ----------------------------------------------------------------------

## PARTS A and B


example <-RProjects[ind.example[3],]
with(example, 
     pValFunPlot(thetao = fiso, thetar = fisr, 
                 seo = se_fiso, ser = se_fisr))

with(example, 
     scepticalCI(thetao = fiso, thetar = fisr, 
            seo = se_fiso, ser = se_fisr))



##############################################################################
##########################  EXERCISE SESSION 3  ##############################
##############################################################################


# Exercise 3.1 (Power of the two-trials rule)
# ----------------------------------------------------------------------

par(mfrow = c(1, 1), las = 1)
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


# Exercise 3.2 (sample size with the two-trials rule)
# ----------------------------------------------------------------------

ss.cond1 <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                                  power = 0.8, designPrior = "conditional")
ss.pred1 <- sampleSizeSignificance(zo = p2z(po.plot, alternative = "one.sided"),
                                  power = 0.8, designPrior = "predictive" )
plot(po.plot, ss.cond1, type = "l", ylim = c(0, 4), col = "red", 
     xlab = "One-sided original p-value", ylab = "Relative sample size")
lines(po.plot, ss.pred1, col = "blue")
legend("topleft", c("Conditional", "Predictive"), col = c("red", "blue"), 
       lty = 1, bty = "n")

# Exercise 3.3 (comparing calculated sample sizes with used ones)
# ----------------------------------------------------------------------

par(mfrow = c(1, 2))
ss_rproj_cond <- with(RProjects, 
                      sampleSizeSignificance(zo = zo, 
                                             power = 0.9, alternative = "one.sided", 
                                             designPrior = "conditional"))

ss_rproj_pred <- with(RProjects, 
                      sampleSizeSignificance(zo = zo, 
                                             power = 0.9, alternative = "one.sided", 
                                             designPrior = "predictive"))

plot(RProjects$c, 
     ss_rproj_cond, pch = 19, 
     xlim = c(0.2, 30), 
     ylim = c(0.2, 30), 
     col = alpha("black", 0.5), 
     log = "xy", 
     xlab = "used in the projects", 
     ylab = "calculated with cond. power", 
     main = "Relative sample size")
abline(a = 0, b = 1, col = "red")


plot(RProjects$c, 
     ss_rproj_pred, pch = 19, 
     xlim = c(0.2, 30), 
     ylim = c(0.2, 30), 
     col = alpha("black", 0.5), 
     log = "xy", 
     xlab = "used in the projects", 
     ylab = "calculated with pred. power", 
     main = "Relative sample size")
abline(a = 0, b = 1, col = "red")

##############################################################################
##########################  EXERCISE SESSION 4  ##############################
##############################################################################


# Exercise 4.1 (power with the sceptical p-value)
# ----------------------------------------------------------------------


par(las = 1)

pow.cond2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "conditional", 
                                     alternative = "one.sided")

pow.pred2 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "predictive", 
                                     alternative = "one.sided")
pow.cond3 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "conditional", 
                                     alternative = "one.sided", 
                                     type = "controlled")

pow.pred3 <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                     c = 1, designPrior = "predictive", 
                                     alternative = "one.sided", 
                                     type = "controlled")

plot(po.plot, pow.cond3*100, type = "l", col = "red", ylim = c(0,100), 
     lwd = 2.1, 
     xlab = "One-sided original p-value", 
     ylab = "Power (in %)", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lty = 1, main = "c = 1", 
     cex.main = 2)
lines(po.plot, pow.cond1*100, col = "blue", lwd = 2.1, lty = 1)
lines(po.plot, pow.cond2*100, col = "orange", lwd = 2.1, lty = 1)
abline(h = 50, lty = 3, lwd = 1.5)

axis(2, at = 50, label = "50", cex.axis = 1.5,
     col = "gray40")

legend("bottomleft", 
       c("two-trials rule", expression(paste("controlled ", p[S])), 
         expression(paste("golden ", p[S]))), 
       col = c("blue", "red", "orange"), 
       lty = 1.5, bty = "n")



## for c = 2

par(las = 1)
pow.cond1_bis <- powerSignificance(zo = p2z(po.plot, alternative = "one.sided"), 
                                   c = 2, designPrior = "conditional", 
                                   alternative = "one.sided")
pow.cond2_bis <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         c = 2, designPrior = "conditional", 
                                         alternative = "one.sided")

pow.pred2_bis <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         c = 2, designPrior = "predictive", 
                                         alternative = "one.sided")
pow.cond3_bis <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         c = 2, designPrior = "conditional", 
                                         alternative = "one.sided", 
                                         type = "controlled")

pow.pred3_bis <- powerReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         c = 2, designPrior = "predictive", 
                                         alternative = "one.sided", 
                                         type = "controlled")

plot(po.plot, pow.cond3_bis*100, type = "l", col = "red", ylim = c(0,100), 
     lwd = 2.1, 
     xlab = "One-sided original p-value", 
     ylab = "Power (in %)", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lty = 1, main = "c = 2", 
     cex.main = 2)
lines(po.plot, pow.cond1_bis*100, col = "blue", lwd = 2.1, lty = 1)
lines(po.plot, pow.cond2_bis*100, col = "orange", lwd = 2.1, lty = 1)
abline(h = 50, lty = 3, lwd = 1.5)

axis(2, at = 50, label = "50", cex.axis = 1.5,
     col = "gray40")



legend("bottomleft", 
       c("two-trials rule", expression(paste("controlled ", p[S])), 
         expression(paste("golden ", p[S]))), 
       col = c("blue", "red", "orange"), 
       lty = 1.5, bty = "n")





# Exercise 4.2 (sample size with the sceptical p-value)
# ----------------------------------------------------------------------

par(mfrow = c(1,1),
    las = 1)

ss.cond2 <- sampleSizeReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         power = 0.8,
                                         alternative = "one.sided",
                                         type = "golden",
                                         designPrior = "conditional",
                                         level = 0.025)
ss.cond3 <- sampleSizeReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         power = 0.8,
                                         alternative = "one.sided",
                                         type = "controlled",
                                         designPrior = "conditional",
                                         level = 0.025)

ss.pred2 <- sampleSizeReplicationSuccess(zo = p2z(po.plot, alternative = "one.sided"),
                                         power = 0.8,
                                         alternative = "one.sided",
                                         designPrior = "predictive",
                                         type = "golden",
                                         level = 0.025)

plot(po.plot,
     ss.cond3,
     type = "l",
     col = "red",
     ylim = c(0, 12),
     cex.lab = 1.5,
     cex.axis = 1.5,
     ylab = "Relative sample size",
     xlab = "One-sided original p-value",
     lwd = 2)

lines(po.plot,
      ss.cond1,
      type = "l",
      col = "blue",
      lwd = 2)
lines(po.plot, 
      ss.cond2, 
      type = "l", 
      col = "orange", 
      lwd = 2)

legend("topleft", 
       c("two-trials rule", expression(paste("controlled ", p[S])), 
         expression(paste("golden ", p[S]))),
       col = c("blue", "red", "orange"), 
       lty = 1.5, bty = "n")

