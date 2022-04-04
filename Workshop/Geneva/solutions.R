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

