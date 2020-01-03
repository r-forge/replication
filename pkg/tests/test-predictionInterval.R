library(ReplicationSuccess)

po <- seq(0.01, 0.05, 0.01)
to <- seq(2, 4, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          d = c(0, 0.5),
                          c = c(0.5, 1, 2),
                          alt = c("less", "greater", "two.sided"),
                          stringsAsFactors = FALSE)
for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  pis <- predictionInterval(po = po, 
                            c = apply_grid$c[i],
                            designPrior = apply_grid$priors[i],
                            alternative = apply_grid$alt[i],
                            d = apply_grid$d[i])
  print(round(pis, digits = 5))
  pis <- predictionInterval(to = to, 
                            c = apply_grid$c[i],
                            designPrior = apply_grid$priors[i],
                            alternative = apply_grid$alt[i],
                            d = apply_grid$d[i])
  print(round(pis, digits = 5))
}

# par(mfrow = c(2, 2), las = 1)
# for (p in unique(ReplicationProjects$project)) {
#   data_project <- subset(ReplicationProjects, project == p)
#   to <- data_project$z_O/data_project$z_se_O
#   tr <- data_project$z_R/data_project$z_se_R
#   c <- data_project$z_se_O^2/data_project$z_se_R^2
#   PI <- predictionInterval(to = to, c = c)
#   color <- ifelse(tr < PI$upper & tr > PI$lower, "#333333B3", "#8B0000B3")
#   coverage <- mean(color == "#333333B3")
#   plot(to, tr, las = 2, ylim = c(-4, 15), xlim = c(0, 22), col = color, pch = 20,
#        xlab = expression(italic(t)[o]), ylab = expression(italic(t)[r]), main = p)
#   text(x = 16.5, y = 12, paste0(round(coverage*100, 2), "% coverage"))
#   abline(h = 0, lty = 2, col = "gray70")
#   arrows(to, PI$lower, to, PI$upper, length = 0.02, angle = 90, code = 3, col = color)
# }