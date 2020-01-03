library(ReplicationSuccess)

po <- seq(0.01, 0.05, 0.01)
to <- seq(2, 4, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          d = c(0, 0.1),
                          shrinkage = c(1, 0.75),
                          alt = c("greater", "two.sided"),
                          stringsAsFactors = FALSE)
for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  c <- sampleSizeSignificance(po = po, 
                              power = 0.8, 
                              level = 0.05,
                              designPrior = apply_grid$priors[i],
                              alternative = apply_grid$alt[i],
                              d = apply_grid$d[i],
                              shrinkage = apply_grid$shrinkage[i])
  print(round(c, digits = 5))
  c <- sampleSizeSignificance(to = to, 
                              power = 0.8, 
                              level = 0.05,
                              designPrior = apply_grid$priors[i],
                              alternative = apply_grid$alt[i],
                              d = apply_grid$d[i])
  print(round(c, digits = 5))
}

# priors <- c("conditional", "predictive", "EB")
# plot(po, sampleSizeSignificance(po = po,
#                                 power = 0.8,
#                                 level = 0.05,
#                                 designPrior = priors[2],
#                                 alternative = "two.sided"),
#      type = "l", log = "y", ylab = "c", lty = 2, ylim = c(0.25, 50), las = 1)
# abline(h = 1, lty = 2, col = "lightgrey")
# lines(po, sampleSizeSignificance(po = po,
#                                  power = 0.8,
#                                  level = 0.05,
#                                  designPrior = priors[1],
#                                  alternative = "two.sided"),
#       lty = 1)
# lines(po, sampleSizeSignificance(po = po,
#                                  power = 0.8,
#                                  level = 0.05,
#                                  designPrior = priors[3],
#                                  alternative = "two.sided"),
#       lty = 3)
# legend("bottomleft", legend = priors,
#        lty = c(1, 2, 3), bty = "n", cex = 0.5)
