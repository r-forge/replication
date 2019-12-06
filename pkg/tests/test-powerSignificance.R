library(ReplicationSuccess)

po <- seq(0.01, 0.06, 0.01)
to <- seq(-3, 3, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          c = c(0.5, 1, 1.5, 2),
                          d = c(0, 1),
                          alt = c("two.sided", "greater"),
                          shrinkage = c(1, 0.5),
                          stringsAsFactors = FALSE)

for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  p <- powerSignificance(po = po, 
                         c = apply_grid$c[i],
                         level = 0.05,
                         designPrior = apply_grid$priors[i],
                         alternative = apply_grid$alt[i],
                         d = apply_grid$d[i],
                         shrinkage = apply_grid$shrinkage[i])
  print(round(p, digits = 5))
  p <- powerSignificance(to = to, 
                         c = apply_grid$c[i],
                         level = 0.05,
                         designPrior = apply_grid$priors[i],
                         alternative = apply_grid$alt[i],
                         d = apply_grid$d[i],
                         shrinkage = apply_grid$shrinkage[i])
  print(round(p, digits = 5))
}

# extending plots from Held (2020)
# par(mfrow = c(2, 2), las = 1)
# for (i in seq_along(c)) {
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[1],
#                          alternative = alt[1])
#   print(round(p, digits = 5))
#   plot(po, p, type = "l", ylim = c(0, 1), main = bquote(c == .(c[i])), ylab = "Power")
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[1],
#                          alternative = alt[2])
#   print(round(p, digits = 5))
#   lines(po, p, lty = 1, col = "darkgrey")
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[2],
#                          alternative = alt[1])
#   print(round(p, digits = 5))
#   lines(po, p, lty = 2)
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[2],
#                          alternative = alt[2])
#   print(round(p, digits = 5))
#   lines(po, p, lty = 2, col = "darkgrey")
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[3],
#                          alternative = alt[1])
#   print(round(p, digits = 5))
#   lines(po, p, lty = 3)
#   p <- powerSignificance(po = po, c = c[i], 
#                          designPrior = priors[3],
#                          alternative = alt[2])
#   print(round(p, digits = 5))
#   lines(po, p, lty = 3, col = "darkgrey")
#   legend("bottomleft", legend = priors, lty = c(1, 2, 3), bty = "n", cex = 0.5)
#   legend("topright", legend = c("two-sided", "one-sided"), 
#          lty = 1, col = c(1, "darkgrey"), bty = "n", cex = 0.5)
# }