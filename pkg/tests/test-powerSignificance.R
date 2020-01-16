library(ReplicationSuccess)

zo <- seq(-3, 3, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          c = c(0.5, 1, 1.5, 2),
                          d = c(0, 1),
                          alt = c("two.sided", "one.sided", "greater"),
                          shrinkage = c(1, 0.5),
                          stringsAsFactors = FALSE)

for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  p <- powerSignificance(zo = zo, 
                         c = apply_grid$c[i],
                         level = 0.05,
                         designPrior = apply_grid$priors[i],
                         alternative = apply_grid$alt[i],
                         d = apply_grid$d[i],
                         shrinkage = apply_grid$shrinkage[i])
  print(round(p, digits = 5))
}