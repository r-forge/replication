library(ReplicationSuccess)

zo <- seq(2, 4, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          d = c(0, 0.5),
                          c = c(0.5, 1, 2),
                          stringsAsFactors = FALSE)
for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  pis <- predictionInterval(zo = zo, 
                            c = apply_grid$c[i],
                            designPrior = apply_grid$priors[i],
                            d = apply_grid$d[i])
  print(round(pis, digits = 5))
}