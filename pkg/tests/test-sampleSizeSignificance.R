library(ReplicationSuccess)

zo <- seq(2, 4, 0.5)
apply_grid <- expand.grid(priors = c("conditional", "predictive", "EB"),
                          d = c(0, 0.1),
                          shrinkage = c(1, 0.75),
                          alt = c("greater", "two.sided"),
                          stringsAsFactors = FALSE)
for (i in seq(1, nrow(apply_grid))) {
  print(apply_grid[i,])
  c <- sampleSizeSignificance(zo = zo, 
                              power = 0.8, 
                              level = 0.05,
                              designPrior = apply_grid$priors[i],
                              alternative = apply_grid$alt[i],
                              d = apply_grid$d[i])
  print(round(c, digits = 5))
}