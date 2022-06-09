
## plot representing replication success/significance
## as a function of original $p$-value 
plotLevel2 <- function(c = Inf,
                      SESI = 1.5,
                      sigmao = 1,
                      level = 0.025, 
                      method = "RS", # "RS", "signif", "meta", "BFr", or "BFs"
                      alternative = "one.sided",
                      minRES = TRUE,
                      title = TRUE){ # showing minRES (only for RS)
    type <- "golden"
  require(scales)
  require(ReplicationSuccess)
  par(las = 1, mgp = c(3, 0.6, 0), pty = "s")

  
  myupper <- round(levelSceptical(level)*100)/100
  myplim <- c(0, myupper)
#  myplim <- c(0, 0.06)
  mydlim <- c(0, 3)
  alpha <- level
  eps <- 10e-10
  pos <- seq(10e-200, levelSceptical(level = level,
                                     alternative = alternative,
                                     type = type) - eps, length.out = 1000)
##  pos <- seq(10e-200, myupper - eps, length.out = 1000)
  zos <- p2z(pos, alternative = alternative)

  ## plot skeleton for all methods
  plot(NULL,
       xlim = myplim,
       ylim = mydlim,
       xlab = expression(paste("Original ", italic(p), "-value ", italic(p)[o])),
       ylab = expression(paste("Relative effect size ", italic(d))),
       xaxs = "i", yaxs = "i", cex.axis = 1)
  abline(h = 1,
         col = "grey", ## "lightgrey",
         lty = 2)
  if (title == TRUE) {
    if(method == "signif") title("Two-trials rule")
    if(method == "meta") title("Meta-analysis")
    if (method == "RS") title("Sceptical p-value")
    if (method == "BFs") title("Sceptical Bayes factor")
    if (method == "BFr") title("Replication Bayes factor")
  }

  ## sceptical p-value
  if(method == "RS"){

      sigmao <- 1
      thetao <- sigmao*zos
      minres <- SESI/thetao
      mylevel <- levelEquivalent(dinf=minres, level=pos)

      
    ## compute for finite c
    resultsS <- matrix(NA, nrow = length(zos), ncol = length(c))
    for(i in 1:length(c)){
      resultsS[, i] <- effectSizeReplicationSuccess(zo = zos,
                                                    c = c[i],
                                                    level = mylevel,
                                                    alternative = alternative)
    }

    x <- c(pos, max(pos), min(pos), min(pos))
    z <- resultsS[, length(c)]
    y <- c(z, max(mydlim), max(mydlim), min(z))
    polygon(x, y, col="lightgreen")
    if(type == "nominal")
      text(0.009, 2.5, labels="Success", col="darkgreen")
    if(type == "golden")
      text(0.02, 2.5, labels="Success", col="darkgreen")
    abline(v = c(alpha, levelSceptical(level = level,
                                       alternative = alternative,
                                       type = type)),
           col = "grey",
           lty = 2)


    for(i in 1:length(c)){
      lines(pos,
            resultsS[, i],
            type = "l",
            lty = 1,
            col = scales::alpha("black", 0.5))
      where <- which.min((pos - 0.005)^2)
      text(pos[where], effectSizeReplicationSuccess(zo = zos[where],
                                                    c = c[i],
                                                    level = mylevel[where],
                                                    alternative = alternative,
                                                    type = type),
           paste0("c = ", as.character(c[i])),
           col = scales::alpha("black", 0.5),
           cex = 0.8)

    }

    if(minRES == TRUE){
      lines(pos, minres, lty = 2,
            col = scales::alpha("red", 0.4),
            lwd = 3)

      polygon(c(0, pos, levelSceptical(max(pos)), levelSceptical(max(pos)), 0), c(0, minres, max(minres), 0, 0),
              density = 5,
              border = NA,
              col = scales::alpha("red", 0.5))
      if(type=="nominal")
        text(0.04, 1.5, labels="Never success", col=2)
      if(type=="golden")
        text(0.04, .5, labels="Never success", col=2)

    }
    ## axis(3, at = level, label = as.character(format(level, nsmall=2, digits=2)),
    ##      col = "black",
    ##      col.axis = "black",
    ##      cex.lab = 0.8, cex.axis = 0.8)
      mypo <- c(0.01, 0.025, 0.04, 0.055)
      myzo <- c(NA, p2z(mypo, alternative="one.sided"))
      mythetao <- myzo*sigmao
      axis(3, at=c(NA, mypo), as.character(round(mythetao/SESI, 1)))
      mtext("relevance of original effect", side = 3, line = 1.5)

  }
  
  ## surround with box
  box()
}
