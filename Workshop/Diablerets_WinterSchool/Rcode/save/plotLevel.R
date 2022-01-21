## How large is zr required to achieve replication success?
zrHowLarge <- function(zo, c, level){ 
  z <- p2z(level, alternative = "one.sided")
  result <- z*sqrt(1 + c/(zo^2/z^2 - 1))
  return(result)
}

## How large is d required to achieve replication success with meta-analysis?
effectSizeMeta <- function(zo, c, level){ 
  z <- p2z(level, alternative = "one.sided")
  result <- sqrt(2/c)*z/zo
  return(result)
}

## plot representing replication success/significance
## as a function of original $p$-value 
    plotLevel <- function(c, 
                          level, 
                          type, 
                          method = "RS", # "RS" or "signif" or "meta"
                          alternative = "one.sided", 
                          minRES = T, title = TRUE){ # showing minRES (only for RS)
      
      require(scales)
         par(las = 1, mgp = c(3, 0.6, 0), pty="s")
   
      myplim <- c(0, 0.06)
      mydlim <- c(0, 3)
      alpha = level
      eps = 10e-10
      pos = seq(10e-200, levelSceptical(level = level, 
                                        alternative = alternative, 
                                        type = type) - eps, length.out=1000)
      zos <- p2z(pos, alternative = alternative)
      if(method == "RS"){
      minres = effectSizeReplicationSuccess(zo = zos, 
                                            c = Inf, 
                                            level = level, 
                                            alternative = alternative, 
                                            type = type)
      
      resultsS <- matrix(NA, nrow = length(zos), ncol = length(c))
      
      for(i in 1:length(c)){
        resultsS[, i] <- effectSizeReplicationSuccess(zo = zos, 
                                                      c = c[i], 
                                                      level = level, 
                                                      alternative = alternative, 
                                                      type = type)
      }
      
      
      ## plot 1: po vs d (replication success)

      
      plot(NULL , 
           xlim = myplim, 
           ylim = mydlim, 
           xlab = expression(paste("Original ", italic(p), "-value ", italic(p)[o])), 
           ylab = expression(paste("Relative effect size ", italic(d))),
           xaxs="i", yaxs="i", cex.axis=1) ##, axes=FALSE)
      if(title==TRUE)
          title("Sceptical p-value")
      
      
      
      box()
      abline(h = 1, 
             col = "grey", ## "lightgrey", 
             lty = 2)
      
      x <- c(pos, rep(min(pos), 2))
      z <- resultsS[, length(c)]
      y <- c(z, max(z), min(z))
      polygon(x, y, col="lightgreen")
      if(type=="nominal")
          text(0.009, 2.5, labels="Success", col="darkgreen")
      if(type=="golden")
          text(0.02, 2.5, labels="Success", col="darkgreen")
      abline(v = c(alpha, levelSceptical(level = level, 
                                         alternative = alternative, 
                                         type = type)), 
             col = "grey", ## c("lightgrey", "lightgrey"), 
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
                                                      level = level,
                                                      alternative = alternative,
                                                      type = type),
             paste0("c = ", as.character(c[i])),
             col = scales::alpha("black", 0.5),
             cex = 0.8)
        
      }
      
      if(minRES == T){
        lines(pos, minres, lty = 2, 
              col = scales::alpha("red", 0.4), 
              lwd = 3)
                
        polygon(c(pos[pos >= 0], 100), c(minres[pos >= 0], 0), 
                density = 5, 
                border = NA, 
                col = scales::alpha("red", 0.5))
        if(type=="nominal")
            text(0.04, 1.5, labels="Never success", col=2)
        if(type=="golden")
            text(0.04, .5, labels="Never success", col=2)
        
      }
      axis(3, at = level, label = as.character(format(level, nsmall=2, digits=2)), col = "black", 
           col.axis = "black", 
           cex.lab = 0.8, cex.axis=0.8)
      # pch = 20)
      } else {

          if(method=="signif")
              mymain <- "Two-trials rule"
          if(method=="meta")
              mymain <- "Meta-analysis"
      plot(NULL ,
           xlim = myplim,
           ylim = mydlim,
           xlab = expression(paste("Original ", italic(p), "-value ", italic(p)[o])),
           ylab = expression(paste("Relative effect size ", italic(d))),
           main = mymain, xaxs="i", yaxs="i", cex.axis=1)
            
      abline(h = 1,
             col = "grey", ## "lightgrey",
             lty = 2)

          if(method=="signif"){
              polygon(x = c(alpha, alpha, 0.15, 0.15),
              y = c(-0.15, 3.5, 3.5, -0.15),
              density = 5,
              border = NA,
              col = alpha("red", 0.5))
              text(0.0425, 1.5, labels="Never success", col=2)
              abline(v = alpha,
                     col = alpha("grey"),
                     lty = 2,
                     lwd = c(3, 1))
          }

      dSignifBound <- function(zo, c, alpha, method="signif") {
        zalpha <- qnorm(p = alpha, lower.tail = FALSE)
        if(method=="signif"){
            zalpha <- qnorm(p = alpha, lower.tail = FALSE)
            d <- zalpha/(zo*sqrt(c))
        }
        if(method=="meta"){
            zalpha <- qnorm(p = alpha^2, lower.tail = FALSE)
            d <- (sqrt(2)*zalpha/zo-1)/sqrt(c)
        }
        return(d)
      }
      
      ## compute c-lines for two-trials rule plot
          myc_comp <- c
          if(method=="signif")
              posigplot <- seq(1e-20, 0.025, length.out = 1000)
          if(method=="meta")
              posigplot <- seq(1e-20, 0.06, length.out = 1000)
          zosigplot <- p2z(p = posigplot, alternative = "one.sided")

          x <- c(posigplot, max(posigplot), rep(min(posigplot), 2))
##           z <- resultsS[, length(c)]
          z <- dSignifBound(zo = zosigplot, c = myc_comp[length(myc_comp)],
                            alpha = 0.025, method=method)
      y <- c(z, rep(max(mydlim),2), min(mydlim))
      polygon(x, y, col="lightgreen") 
      text(0.0125, 2.5, labels="Success", col="darkgreen")
      for(i in 1:length(myc_comp)){
          lines(posigplot, 
                dSignifBound(zo = zosigplot, c = myc_comp[i], alpha = 0.025, method=method),
                type = "l",
                lty = 1,
                col = alpha("black", 0.5))
          where <- which.min((posigplot - 0.015)^2)
          text(x = posigplot[where],
               y = dSignifBound(zo = zosigplot[where], c = myc_comp[i], alpha = 0.025, method=method),
               labels = paste0("c=",as.character(myc_comp[i])),
               col = alpha("black", 0.5),
               cex = 0.8)


      }
          
          axis(3, at = level, label = as.character(format(level, nsmall=2, digits=2)), col = "black", 
           col.axis = "black", 
           cex.axis=0.8, cex.lab = 0.8)
      }
    }
    
