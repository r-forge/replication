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

## range for d based on sceptical BF
## =============================================================================
## helper function to compute sufficiently sceptical relative prior variance
vss <- function (x, gamma, jeffreys = FALSE) {
  vssVec <- mapply(FUN = function(x, gammma, jeffreys) {
    if (!is.numeric(x))
      stop("x must be numeric")
    if (!is.numeric(gamma))
      stop("gamma must be numeric")
    if (gamma <= 0 || gamma > 1)
      stop("gamma must be in (0, 1]")
    if (!is.logical(jeffreys))
      stop("jeffreys must be logical")
    if (abs(x) <= 1) {
      if (gamma == 1)
        out <- 0
      else out <- NaN
    }
    else {
      y <- -x^2 * exp(-x^2)/gamma^2
      if (jeffreys == FALSE) {
        out <- -x^2/lamW::lambertWm1(x = y) - 1
      }
      else {
        out <- -x^2/lamW::lambertW0(x = y) - 1
      }
      if (!is.nan(out) && out < 0)
        out <- NaN
    }
    return(out)
  }, x, gamma, jeffreys)
  return(vssVec)
}

## function to compute range for d
dBFs <- function (zo, c = Inf, gamma, paradox = FALSE) {
  drangeVec <- mapply(FUN = function(zo, c, gamma, paradox) {
    if (!is.numeric(gamma))
      stop("gamma must be numeric")
    if (gamma <= 0 | gamma > 1)
      stop("gamma must be in (0, 1]")
    if (!is.numeric(zo))
      stop("zo must be numeric")
    if (!is.numeric(c))
      stop("c must be numeric")
    if (c < 0)
      stop("c cannot be smaller than 0")
    if (!is.logical(paradox))
      stop("paradox must be logical")
    g <- vss(x = zo, gamma = gamma)
    if (is.nan(g)) {
      dmin <- NaN
      dmax <- NaN
    }
    else {
      if (isTRUE(all.equal(g, 1))) {
        if (paradox == FALSE) {
          dmin <- ((-log(2)/zo^2 + 0.5) * (1/c + 1) +
                     1) * 0.5
          dmax <- Inf
        }
        else {
          dmin <- NaN
          dmax <- NaN
        }
      }
      else {
        A <- log((1/c + 1)/(1/c + g)/(1 + g))/zo^2 +
          g/(1 + g) + 1/(1 - g)
        B <- (1 - g)/(1/c + g)/(1/c + 1)
        M <- (1/c + g)/(g - 1)
        lower <- M - sqrt(A/B)
        upper <- M + sqrt(A/B)
        if (g < 1) {
          if (paradox == FALSE) {
            dmin <- upper
            dmax <- Inf
          }
          else {
            dmin <- -Inf
            dmax <- lower
          }
        }
        else {
          if (paradox == FALSE) {
            dmin <- lower
            dmax <- upper
          }
          else {
            dmin <- NaN
            dmax <- NaN
          }
        }
      }
    }
    drange <- c("dmin" = dmin, "dmax" = dmax)
    return(drange)
  }, zo, c, gamma, paradox)
  return(t(drangeVec))
}

## plot representing replication success/significance
## as a function of original $p$-value 
plotLevel <- function(c = Inf,
                      level = 0.025,
                      type = "golden",
                      method = "RS", # "RS", "signif", "meta", or "BFs"
                      alternative = "one.sided",
                      minRES = TRUE,
                      title = TRUE){ # showing minRES (only for RS)

  require(scales)
  require(ReplicationSuccess)
  par(las = 1, mgp = c(3, 0.6, 0), pty = "s")

  myplim <- c(0, 0.06)
  mydlim <- c(0, 3)
  alpha <- level
  eps <- 10e-10
  pos <- seq(10e-200, levelSceptical(level = level,
                                     alternative = alternative,
                                     type = type) - eps, length.out = 1000)
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
  }

  ## sceptical p-value
  if(method == "RS"){

    ## compute for c = infty
    minres <- effectSizeReplicationSuccess(zo = zos,
                                           c = Inf,
                                           level = level,
                                           alternative = alternative,
                                           type = type)
    ## compute for finite c
    resultsS <- matrix(NA, nrow = length(zos), ncol = length(c))
    for(i in 1:length(c)){
      resultsS[, i] <- effectSizeReplicationSuccess(zo = zos,
                                                    c = c[i],
                                                    level = level,
                                                    alternative = alternative,
                                                    type = type)
    }

    x <- c(pos, rep(min(pos), 2))
    z <- resultsS[, length(c)]
    y <- c(z, max(z), min(z))
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
                                                    level = level,
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

      polygon(c(0, pos, max(pos), max(pos), 0), c(0, minres, max(minres), 0, 0),
              density = 5,
              border = NA,
              col = scales::alpha("red", 0.5))
      if(type=="nominal")
        text(0.04, 1.5, labels="Never success", col=2)
      if(type=="golden")
        text(0.04, .5, labels="Never success", col=2)

    }
    axis(3, at = level, label = as.character(format(level, nsmall=2, digits=2)),
         col = "black",
         col.axis = "black",
         cex.lab = 0.8, cex.axis = 0.8)
  } else if (method == "BFs"){
    ## compute BF level based on significance level (with unit variance
    ## recalibration)
    zalpha <- qnorm(p = 1 - level)
    gammaS <- exp(1/2)*zalpha*exp(-0.5*zalpha^2) ## gamma_S = BF0:S(zalpha, g = 1
##     gammaS <- sqrt(2)*exp(-0.25*zalpha^2) ## gamma_S = BF0:S(zalpha, g = 1

    ## compute success region for c = infty
    if(minRES == TRUE){
      dminmax_cInf <- dBFs(zo = zos, c = Inf, gamma = gammaS, paradox =  FALSE)
      limitInd <- which.min(dminmax_cInf[,2]) ## last p-value where success possible
      x <- c(pos[1:limitInd], rev(pos[1:limitInd]))
      y <- c(dminmax_cInf[1:limitInd,2], rev(dminmax_cInf[1:limitInd,1]))
      lines(x, y,
            lty = 2,
            col = scales::alpha("red", 0.4),
            lwd = 3)

      xpoly <- c(0, pos[1:limitInd], pos[limitInd+1], max(pos), max(pos), 0)
      ypoly <- c(0, 0, dminmax_cInf[2:limitInd,1], 100, 100, 0, 0)
      polygon(xpoly, ypoly,
              density = 5,
              border = NA,
              col = scales::alpha("red", 0.5))
        text(0.04, .5, labels="Never success", col=2)
    }

    ## compute success regions for finite c
    where <- which.min((pos - 0.005)^2)
    for(i in 1:length(c)){
      dminmax_c <- dBFs(zo = zos, c = c[i], gamma = gammaS, paradox =  FALSE)
      limitIndc <- which.min(dminmax_c[,2]) ## last p-value where success possible

      x <- c(pos[1:limitIndc], rev(pos[1:limitIndc]))
      y <- c(dminmax_c[1:limitIndc,2], rev(dminmax_c[1:limitIndc,1]))
      lines(x, y,
            type = "l",
            lty = 1,
            col = scales::alpha("black", 0.5))
      text(pos[where], dBFs(zo = zos[where], c = c[i], gamma = gammaS, paradox =  FALSE),
           paste0("c = ", as.character(c[i])),
           col = scales::alpha("black", 0.5),
           cex = 0.8)

    }

    ## first p-value where finite upper limit on replication success
    limitInd2 <- which.max(pos[!is.nan(dminmax_c[,2]) & dminmax_c[,2] < Inf])
    xpoly <- c(0, pos[1:limitInd], rev(pos[limitInd2:limitIndc]), 0, 0)
    ypoly <- c(0, 0, dminmax_c[2:limitIndc,1],
               rev(dminmax_c[limitInd2:limitIndc,2]), 100, 0)
    polygon(x, y, col="lightgreen")
    lines(x, y, type = "l")

    } else {

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
      posigplot <- seq(1e-20, level, length.out = 1000)
    if(method=="meta")
      posigplot <- seq(1e-20, 0.06, length.out = 1000)
    zosigplot <- p2z(p = posigplot, alternative = "one.sided")

    x <- c(posigplot, max(posigplot), rep(min(posigplot), 2))
    ##           z <- resultsS[, length(c)]
    z <- dSignifBound(zo = zosigplot, c = myc_comp[length(myc_comp)],
                      alpha = level, method=method)
    y <- c(z, rep(max(mydlim),2), min(mydlim))
    polygon(x, y, col="lightgreen")
    text(0.0125, 2.5, labels="Success", col="darkgreen")
    for(i in 1:length(myc_comp)){
      lines(posigplot,
            dSignifBound(zo = zosigplot, c = myc_comp[i], alpha = level, method=method),
            type = "l",
            lty = 1,
            col = alpha("black", 0.5))
      where <- which.min((posigplot - 0.015)^2)
      text(x = posigplot[where],
           y = dSignifBound(zo = zosigplot[where], c = myc_comp[i], alpha = level, method=method),
           labels = paste0("c=",as.character(myc_comp[i])),
           col = alpha("black", 0.5),
           cex = 0.8)

    }

    axis(3, at = level, label = as.character(format(level, nsmall=2, digits=2)), col = "black",
         col.axis = "black",
         cex.axis=0.8, cex.lab = 0.8)
  }

  ## surround with box
  box()
}
