# functions used in workshop



# define ForestPlot function (nearly unchanged from hMeanChiSqCI.Rnw) ----------------------------------------------------

ForestPlot <- function(thetahat, se, 
                       barHeight = 0.5,
                       arrowHeight = 0.2,
                       lwd = 1.1,
                       cex = 2,
                       level = 0.95, 
                       studyNames = c("original", "replication"),
                       TTR = FALSE, 
                       title = NA, 
                       scepticalCI = FALSE, 
                       metaAn = FALSE){
  
  diamond <- function(center, height, width) {
    base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) 
    trans <- rbind(base[1,] * width, base[2,] * height) + center
    geom_polygon(data = as.data.frame(t(trans)), mapping = aes(x = V1, y = V2)) 
  }
  
  
  
  #############################################
    if(TTR == FALSE)
        c <- (se[1]/se[2])^2
    if(TTR == TRUE)
        c <- 0
    meta <- metagen(TE = thetahat, seTE = se, random = TRUE, level.ci = level)
    CIs <- scepticalCI(thetao = thetahat[1], thetar = thetahat[2], 
                       seo = se[1], ser = se[2], c = c,
                        alternative = "two.sided", levelCI = level)

##    if(is.na(studyNames))
##        studies <- tibble(y = thetahat,
##                          lower = thetahat - 1.96 * se,
##                          upper = thetahat + 1.96 * se,
##                          names = paste("study", seq_along(y)))
##    else
        studies <- tibble(y = thetahat,
                          lower = thetahat - 1.96 * se,
                          upper = thetahat + 1.96 * se,
                          names = studyNames)
    
    
  hMean <- tibble(lower = CIs$CI[,1],
                  upper = CIs$CI[,2])
  
  fixedEffect <- tibble(lower = meta$lower.fixed,
                         upper = meta$upper.fixed)
  randomEffect <- tibble(lower = meta$lower.random,
                         upper = meta$upper.random)

    mmin <- min(studies$lower, fixedEffect$lower, hMean$lower)
    mmax <- max(studies$upper, fixedEffect$upper, hMean$upper)
    plotMargin <- mmin - 2.125 * (mmax - mmin)/3
    textStart <- mmin - 2 * (mmax - mmin)/3
  
  p <- ggplot(data = studies,
              mapping = aes(x=y, y=rev(seq_along(y)))) + 
    geom_errorbar(width = barHeight, xmin = studies$lower, xmax = studies$upper, size=lwd, col=c(5,6)) +
    geom_point(shape = 15, col=c(5,6)) +
    geom_vline(xintercept = 0, linetype = "dashed", col="darkgrey") +
    xlim(c(plotMargin, mmax)) +
    geom_text(aes(label = names), x = textStart, hjust = 0) + 
    xlab(expression(theta)) + 
    scale_y_continuous(expand = c(.1, .1)) +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust=0.5)) +
      ggtitle(title) 
  
  if(scepticalCI == TRUE){
  if(!is.null(hMean)){
    p <- p + geom_segment(data = hMean,
                          mapping = aes(x = lower, xend = upper),
                          y= 0, yend=0,
                          arrow = arrow(ends="both", length = unit(arrowHeight, "in")),
                          color="red", size=lwd) 
    p <- p + geom_text(label = "sceptical", x = textStart, y = 0, hjust = 0) +
      geom_point(data = studies, mapping = aes(x=y), y=0, color="blue", size=cex)
  }
  }
  
  if(metaAn == TRUE){
  if(!is.null(fixedEffect)){
    p <- p + diamond(center = c(mean(c(fixedEffect$lower, fixedEffect$upper)), ifelse(scepticalCI == TRUE, -1, 0)),
                     height = .2, width = (fixedEffect$upper-fixedEffect$lower)/2) 
    p <- p + geom_text(label = "meta-analysis", x = textStart, y =  ifelse(scepticalCI == TRUE, -1, 0) , hjust = 0)
  }
}
  
    return(p)
}
