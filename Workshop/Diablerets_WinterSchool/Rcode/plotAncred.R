
## computes posterior expectation in normal-normal model
postE <- function(Edat, VARdat, Eprior, VARprior){
  PRECdat <- 1/VARdat
  PRECprior <- 1/VARprior
  PRECposterior <- PRECdat + PRECprior
  Eposterior <- weighted.mean(x=c(Edat, Eprior), w=c(PRECdat, PRECprior))
  return(Eposterior)
}

## computes posterior variance in normal-normal model
postVAR <- function(VARdat, VARprior){
  PRECdat <- 1/VARdat
  PRECprior <- 1/VARprior
  PRECposterior <- PRECdat + PRECprior
  return(1/PRECposterior)
}

## computes lower and upper limit of posterior credible interval in normal-normal model
ULpost <- function(Edat, VARdat, Eprior, VARprior, alpha=0.05){
  Eposterior <- postE(Edat, VARdat, Eprior, VARprior)
  SDposterior <- sqrt(postVAR(VARdat, VARprior))
  t <- qnorm(1-alpha/2)
  U <- Eposterior + t*SDposterior
  L <- Eposterior - t*SDposterior
  return(cbind(U, L))
}


###################################################  
## format p-value compactly
fp <- function(p) formatPval(p)
###################################################
## computes sceptical limit
SL <- function(U, L, ratio=FALSE){
  if(ratio==TRUE){
    L <- log(L)
    U <- log(U)
  }        
  res <- (U-L)^2/(4*sqrt(U*L))
  if(ratio==TRUE)
    res <- exp(res)
  return(res)
}

## computes sceptical prior variance 
SLvar <- function(U, L, ratio=FALSE, alpha=0.05){
  mySL <- SL(U, L, ratio=FALSE)
  t <- qnorm(1-alpha/2)
  VARprior <- (mySL/t)^2
  return(VARprior)
}

## computes posterior variance in normal-normal model
postVAR <- function(VARdat, VARprior){
  PRECdat <- 1/VARdat
  PRECprior <- 1/VARprior
  PRECposterior <- PRECdat + PRECprior
  return(1/PRECposterior)
}
## computes upper and lower limit of confidence interval
UL <- function(theta, se, alpha=0.05){
  t <- qnorm(1-alpha/2)
  U <- theta + t*se
  L <- theta - t*se
  return(c(U, L))
}

library(biostatUZH, warn.conflicts = FALSE, quietly = TRUE)
#######################################################
#######################################################
#######################################################
plotAnCred <- function(ULori, 
                       ULrep, 
                       alpha = 0.05, 
                       myylim, 
                       corr = FALSE, 
                       pS = TRUE, 
                       showc = FALSE, 
                       colo = c(1, 3, 4, 2), # color of ori, post, prior, rep
                                             # 0 for white
                       colo2 = 1){ #color or arrows, text
  
  Uori <- ULori[1]
  Lori <- ULori[2]
  Urep <- ULrep[1]
  Lrep <- ULrep[2]
  
  Erep <- (Lrep+Urep)/2
  mySL <- SL(Uori, Lori)
  
  Lprior <- -mySL
  Uprior <- mySL
  Eprior <- 0
  
  t <- qnorm(1-alpha/2)
  VARprior <- (mySL/t)^2
  VARrep <- ((Urep-Lrep)/(2*t))^2
  VARori <- ((Uori-Lori)/(2*t))^2
  VARpost <- postVAR(VARdat=VARori, VARprior=VARprior)
  pRep <- z2p(abs(Erep/sqrt(VARrep)), alternative="greater") 
  mylen <- 21
  n <- 2*mylen-1
  al <- seq(.8,.4,len=mylen)
  cols <- rainbow(mylen, start = 0, end = 1/6, alpha = al)
  colours <- c(rev(cols[-mylen]), cols)
  fanP <- seq(alpha/2,1-alpha/2,len=n)
  fanQuant <- qnorm(fanP, Eprior, sqrt(VARprior))
  fanMat <- matrix(fanP, ncol = n, nrow = 2, byrow = TRUE)
  
  if(Erep>0){
    Lpost <- 0
    Epost <- t*sqrt(VARpost)
    Upost <- 2*Epost
  }
  if(Erep<0){
    Upost <- 0
    Epost <- -t*sqrt(VARpost)
    Lpost <- 2*Epost
  }
  
  Eori <- (Lori+Uori)/2
  VARori <- ((Uori-Lori)/(2*t))^2
  pOri <- z2p(abs(Eori/sqrt(VARori)), alternative="greater")
  myc <- VARori/VARrep
  mypS <- pSceptical(zo=Eori/sqrt(VARori), zr=Erep/sqrt(VARrep), c=myc)
  
  
  par(las=1)
  myylab <- ifelse(corr==FALSE, "Effect Size", "Correlation")
  
  plot(0, 0, 
       xlim = c(0.6,4.4), 
       ylim = myylim, 
       type = "n", 
       ylab = myylab, 
       xlab = "", 
       axes = FALSE)
  
  x1 <- seq(1, 2, len=n)
  y1 <- seq(0, Uori, len=n)
  
  fanMat <- matrix(0, ncol = length(y1), nrow = length(x1), byrow = TRUE)
  eps <- 1
  cols2 <- gray(c(1:n)/n)
  colours2 <- c(rev(cols2[-mylen]), cols2)
  
  for(i in 1:length(x1)){
    myvar <- exp(seq(log(1/eps), log(VARprior), len=length(x1)))
    myE <- postE(Edat=Eori, VARdat=VARori, Eprior, VARprior=myvar[i])
    myVAR <- postVAR(VARdat=VARori, VARprior=myvar[i])
    fanMat[i,] <- pnorm(y1, myE, sqrt(myVAR))
  }
  
  abline(v=3, lty=1, col= ifelse(colo[3] == 2, "grey", "white"))
  
  ## fanplot for posterior in interval 2-3
  x2 <- seq(2, 3, len=n)
  y2 <- seq(min(Lprior, Lpost), max(Uprior, Upost), len=n)
  for(i in 1:length(x2)){
    myvar <- exp(seq(log(VARori), log(1/eps), len=length(x2)))
    myE <- postE(Edat=Eori, VARdat=myvar[i], Eprior=Eprior, VARprior=VARprior)
    myVAR <- postVAR(VARdat=myvar[i], VARprior=VARprior)
    fanMat[i,] <- pnorm(y2, myE, sqrt(myVAR))
  }
  
  if(corr==FALSE)
    axis(2)
  if(corr==TRUE){
    where <- seq(-.9, .9, .2)
    axis(2, at=fisher(where), labels=(where))
  }
  axis(1, at=1, cex.axis=0.85, 
       labels=(c("Original Study")), 
       col.axis = colo2[1], 
       col = colo2[1])
  axis(1, at = 2, cex.axis = 0.85, 
       labels = "Posterior", 
       col.axis = colo2[2], 
       col = colo2[2])
  axis(1, at = 3, cex.axis = 0.85, 
       labels = "Sufficiently Sceptical Prior",
       col.axis = colo2[3], 
       col = colo2[3])
      
  axis(1, at = 4, cex.axis = 0.85, 
       labels = "Replication Study", 
       col.axis = colo2[4], 
       col = colo2[4])
  box()
  abline(h=0, lty=2)

  if(showc)
    text(0.45, -0.65, cex=1.25, paste("c=", as.character(myc), sep=""), pos=4)
  if(pS)
    text(2, 1.0, cex=1.5, substitute(paste(italic(p)[S],"=",pval), 
                                     list(pval=fp(mypS))), pos=3)
  
  text(1, Lori - 0.1, substitute(paste(italic(p)[o],"=",pval), 
                                 list(pval=fp(pOri))), pos=1)
  text(1, Eori, substitute(paste(hat(theta)[o],"=",thetaoval), 
                           list(thetaoval=round(Eori, 2))), pos=2)
  text(4, Lrep - 0.1, substitute(paste(italic(p)[r],"=",pval), 
                                 list(pval=fp(pRep))), pos=1)
  text(4, Erep, substitute(paste(hat(theta)[r],"=",thetarval), 
                           list(thetarval=round(Erep, 2))), pos=2)
  mygrey <- "gray30"
  v <- 1
  gplots::plotCI(x = (c(Eori,Epost, Eprior,Erep)), 
                 li =(c(Lori,Lpost, Lprior,Lrep)), 
                 ui = (c(v,v,1,v)*c(Uori,Upost, Uprior,Urep)), 
                 col = colo, 
                 lwd = 2, 
                 add = TRUE, 
                 sfrac = 0)
  points(1, Lori, pch = "-", cex = 2, col = colo[1])
  points(2, Lpost, pch = "-", cex = 2, col = colo[2])
  points(3, Lprior, pch = "-", cex = 2, col = colo[3])
  points(4, Lrep, pch = "-", cex = 2, col = colo[4])
  points(1, Uori, pch = "-", cex = 2, col = colo[1])
  points(2, Upost, pch="-", cex = 2, col = colo[2])
  points(3, Uprior, pch="-", cex = 2, col = colo[3])
  points(4, Urep, pch="-", cex = 2, col = colo[4])
  

  eps <- .07
  arrows(1.9, -0.2, 2, -0.05, cex=0.8, length=0.05, 
         col = colo2[2])
  text(1.55, -0.2, "fixed at zero", cex=0.9, 
       col = colo2[2])
  text(2.4, 1.05, "Reverse-Bayes analysis", cex=0.8, 
       col = colo2[3])
  arrows(1.9, 0.9, 2.9, 0.9,  cex=1.2, length=0.1, col = colo2[3])
  text(3.675, 1.05, "Assessing prior-data conflict", cex=0.8,
       col = colo2[5])
  arrows(3.1, 0.9, 4.25, 0.9,  cex=1.2, length=0.1, 
         col = colo2[5])
  
}
