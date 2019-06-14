

pIntrinsic <- function(p=t2p(t, alternative=alternative),
                       t=NULL, alternative="two.sided", type="Held"){
    if(min(p) <=0 || max(p) >1)
        stop("All elements of p must lie in (0,1]!")
    if(type=="Held"){
        t <- p2t(p, alternative=alternative)/sqrt(2)
        iP <- t2p(t, alternative=alternative)
    }
    if(type=="Matthews"){
        t <- p2t(p, alternative=alternative)/sqrt(2)*sqrt(sqrt(5)-1)
        iP <- t2p(t, alternative=alternative)
    }
    return(iP)
}
