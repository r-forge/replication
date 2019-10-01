
zBox <- function(to, tr, c, sig.level){
    z <- qnorm(1-sig.level/2, lower.tail=TRUE)
    den <- ifelse((to^2 > z^2), c/(to^2/z^2-1)+1, NA)
    res <- tr^2/den
    return(sqrt(res))
}

pBox <- function(po=NULL, pr=NULL, to=p2t(po, alternative=alternative),
                 tr=p2t(pr, alternative=alternative), c, sig.level=0.05,
                 alternative="two.sided"){
    t <- zBox(to, tr, c, sig.level)
    res <- t2p(t)
    if(alternative=="one.sided")
        res <- ifelse(sign(to)==sign(tr), res/2, 1-res/2)
    return(res)
}
