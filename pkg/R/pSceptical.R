

zSceptical <- function(to, tr, c){
    ## arithmetic mean
    arit.mean <- function(x, y)
        return((x+y)/2)
    ## harmonic mean
    harm.mean <- function(x, y)
        return(2/(1/x+1/y))
    
    t2H <- harm.mean(to^2, tr^2)
    t2A <- arit.mean(to^2, tr^2)
    if(length(c)==1){
        if(c==1)
            t2 <- t2H/2
        else
            t2 <- (sqrt(t2A*(t2A+(c-1)*t2H))-t2A)/(c-1)
    }
    if(length(c)>1)
        t2 <- ifelse(c==1, t2H/2, (sqrt(t2A*(t2A+(c-1)*t2H))-t2A)/(c-1))
    return(sqrt(t2))
}

pSceptical <- function(po=NULL, pr=NULL, to=p2t(po, alternative=alternative),
                       tr=p2t(pr, alternative=alternative),
                       c, alternative="two.sided"){
    t <- zSceptical(to, tr, c)
    res <- t2p(t)
    if(alternative=="one.sided")
        res <- ifelse(sign(to)==sign(tr), res/2, 1-res/2)
    return(res)
}
