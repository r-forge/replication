
deltaEB <- function(po=NULL, pr=NULL,
                    to=p2t(po, alternative=alternative),
                    tr=p2t(pr, alternative=alternative),
                    c, alternative="two.sided"){
    term1 <- (tr/sqrt(c)-to)^2
    term2 <- 1+1/c
    res <- pmax(term1, term2) - 1/c
    return(1/res)
}
    
