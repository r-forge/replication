
## function powerSignificant computes a generalisation of Goodman's pSig
## the probability of significance of a replication experiment
## to: test statistic of original study
## c: variance ratio sigma_o^2/sigma_r^2 = n_r/n_o

powerSignificance <- function(po=NULL, to=p2t(po, alternative=alternative),
                              c=1, level=0.05,
                              designPrior="conditional",
                              alternative="two.sided"){
    v <- p2t(level, alternative=alternative)
    to <- abs(to)
    if(designPrior=="conditional")
        pSig <- pnorm(to*sqrt(c) - v)
    if(designPrior=="predictive")
        pSig <- pnorm(v, mean=sqrt(c)*to, sd=sqrt(1+c), lower.tail=FALSE)
    return(pSig)
}


