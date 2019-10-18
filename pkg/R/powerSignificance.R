
## function powerSignificant computes a generalisation of Goodman's pSig
## the probability of significance of a replication experiment
## to: test statistic of original study
## c: variance ratio sigma_o^2/sigma_r^2 = n_r/n_o
## d: relative between study heterogeneity tau^2/sigma_o^2

# powerSignificance <- function(po=NULL, to=p2t(po, alternative=alternative),
#                               c=1, level=0.05,
#                               designPrior="conditional",
#                               alternative="two.sided"){
#     v <- p2t(level, alternative=alternative)
#     to <- abs(to)
#     if(designPrior=="conditional")
#         pSig <- pnorm(to*sqrt(c) - v)
#     if(designPrior=="predictive")
#         pSig <- pnorm(v, mean=sqrt(c)*to, sd=sqrt(1+c), lower.tail=FALSE)
#     return(pSig)
# }

powerSignificance <- function(po=NULL, to=p2t(po, alternative=alternative),
                              c=1, level=0.05,
                              designPrior="conditional",
                              alternative="two.sided",
                              shrinkage=FALSE,
                              d=0){
    v <- p2t(level, alternative=alternative)
    to <- abs(to)
    if(designPrior=="conditional")
        pSig <- pnorm(to*sqrt(c) - v)
    if(designPrior=="predictive"){
        if (shrinkage==TRUE) s <- pmax(1 - (1 + d)/to^2, 0)
        else s <- 1
        pSig <- pnorm(q = v, mean= s*to*sqrt(c), 
                      sd=sqrt(s*c*(1 + d) + 1 + d*c), 
                      lower.tail=FALSE)
    }
    return(pSig)
}
