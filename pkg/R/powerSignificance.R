
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

powerSignificance <- function(po = NULL, to = p2t(po, alternative = alternative),
                              c = 1, level = 0.05,
                              designPrior = "conditional",
                              alternative = "two.sided",
                              d = 0,
                              shrinkage = 1,
                              EB = FALSE){
    v <- p2t(level, alternative = alternative)
    lowertail <- ifelse(alternative == "less", TRUE, FALSE)
    # to <- abs(to)
    if (EB == TRUE){
        s <- pmax(1 - (1 + d)/to^2, 0)
        pSig <- pnorm(q = v, mean = s*to*sqrt(c),
                      sd = sqrt(s*c*(1 + d) + 1 + d*c),
                      lower.tail = lowertail)
        return(pSig)
    }
    if(designPrior == "conditional")
        pSig <- pnorm(q = v, mean = shrinkage*to*sqrt(c), 
                      lower.tail = lowertail)
    if(designPrior == "predictive"){
        pSig <- pnorm(q = v, mean = shrinkage*to*sqrt(c),
                      sd = sqrt(c + 1 + 2*d*c),
                      lower.tail = lowertail)
    }
    return(pSig)
}