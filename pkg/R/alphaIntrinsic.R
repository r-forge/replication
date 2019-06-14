
alphaIntrinsic <- function(alpha, alternative="two.sided", type="Held"){
    t <- p2t(alpha, alternative=alternative)
    if(type=="Held")
        result <- t2p(sqrt(2)*t, alternative=alternative)
    if(type=="Matthews")
        result <- t2p(sqrt(2)*t/sqrt(sqrt(5)-1), alternative=alternative)
    return(result)
}
