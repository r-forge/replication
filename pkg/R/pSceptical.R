pSceptical <- function(zo,
                       zr,
                       c, 
                       alternative = "one.sided"){
    z <- zSceptical(zo, zr, c)
    res <- z2p(z = z)
    if(alternative == "one.sided")
        res <- ifelse(sign(zo) == sign(zr), res/2, 1 - res/2)
    return(res)
}
