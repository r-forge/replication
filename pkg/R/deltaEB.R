
deltaEB <- function(zo,
                    zr,
                    c){
    term1 <- (zr/sqrt(c) - zo)^2
    term2 <- 1 + 1/c
    res <- pmax(term1, term2) - 1/c
    return(1/res)
}
    
