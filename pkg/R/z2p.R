
z2p <- function(z, 
                alternative = "two.sided"){
    if (!(alternative %in% c("less", "greater", "two.sided", "one.sided")))
        stop('alternative must be either "less", "greater", "two.sided", or "one.sided"')
    
    if(alternative=="two.sided")
        p <- 2*pnorm(abs(z), lower.tail = FALSE)
    # if(alternative=="one.sided")
    #     p <- pnorm(q = z, lower.tail=FALSE)
    
    if (alternative == "less")
        p <- pnorm(q = z, lower.tail = TRUE)
    
    if (alternative == "greater" | alternative == "one.sided")
        p <- pnorm(q = z, lower.tail = FALSE)
    return(p)
}
