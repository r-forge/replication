
t2p <- function(t, alternative="two.sided"){
    if(alternative=="two.sided")
        p <- 2*pnorm(abs(t), lower.tail=FALSE)
    # if(alternative=="one.sided")
    #     p <- pnorm(t, lower.tail=FALSE)
    if (alternative == "less")
        p <- pnorm(q = t, lower.tail = TRUE)
    if (alternative == "greater" | alternative == "one.sided")
        p <- pnorm(q = t, lower.tail = FALSE)
    return(p)
}
