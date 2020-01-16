
thresholdSceptical <- function(level, 
                               alternative = "one.sided", 
                               type = "nominal"){
    if (!(type %in% c("nominal", "liberal", "controlled")))
        stop('type must be either "nominal", "liberal", or "controlled"')
    if (!(alternative %in% c("one.sided", "two.sided")))
        stop('type must be either "one.sided" or "two.sided"')
    
    if(type == "nominal")
        res <- level
    if(type == "liberal")
        res <- pIntrinsic(level, alternative = "one.sided")
    if(type == "controlled"){
        pI2 <- pIntrinsic(p = 2*level^2, alternative = "one.sided")
        res <- pIntrinsic(p = pI2, alternative = "one.sided")
    }
    if(alternative == "two.sided")
        res <- 2*res
    
    return(res)
}
