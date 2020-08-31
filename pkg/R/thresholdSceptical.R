
thresholdSceptical <- function(level, 
                               alternative = "one.sided", 
                               type = "nominal"){
    if (!(type %in% c("nominal", "liberal", "controlled", "golden")))
        stop('type must be either "nominal", "liberal", "controlled", or "golden"')
    if (!(alternative %in% c("one.sided", "two.sided")))
        stop('type must be either "one.sided" or "two.sided"')
    
    if(type == "nominal")
        res <- level
    if(type == "liberal")
        res <- pIntrinsic(p = level, alternative = "one.sided", type = "Held")
    if(type == "controlled"){
        pI2 <- pIntrinsic(p = 2*level^2, alternative = "one.sided", type = "Held")
        res <- pIntrinsic(p = pI2, alternative = "one.sided", type = "Held")
    }
    if(type == "golden"){
        res <- pIntrinsic(p = level, alternative = "one.sided", type = "Matthews")
    }
    if(alternative == "two.sided")
        res <- 2*res
    
    return(res)
}

