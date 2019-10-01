
library(ReplicationSuccess)
p <- seq(0.001, 0.075, 0.001)
## powerSignificance
pS <- powerSignificance(po=p, c=2)
pS2 <- powerSignificance(po=p, c=1)
pS3 <- powerSignificance(po=p, c=1/2)

## pReplicate
pR <- pReplicate(p=p, c=2)
pR2 <- pReplicate(p=p, c=1)
pR3 <- pReplicate(p=p, c=1/2)

## powerReplicationSuccess
pRS <- powerReplicationSuccess(po=p, c=2)
pRS2 <- powerReplicationSuccess(po=p, c=1)
pRS3 <- powerReplicationSuccess(po=p, c=1/2)



y <- cbind(pS,pS2,pS3,pR,pR2,pR3,pRS,pRS2,pRS3)

if(1){
    par(las=1, mfrow=c(1,1))
    matplot(p, y, type="l", ylim=c(0,1), lwd=2,
            lty=c(1,1,1,5,5,5,6,6,6), col=c(1,2,3,1,2,3),
            xlab="p-value", ylab="Power")
    where <- seq(0, 1, .1)
    for(i in seq_len(length(where)))
        abline(col="grey", lty=2, h=where[i])
    abline(col="grey", lty=2, v=0.05)
    axis(1, at=0.05,as.character(0.05)) 
    legend("bottomleft", legend=c("c=2", "c=1", "c=1/2"), lwd=2,
           col=c(1,2,3), lty=1, bg="white")
    legend("bottomright",
           legend=c("powerSignificance",
                    "pReplicate",
                    "powerReplicationSuccess"),
           lwd=2, col=c(1), lty=c(1,2,6), bg="white")
}



## powerSignificance one-sided
pS <- powerSignificance(po=p, c=2, alternative="one.sided")
pS2 <- powerSignificance(po=p, c=1, alternative="one.sided")
pS3 <- powerSignificance(po=p, c=1/2, alternative="one.sided")

## pReplicate one-sided
pR <- pReplicate(p=p, c=2, alternative="one.sided")
pR2 <- pReplicate(p=p, c=1, alternative="one.sided")
pR3 <- pReplicate(p=p, c=1/2, alternative="one.sided")

## powerReplicationSuccess one-sided
pRS <- powerReplicationSuccess(po=p, c=2, alternative="one.sided")
pRS2 <- powerReplicationSuccess(po=p, c=1, alternative="one.sided")
pRS3 <- powerReplicationSuccess(po=p, c=1/2, alternative="one.sided")



y <- cbind(pS,pS2,pS3,pR,pR2,pR3,pRS,pRS2,pRS3)

if(1){
    par(las=1)
    matplot(p, y, type="l", ylim=c(0,1), lwd=2,
            lty=c(1,1,1,5,5,5,6,6,6),
            col=c(1,2,3,1,2,3), xlab="p-value",
            ylab="Power")
    where <- seq(0, 1, .1)
    for(i in seq_len(length(where)))
        abline(col="grey", lty=2, h=where[i])
    abline(col="grey", lty=2, v=0.05)
    axis(1, at=0.05,as.character(0.05)) 
    legend("bottomleft", legend=c("c=2", "c=1", "c=1/2"),
           lwd=2, col=c(1,2,3), lty=1, bg="white")
    legend("bottomright",
           legend=c("powerSignificance",
                    "pReplicate",
                    "powerReplicationSuccess"),
           lwd=2, col=c(1), lty=c(1,2,6), bg="white")
}
