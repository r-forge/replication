\name{powerReplicationSuccess}
\alias{powerReplicationSuccess}
\title{Computes the power for replication success}
\description{
  The power for replication success is computed based on the result of the original study, the corresponding variance ratio and the design prior. 
}
\usage{
powerReplicationSuccess(zo, c = 1, level = 0.025, 
                        designPrior = "conditional", 
                        alternative = "one.sided")
}
\arguments{
  \item{zo}{
    A vector of z-values from original studies. 
  }
  \item{c}{
    A vector of variance ratios of the original and replication
    effect estimates. This is usually the ratio of the sample
    size of the replication study to the sample size of the
    original study.
    }
  \item{level}{
    Replication success level. Default is 0.025.
    }
  \item{designPrior}{
    Either \code{"conditional"} or \code{"predictive"}.  
    Defaults to \code{"conditional"}. 
}
  \item{alternative}{
  Either \code{"one.sided"} or \code{"two.sided"}.  Defaults to \code{"one.sided"}. Specifies if the p-value is one-sided or two-sided. If \code{"one.sided"}, then power calculations are based on a one-sided assessment of replication success in the direction of the original effect estimate.
}
}
\value{
  The power for replication success. 
}
\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 183(2):431 - 448.
  \url{https://doi.org/10.1111/rssa.12493}
}
\author{
  Leonhard Held
}
\seealso{
  \code{\link{sampleSizeReplicationSuccess}}, \code{\link{pSceptical}}, \code{\link{thresholdSceptical}}
}
\examples{
thresControlled <- thresholdSceptical(0.025, type = "controlled")

powerReplicationSuccess(zo = p2z(0.005), c = 2)
powerReplicationSuccess(zo = p2z(0.005), c = 2, designPrior = "predictive")
powerReplicationSuccess(zo = p2z(0.005), c = 2, level = thresControlled)
powerReplicationSuccess(zo = p2z(0.005), c = 2, level = thresControlled, 
                        designPrior = "predictive")

powerReplicationSuccess(zo = p2z(0.005), c = 1/2)
powerReplicationSuccess(zo = p2z(0.005), c = 1/2, designPrior = "predictive")
powerReplicationSuccess(zo = p2z(0.005), c = 1/2, level = thresControlled)
powerReplicationSuccess(zo = p2z(0.005), c = 1/2, level = thresControlled, 
                        designPrior = "predictive")

}
