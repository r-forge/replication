\name{powerReplicationSuccess}
\alias{powerReplicationSuccess}
\title{Computes the power for replication success}
\description{
  The power for replication success is
  computed based on the result of the
  original study, the corresponding variance ratio and the design prior.
}
\usage{
powerReplicationSuccess(po=NULL, to=p2t(po, alternative=alternative),
c=1, level=0.05, designPrior="conditional", alternative="two.sided")
}
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. Is calculated from \code{po}, if necessary. 
  }
  \item{c}{A vector of variance ratios of the original and replication
    effect estimates. This is usually the ratio of the sample
    size of the replication study to the sample size of the
    original study.}
  \item{level}{Replication success level. Default is 0.05.}
  \item{designPrior}{either \code{"conditional"} or \code{"predictive"}.  Defaults to
    \code{"conditional"}.}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
          \code{"two.sided"}.  Specifies if the P-value is one-sided or
two-sided. If the P-value is one-sided, then power is based on a
one-sided assessment of replication success.}
}
\value{
  The power for replication success. 
}
\references{
Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 
  \url{https://doi.org/10.1111/rssa.12493}
}
\author{
  Leonhard Held
}
\seealso{
  \code{\link{sampleSizeReplicationSuccess}}, \code{\link{pSceptical}}, \code{\link{zSceptical}}
}
\examples{
powerReplicationSuccess(po=0.005, c=2)
powerReplicationSuccess(po=0.005, c=2, designPrior="predictive")
powerReplicationSuccess(po=0.005, c=2, alternative="one.sided")
powerReplicationSuccess(po=0.005, c=2, designPrior="predictive", alternative="one.sided")

powerReplicationSuccess(po=0.005, c=1/2)
powerReplicationSuccess(po=0.005, c=1/2, designPrior="predictive")
powerReplicationSuccess(po=0.005, c=1/2, alternative="one.sided")
powerReplicationSuccess(po=0.005, c=1/2, designPrior="predictive", alternative="one.sided")

}
