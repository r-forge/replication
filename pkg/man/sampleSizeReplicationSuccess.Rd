\name{sampleSizeReplicationSuccess}
\alias{sampleSizeReplicationSuccess}
\title{Computes the required relative sample size to achieve replication success
  with a certain power}
\description{
  The relative sample size for replication success is
  computed based on the test statistic of the
  original study, the power and the design prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizeReplicationSuccess(po=NULL, to=p2t(po, alternative=alternative),
power, level=0.05, designPrior="conditional", alternative="two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. Is calculated from \code{po}, if necessary. 
  }
  \item{power}{The power to achieve replication success.}
  \item{level}{Replication success level.}
  \item{designPrior}{either \code{"conditional"} or \code{"predictive"}.  Defaults to
          \code{"conditional"}.}
	\item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
          \code{"two.sided"}.  Specifies if the P-value is one-sided or
	  two-sided. If the P-value is one-sided, then sample size
	  calculations are based on a
	  one-sided assessment of replication success.}
      }
\value{
The relative sample size for replication success. If larger that 1000
then \code{Inf} is returned. 
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
  \code{\link{pSceptical}}, \code{\link{powerReplicationSuccess}}, \code{\link{zSceptical}}
}
\examples{
sampleSizeReplicationSuccess(po=0.005, power=0.8)
sampleSizeReplicationSuccess(po=0.005, power=0.8, designPrior="predictive")
sampleSizeReplicationSuccess(po=0.005, power=0.8, alternative="one.sided")
sampleSizeReplicationSuccess(po=0.005, power=0.8, designPrior="predictive", alternative="one.sided")

}
