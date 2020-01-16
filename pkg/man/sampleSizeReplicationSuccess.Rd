\name{sampleSizeReplicationSuccess}
\alias{sampleSizeReplicationSuccess}
\title{Computes the required relative sample size to achieve replication success with a certain power}
\description{
  The relative sample size for replication success is computed based on the z-value of the original study, the power, and the design prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizeReplicationSuccess(zo, power, level = 0.025, 
                             designPrior = "conditional",
                             alternative = "one.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zo}{
  A vector of z-values from original studies. 
  }
  \item{power}{
  The power to achieve replication success.
  }
  \item{level}{
  Replication success level. Default is 0.025.
  }
  \item{designPrior}{
  Either \code{"conditional"}, \code{"predictive"}, or \code{"EB"}.  
  Defaults to \code{"conditional"}.
}
	\item{alternative}{
	Either \code{"one.sided"} or \code{"two.sided"}.  
	Defaults to \code{"one.sided"}. If \code{"one.sided"}, then sample size calculations are based on a one-sided assessment of replication success in the direction of the original effect estimate.
	 }
}
\value{
The relative sample size for replication success. If larger than 1000
\code{Inf} is returned. 
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
  \code{\link{pSceptical}}, \code{\link{powerReplicationSuccess}}, \code{\link{thresholdSceptical}}
}
\examples{
thresControlled <- thresholdSceptical(0.025, type = "controlled")

sampleSizeReplicationSuccess(zo = p2z(0.005), power = 0.8)
sampleSizeReplicationSuccess(zo = p2z(0.005), power = 0.8, 
                             designPrior = "predictive")
sampleSizeReplicationSuccess(zo = p2z(0.005), power = 0.8, 
                             level = thresControlled,
                             alternative = "one.sided")
sampleSizeReplicationSuccess(zo = p2z(0.005), power = 0.8, 
                             level = thresControlled,
                             alternative = "one.sided", 
                             designPrior = "predictive")
}
 