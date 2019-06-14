\name{sampleSizeSignificance}
\alias{sampleSizeSignificance}
\title{Computes the required relative sample size to achieve significance
  with a certain power}
\description{
  The relative sample size to achieve significance is
  computed based on the test statistic of the
  original study, the power and the design prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizeSignificance(po=NULL, to=p2t(po, alternative=alternative),
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
  \item{level}{Significance level.}
  \item{designPrior}{either \code{"conditional"} or \code{"predictive"}.  Defaults to
          \code{"conditional"}.}
    \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
      \code{"two.sided"}.  Specifies if the P-value is one-sided or
      two-sided. If the P-value is one-sided, then sample size
      calculations are based on a
      one-sided assessment of significance.}
      }
\value{
The relative sample size to achieve significance. If larger that 1000
then Inf is returned. 
}
\references{
    Held, L. (2019). A new standard for the analysis and design of
    replication studies (with discussion). To appear in
    \emph{Journal of the Royal Statistical Society, Series A}. 
    \url{https://arxiv.org/abs/1811.10287}
}
\author{
Leonhard Held
}
\examples{
sampleSizeSignificance(po=0.005, power=0.8)
sampleSizeSignificance(po=0.005, power=0.8, designPrior="predictive")
sampleSizeSignificance(po=0.005, power=0.8, alternative="one.sided")
sampleSizeSignificance(po=0.005, power=0.8, designPrior="predictive", alternative="one.sided")

}
