\name{sampleSizeSignificance}
\alias{sampleSizeSignificance}
\title{
  Computes the required relative sample size to achieve significance with a certain power}
\description{
  The relative sample size to achieve significance is computed based on the test statistic of the original study, the power, the design prior, and the realative heterogeneity.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  sampleSizeSignificance(po = NULL, to = p2t(po, alternative = alternative),
  power, level = 0.05, designPrior = "conditional", alternative = "two.sided",
  shrinkage = FALSE, d = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. 
    Is calculated from \code{po}, if necessary. 
  }
  \item{power}{
    The power to achieve replication success.
  }
  \item{level}{
    Significance level.
  }
  \item{designPrior}{
    Either \code{"conditional"} or \code{"predictive"}.  
    Defaults to \code{"conditional"}.
  }
  \item{alternative}{
    Either \code{"one.sided"} or \code{"two.sided"}.  
    Defaults to \code{"two.sided"}.  
    Specifies if the P-value is one-sided or two-sided. 
    If the P-value is one-sided, then sample size calculations are based on a one-sided assessment of significance.
  }
  \item{shrinkage}{
  Either \code{TRUE} or \code{FALSE}. 
  Defaults to \code{FALSE}. If \code{TRUE} and \code{designPrior="predictive"}, the probability of significance is computed under a predictive distribution where the original effect is shrunken towards zero based on the evidence in the original study (with an empirical Bayes shrinkage estimator).
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  }
}
\value{
  The relative sample size to achieve significance. 
  If larger than 1000 then Inf is returned. 
}
\references{
  Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  To appear in \emph{Journal of the Royal Statistical Society, Series A}. 
  \url{https://arxiv.org/abs/1811.10287}
}
\author{
  Leonhard Held, 
  Samuel Pawel (shrinkage and heterogeneity options)
}
\examples{
sampleSizeSignificance(po = 0.005, power = 0.8)
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive")
sampleSizeSignificance(po = 0.005, power = 0.8, alternative = "one.sided")
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive", 
                       alternative = "one.sided", shrinkage = FALSE, d = 1)
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive", 
                       shrinkage = TRUE, d = 1)

}
