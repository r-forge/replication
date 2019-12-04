\name{sampleSizeSignificance}
\alias{sampleSizeSignificance}
\title{
  Computes the required relative sample size to achieve significance with a certain power}
\description{
  The relative sample size to achieve significance is computed based on the test statistic of the original study, the power, the design prior, and the relative heterogeneity.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  sampleSizeSignificance(po = NULL, to = p2t(po, alternative = alternative),
  power, level = 0.05, designPrior = "conditional", alternative = "two.sided",
  d=0, shrinkage=1)
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
    Either \code{"conditional"}, \code{"predictive"}, or \code{"EB"}.  
    Defaults to \code{"conditional"}. If \code{"EB"}, the power is computed under a predictive distribution where the contribution of the original study is shrunken towards zero based on the evidence in the original study (with an empirical Bayes shrinkage estimator).
  }
  \item{alternative}{
    Either \code{"less"}, \code{"greater"}, or \code{"two.sided"}.  
    Defaults to \code{"two.sided"}.  
    Specifies if the alternative is one-sided or two-sided. 
    If the alternative is one-sided, then sample size calculations are based on a one-sided assessment of significance.
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  Is only taken into account when \code{designPrior = "predictive"} or \code{designPrior = "EB"}.
  } 
  \item{shrinkage}{
  A number in [0,1]. 
  Defaults to \code{1}. 
  Specifies how much the original effect estimate is shrunken towards zero (e.g. the effect is halfed for \code{shrinkage = 0.5}). 
  Is only taken into account when \code{designPrior = "conditional"} or \code{designPrior = "predictive"}.
  }
}
\value{
  The relative sample size to achieve significance in the specified direction. 
  If larger than 1000 then Inf is returned. 
}
\references{
  Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  To appear in \emph{Journal of the Royal Statistical Society, Series A}. 
  \url{https://arxiv.org/abs/1811.10287}
  
  Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. \url{https://dx.doi.org/10.31234/osf.io/fhwb7}
}
\author{
  Leonhard Held, 
  Samuel Pawel (shrinkage and heterogeneity options)
}
\examples{
sampleSizeSignificance(po = 0.005, power = 0.8)
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive")
sampleSizeSignificance(po = 0.005, power = 0.8, alternative = "greater")
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive", 
                       alternative = "greater", d = 1)
sampleSizeSignificance(po = 0.005, power = 0.8, designPrior = "predictive", 
                       shrinkage = 0.5, d = 0.25)
sampleSizeSignificance(to = 3, power = 0.8, designPrior = "EB", 
                       d = 0.5)
}
