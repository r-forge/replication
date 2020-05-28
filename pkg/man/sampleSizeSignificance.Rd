\name{sampleSizeSignificance}
\alias{sampleSizeSignificance}
\title{
  Computes the required relative sample size to achieve significance with a certain power}
\description{
  The relative sample size to achieve significance is computed based on the z-value of the original study, the power, the design prior, and the relative heterogeneity.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizeSignificance(zo, power, level = 0.025, designPrior = "conditional", 
                       alternative = "one.sided", d = 0, shrinkage = 0)
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
    Significance level. Default is 0.025.
  }
  \item{designPrior}{
    Either \code{"conditional"}, \code{"predictive"}, or \code{"EB"}.  
    Defaults to \code{"conditional"}. If \code{"EB"}, the power is computed under a predictive distribution where the contribution of the original study is shrunken towards zero based on the evidence in the original study (with an empirical Bayes shrinkage estimator).
  }
  \item{alternative}{
    Either \code{"two.sided"}, \code{"one.sided"}, \code{"less"}, or \code{"greater"}.  
    Specifies direction of the alternative. 
    Defaults to \code{"one.sided"}, the same direction as the original estimate.  
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  Default is \code{0} (no heterogeneity).
  Is only taken into account when \code{designPrior = "predictive"} or \code{designPrior = "EB"}.
  } 
  \item{shrinkage}{
  A number in [0,1]. 
  Defaults to \code{0}. 
   Specifies how much the original effect estimate is shrunken towards zero (e.g. the effect is shrunken by a factor of 25\% for \code{shrinkage = 0.25}).
  Is only taken into account when \code{designPrior = "conditional"} or \code{designPrior = "predictive"}.
  }
}
\value{
  The relative sample size to achieve significance in the specified direction. 
  If larger than 1000 then NA is returned. 
}
\references{
  Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 183(2):431 - 448.
  \url{https://doi.org/10.1111/rssa.12493},
  
  Pawel, S., Held, L. (2020). Probabilistic forecasting of replication studies. PLoS ONE 15(4):e0231416. \url{https://doi.org/10.1371/journal.pone.0231416}
}
\author{
  Leonhard Held, 
  Samuel Pawel
}
\seealso{
  \code{\link{powerSignificance}}
}
\examples{
sampleSizeSignificance(zo = p2z(0.005), power = 0.8)
sampleSizeSignificance(zo = p2z(0.005, alternative = "greater"), power = 0.8)
sampleSizeSignificance(zo = p2z(0.005), power = 0.8, designPrior = "predictive")

sampleSizeSignificance(zo = 3, power = 0.8, designPrior = "predictive", 
                       shrinkage = 0.5, d = 0.25)
sampleSizeSignificance(zo = 3, power = 0.8, designPrior = "EB", 
                       d = 0.5)
                       
# required relative sample size for 0.8 power as function of original p-value
zo <- p2z(seq(0.0001, 0.05, 0.0001))
plot(z2p(zo), sampleSizeSignificance(zo = zo, designPrior = "conditional", power = 0.8),
     type = "l", ylim = c(0.5, 10), log = "y", lwd = 1.5, ylab = "Relative sample size",
     xlab = expression(italic(p)[o]), las = 1)
lines(z2p(zo), sampleSizeSignificance(zo = zo, designPrior = "predictive", power = 0.8), 
      lwd = 2, lty = 2)
lines(z2p(zo), sampleSizeSignificance(zo = zo, designPrior = "EB", power = 0.8),
      lwd = 1.5, lty = 3)
legend("topleft", legend = c("conditional", "predictive", "EB"), 
       title = "Design prior", lty = c(1, 2, 3), lwd = 1.5, bty = "n")
}
