\name{sampleSizePI}
\alias{sampleSizePI}
\title{
  Computes the required relative sample size for the prediction interval of the replication effect estimate not to include zero
}
\description{
  The relative sample size for the prediction interval of the replication effect estimate not to include zero is computed based on the z-value of the original study, the confidence level, the design prior, and the relative heterogeneity. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizePI(zo, conf.level = 0.95, designPrior = "predictive", d = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zo}{
    A vector of z-values from original studies. 
  }
  \item{conf.level}{
    The confidence level of the prediction intervals. Default is 0.95.
  }
  \item{designPrior}{
    Either \code{"conditional"}, \code{"predictive"}, or \code{"EB"}.  
    Defaults to \code{"predictive"}. If \code{"EB"}, the contribution of the original study to the predictive distribution is shrunken towards zero based on the evidence in the original study (with empirical Bayes).
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  Default is \code{0} (no heterogeneity).
  Is only taken into account when \code{designPrior = "predictive"} or \code{designPrior = "EB"}.
  } 
}
\value{
  The relative sample size for a prediction interval of the replication effect estimate
  to have a certain width (\code{method = "width"}) or not to include zero (\code{method = "fixzero"}). 
  If this is impossible at the specified level, NA is returned.
}
\references{
 Patil, P., Peng, R. D., Leek, J. T. (2016). What should researchers expect when they replicate studies? A statistical view of replicability in psychological science. \url{https://dx.doi.org/10.1177/1745691616646366}
  
  Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. \url{https://dx.doi.org/10.31234/osf.io/fhwb7}
}
\author{
  Samuel Pawel
}
\seealso{
  \code{\link{predictionInterval}}, \code{\link{sampleSizePIwidth}}
}
\examples{
sampleSizePI(zo = c(1, 2, 3), conf.level = 0.99)
sampleSizePI(zo = c(1, 2, 3), designPrior = "EB")

# required relative sample size for 0.95 PI not to include 0 as function of 
# z-value
zo <- seq(0, 5, 0.01)
plot(zo, sampleSizePI(zo = zo), type = "l", log = "y", lwd = 1.5, las = 2,
     ylab = expression(paste("Relative sample size ", n[r]/n[o])),
     xlab = expression(italic(z)[o]))
abline(h = 1, lty = 3)
lines(zo, sampleSizePI(zo = zo, designPrior = "EB"), lty = 2, lwd = 1.5)
legend("topright", title = "designPrior", lty = c(1, 2, 4), lwd = 1.5,
       legend = c("predictive", "EB"), bty = "n")
}
