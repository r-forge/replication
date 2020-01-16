\name{sampleSizePIwidth}
\alias{sampleSizePIwidth}
\title{
  Computes the required relative sample size for the prediction interval of the replication effect estimate to have a desired width}
\description{
  The relative sample size for the prediction interval of the replication effect estimate to have a certain width is computed. 
  %based on the test statistic of the original study, the confidence level, the design prior, and the relative heterogeneity. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizePIwidth(zo, w, conf.level = 0.95, 
                  designPrior = "predictive", d = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{zo}{
    A vector of z-values from original studies. 
    Only required if \code{designPrior = "EB"}.
  }
  \item{w}{
    Desired relative width of prediction interval (relative to width of
    confidence interval from original estimate).
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
  to have a certain width. 
  If this is impossible at the specified level, NA is returned.
}
\author{
  Samuel Pawel
}
\references{
 Patil, P., Peng, R. D., Leek, J. T. (2016). What should researchers expect when they replicate studies? A statistical view of replicability in psychological science. \url{https://dx.doi.org/10.1177/1745691616646366}
  
  Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. \url{https://dx.doi.org/10.31234/osf.io/fhwb7}
}
\seealso{
  \code{\link{predictionInterval}}, \code{\link{sampleSizePI}}
}
\examples{
sampleSizePIwidth(w = c(1.1, 1.5, 2), d = 0.5)
sampleSizePIwidth(zo = c(1, 2, 3), w = 1.5, designPrior = "EB")

w <- seq(0.5, 3, 0.01)
plot(w, sampleSizePIwidth(zo = NULL, w = w), type = "l", log = "y", lwd = 1.5, las = 1,
     ylab = expression(paste("Relative sample size ", n[r]/n[o])),
     xlab = expression(paste("Relative width ", tilde(italic(w)))))
abline(h = 1, lty = 3)
lines(w, sampleSizePIwidth(zo = 2, w = w, designPrior = "EB"),
      lty = 2, lwd = 1.5)
legend("topright", title = "designPrior", lty = c(1, 2, 4), lwd = 1.5, bty = "n",
       legend = c("predictive", expression(paste("EB, ", italic(z)[o] == 2))))
}
