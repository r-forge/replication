\name{levelEquivalent}
\alias{levelEquivalent}
\title{Computes the level for the sceptical p-value}

\description{
  The level for the sceptical p-value is computed based on the specified
  alternative and calibration type.
}

\usage{
levelEquivalent(dinf, level=.025, alternative = "one.sided")
}

\arguments{
\item{dinf}{
  The relative effect size bound below which replication success is
  impossible for borderline significant original studies. 
}
\item{level}{
  The required significance level. This can be a vector. Defaults to
  0.025. 
}
\item{alternative}{
  Either \code{"one.sided"} or \code{"two.sided"}.
  Defaults to \code{"one.sided"}. Specifies if the significance level
  is one-sided or two-sided. If the significance level is one-sided, 
  then a one-sided level for the sceptical p-value is computed.
}
}

\value{
The equivalent level where \code{dmin}=1, so at which replication success is
  impossible for borderline significant original studies and
  shrinkage of the replication effect estimate. 
}

\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics
    in Society)}, \bold{183}, 431-448.  
  \url{https://doi.org/10.1111/rssa.12493}

Held, L. (2020). The harmonic mean chi-squared test to substantiate scientific findings. 
  \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \bold{69}, 697-708.  
\url{https://doi.org/10.1111/rssc.12410}

Held, L., Micheloud, C. & Pawel, S. (2020). The assessment of replication success 
based on relative effect size. \url{https://arxiv.org/abs/2009.07782}
}

\seealso{
  \code{\link{pSceptical}}, \code{\link{T1EpSceptical}}
}

\author{
Leonhard Held
}

\examples{
levelEquivalent(dinf = 0.8, level = 0.025)
levelEquivalent(dinf = 0.8, level = 0.05, alternative="two.sided")
}
