\name{pScepticalMu}
\alias{pScepticalMu}
\title{Sceptical p-value for a general point null hypothesis}

\description{
The sceptical p-value is computed based on the effect estimates of the original 
and replication studies and their respective standard errors for the general 
point null hypothesis \eqn{H_0: \theta = \mu}. If specified, a recalibration is
applied.
}

\usage{
pScepticalMu(thetao, thetar, se_thetao, se_thetar, mu, 
             alternative = "two.sided", type = "golden")
}


\arguments{
\item{thetao}{
    A vector of effect estimates from original studies.
}
\item{thetar}{
    A vector of effect estimates from replication studies.
}

\item{se_thetao}{
    Standard errors of \code{thetao}
}
\item{se_thetar}{
    Standard errors of \code{thetar}
}
\item{mu}{
    A vector of values for the null hypothesis
}
\item{alternative}{
    Either \code{"one.sided"} or \code{"two.sided"}. Defaults to 
    \code{"two.sided"}. If \code{"one.sided"} then the
    sceptical p-value is based on a one-sided assessment of replication
    success in the direction of the original effect estimate. If \code{"two.sided"} then the
    sceptical p-value is based on a two-sided assessment of replication
    success regardless of the direction of the original and replication effect estimate.
  }
\item{type}{
  Type of recalibration. Can be either \code{"nominal"}, 
  \code{"liberal"}, \code{"controlled"} or \code{"golden"}.
  Defaults to \code{"golden"}.
  Setting  \code{type} to \code{"nominal"} corresponds to no recalibration 
  as in Held et al. (2020). A recalibration is
  applied if \code{type} is \code{"liberal"}, \code{"controlled"} or
  \code{"golden"} and the sceptical p-value can then be interpreted on 
  the same scale as an ordinary p-value (e.g. a one-sided sceptical p-value
  can be thresholded at the conventional 0.025 level). 
  See \code{?levelSceptical} for details about recalibration types. 
}
}

\value{
The sceptical p-value. 
}

\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics
    in Society)}, \bold{183}, 431-448.  
  \url{https://doi.org/10.1111/rssa.12493}
  
Held, L. (2020). The harmonic mean chi-squared test to substantiate scientific findings. 
  \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \bold{69}, 697-708.  
\url{https://doi.org/10.1111/rssc.12410}
}

\author{
Charlotte Micheloud
}

\seealso{
  \code{\link{pSceptical}}, \code{\link{levelSceptical}}
}