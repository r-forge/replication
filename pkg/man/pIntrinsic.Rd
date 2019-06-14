\name{pIntrinsic}
\alias{pIntrinsic}
\title{Computes the P-value for intrinsic credibility}
\description{
The P-value for intrinsic credibility is computed. 
}
\usage{
pIntrinsic(p=t2p(t, alternative=alternative), t=NULL,
alternative="two.sided", type="Held")
}
\arguments{
  \item{p}{a vector of P-values, is calculated from \code{z}, if necessary.}
  \item{t}{a vector of t-values, default is \code{NULL}}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
          \code{"two.sided"}.  Specifies if the P-value is one-sided or
two-sided. If the P-value is one-sided, then a one-sided P-value for
intrinsic credibility is computed.}
\item{type}{Type of intrinsic p-value. Default is \code{"Held"} as in
Held (2019). The other option is \code{"Matthews"} as in Matthews
(2018).}
}
\value{
The P-value for intrinsic credibility. 
}
\references{Held, L. (2019). The assessment of intrinsic
credibility and a new argument for p < 0.005. \emph{Royal Society Open
Science}, 6:181534. \url{https://dx.doi.org/10.1098/rsos.181534}

Matthews, R. A. J. (2018).  Beyond 'significance': principles and
practice of the analysis of credibility.
\emph{Royal Society Open Science}, 5:171047. \url{https://doi.org/10.1098/rsos.171047}
}
\author{
Leonhard Held
}
\examples{
p <- c(0.005, 0.01, 0.05)
t <- c(2,3,4)
pIntrinsic(p)
pIntrinsic(p, type="Matthews")
pIntrinsic(t=t)
pIntrinsic(t=t, type="Matthews")
pIntrinsic(p, alternative="one.sided")
pIntrinsic(p, alternative="one.sided", type="Matthews")
}
