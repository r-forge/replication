\name{alphaIntrinsic}
\alias{alphaIntrinsic}
\title{Computes the threshold for intrinsic credibility}
\description{
The P-value threshold for intrinsic credibility is computed. 
}
\usage{
alphaIntrinsic(alpha, alternative="two.sided", type="Held")
}
\arguments{
  \item{alpha}{The intrinsic credibility level. This can be a vector.}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.
    Defaults to \code{"two.sided"}.  Specifies if the threshold is for
    one-sided or two-sided P-values.}
  \item{type}{Type of intrinsic p-value threshold. Default is \code{"Held"} as in
Held (2019). The other option is \code{"Matthews"} as in Matthews
(2018).}

}
\value{
The threshold for intrinsic credibility. 
}
\references{Held, L. (2019). The assessment of intrinsic
credibility and a new argument for p < 0.005.  \emph{Royal Society Open
  Science}, 6:181534. \url{https://dx.doi.org/10.1098/rsos.181534}

Matthews, R. A. J. (2018).  Beyond 'significance': principles and practice of the analysis
of credibility. \emph{Royal Society Open Science}, 5:171047. \url{https://doi.org/10.1098/rsos.171047}
}
\author{
Leonhard Held
}
\examples{
alphaIntrinsic(alpha=c(0.05, 0.1))
alphaIntrinsic(alpha=c(0.05, 0.1), alternative="one.sided")
}
