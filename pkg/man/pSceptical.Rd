\name{pSceptical}
\alias{pSceptical}
\title{Computes the sceptical p-value}
\description{
  The sceptical p-value is computed based on the z-values of the
  original and the replication study and the corresponding variance ratio.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pSceptical(zo, zr, c, alternative = "one.sided", type="nominal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zo}{
    A vector of z-values from original studies.
  }
  \item{zr}{
    A vector of z-values from replication studies.
  }
  \item{c}{A vector of variance ratios of the original and replication
          effect estimates. This is usually the ratio of the sample
          size of the replication study to the sample size of the
          original study.}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to \code{"one.sided"}. If \code{"one.sided"} then the
    sceptical p-value is based on a one-sided assessment of replication
    success in the direction of the original effect estimate. If \code{"two.sided"} then the
    sceptical p-value is based on a two-sided assessment of replication
    success regardless of the direction of the original and replication
    effect estimate.
  }
\item{type}{
Type of recalibration. Default is \code{"nominal"}, so no recalibration as in 
  Held (2020). This ensures that the sceptical p-value is always larger
  than both the original and replication p-values. A recalibration is
  applied of \code{type} is \code{"liberal"}, \code{"controlled"} or
  \code{"golden"}.
}
}
\value{
The sceptical p-value. 
}
\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 183(2):431 - 448.
  \url{https://doi.org/10.1111/rssa.12493}
}
\author{
Leonhard Held
}
\seealso{
  \code{\link{sampleSizeReplicationSuccess}}, \code{\link{powerReplicationSuccess}}, \code{\link{thresholdSceptical}}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## two-sided p-values 0.01 and 0.02, relative sample size 2
pSceptical(zo = p2z(0.01), zr = p2z(0.02), c = 2, alternative = "one.sided")
## reverse the studies
pSceptical(zo = p2z(0.02), zr = p2z(0.01), c = 1/2, alternative = "one.sided")
## both p-values 0.01, relative sample size 2
pSceptical(zo = p2z(0.01), zr = p2z(0.01), c = 2, alternative = "two.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
