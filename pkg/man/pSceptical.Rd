\name{pSceptical}
\alias{pSceptical}
\title{Computes the sceptical p-value}
\description{
  The sceptical p-value is computed based on the z-values of the
  original and the replication study and the corresponding variance ratio.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pSceptical(zo, zr, c, alternative = "two.sided")
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
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to \code{"two.sided"}. If \code{"one.sided"} then the
    sceptical p-value is based on a one-sided assessment of replication success in the direction of the original effect estimate.}
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
pSceptical(zo = p2z(0.01), zr = p2z(0.02), c = 2)
pSceptical(zo = p2z(0.02), zr = p2z(0.01), c = 1/2)
pSceptical(zo = p2z(0.01), zr = p2z(0.01), c = 2, alternative = "one.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
