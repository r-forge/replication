\name{zSceptical}
\alias{zSceptical}
\title{Computes the test statistic of the sceptical P-value}
\description{
  The test statistic of the sceptical P-value is
  computed based on the test statistic of the
  original and the replication study and the corresponding variance ratio.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zSceptical(to, tr, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{to}{
    A vector of test statistics from original studies. 
  }
  \item{tr}{
    A vector of test statistics from replication studies. 
  }
  \item{c}{A vector of variance ratios of the original and replication
          effect estimates. This is usually the ratio of the sample
          size of the replication study to the sample size of the
          original study.}
}
\value{
The test statistic of the sceptical P-value. 
}
\references{
Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 
  \url{https://doi.org/10.1111/rssa.12493}
}
\author{
Leonhard Held
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\seealso{
  \code{\link{sampleSizeReplicationSuccess}}, \code{\link{powerReplicationSuccess}}, \code{\link{pSceptical}}
}
\examples{
zSceptical(to=2, tr=3, c=2)
zSceptical(to=3, tr=2, c=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
