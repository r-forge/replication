\name{pSceptical}
\alias{pSceptical}
\title{Computes the sceptical P-value}
\description{
  The sceptical P-value is
  computed based on the test statistic of the
  original and the replication study and the corresponding variance ratio.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pSceptical(po=NULL, pr=NULL, to=p2t(po, alternative=alternative),
tr=p2t(pr, alternative=alternative), c, alternative="two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{pr}{
    A vector of P-values from replication studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. Is calculated from \code{po}, if necessary. 
  }
  \item{tr}{
    A vector of test statistics from replication studies. Is calculated from \code{po}, if necessary. 
  }
  \item{c}{A vector of variance ratios of the original and replication
          effect estimates. This is usually the ratio of the sample
          size of the replication study to the sample size of the
          original study.}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
    \code{"two.sided"}.  Specifies if the P-values
    are one-sided or two-sided. If the P-values are one-sided, then the
    sceptical P-value is based on a one-sided assessment of replication success.}
}
\value{
The sceptical P-value. 
}
\references{
Held, L. (2019). A new standard for the analysis and design of
replication studies (with discussion). To appear in
\emph{Journal of the Royal Statistical Society, Series A}. 
\url{https://arxiv.org/abs/1811.10287}
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
\examples{
pSceptical(po=0.01, pr=0.02, c=2)
pSceptical(po=0.02, pr=0.01, c=1/2)
pSceptical(po=0.01, pr=0.01, c=2, alternative="one.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
