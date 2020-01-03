\name{pBox}
\alias{pBox}
\title{Computes Box's tail probability}
\description{
  Box's tail probability is
  computed based on the test statistic of the
  original and the replication study, the corresponding variance ratio
  and the significance level. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pBox(po=NULL, pr=NULL, to=p2t(po, alternative=alternative),
tr=p2t(pr, alternative=alternative), c, sig.level=0.05,
alternative="two.sided")
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
  \item{sig.level}{Significance level. Default is 0.05.}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to
    \code{"two.sided"}.  Specifies if one-sided or two-sided Box's tail
    probability is computed. If \code{"one.sided"} then the p-values 
    given in \code{po} and \code{pr} are assumed to be
    one-sided.}
}
\value{
Box's tail probability. 
}
\references{
Box, G.E.P. (1980). Sampling and Bayes' inference in scientific
modelling and robustness (with discussion). \emph{Journal of the Royal Statistical
  Society, Series A}, \bold{143}, 383-430. 


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
\examples{
pBox(po=0.01, pr=0.02, c=2)
pBox(po=0.02, pr=0.01, c=1/2)
pBox(po=0.02, pr=0.01, c=1/2, alternative="one.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
