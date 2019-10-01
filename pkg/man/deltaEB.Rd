\name{deltaEB}
\alias{deltaEB}
\title{Computes the weight for the original study with empirical Bayes}
\description{
The weight parameter delta for the original study is computed based on
the potential conflict with the replication study. A power prior
approach is combined with empirical Bayes. 
}
\usage{
deltaEB(po=NULL, pr=NULL, to=p2t(po, alternative=alternative),
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
    are one-sided or two-sided. }
}
\value{
The weight delta. 
}
\references{
  Gravestock, I. and Held, L. (2017).
  Adaptive power priors with empirical Bayes for clinical trials
  \emph{Pharmaceutical Statistics}, 16(5), 349-360, 2017.
\url{https://dx.doi.org/10.1002/pst.1814}
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
deltaEB(po=0.01, pr=0.02, c=2)
deltaEB(po=0.02, pr=0.01, c=1/2)
deltaEB(po=0.01, pr=0.01, c=2, alternative="one.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
