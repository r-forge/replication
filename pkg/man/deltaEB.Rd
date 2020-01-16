\name{deltaEB}
\alias{deltaEB}
\title{Computes the weight for the original study with empirical Bayes}
\description{
The weight parameter delta for the original study is computed based on
the potential conflict with the replication study. A power prior
approach is combined with empirical Bayes. 
}
\usage{
deltaEB(zo, zr, c)
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
deltaEB(zo = p2z(0.01), zr = p2z(0.02), c = 2)
deltaEB(zo = p2z(0.02), zr =  p2z(0.01), c = 1/2)
deltaEB(zo = p2z(0.01, alternative = "one.sided"), 
        zr = p2z(0.01, alternative = "one.sided"), c = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
