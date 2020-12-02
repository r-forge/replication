\name{scepticalCI}
\alias{scepticalCI}
\title{Sceptical confidence interval}
\description{
The sceptical confidence interval is computed based on the effect estimates of the original and replication studies and their respective standard errors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scepticalCI(thetao, thetar, se_thetao, se_thetar, alphas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thetao}{
   Effect estimate from original study.
  }
  \item{thetar}{
    Effect estimate from replication study.
  }
  
   \item{se_thetao}{
    Standard error of \code{thetao}
  }
  \item{se_thetar}{
   Standard error of \code{thetar}
  }
  \item{alphas}{Threshold for replication success}
}
\value{
The sceptical confidence interval. 
}
\references{
Held, L. (2020). The harmonic mean chi-squared test to substantiate scientific findings. 
  \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \bold{69}, 697-708.  
\url{https://doi.org/10.1111/rssc.12410}
}
\author{
Charlotte Micheloud
}
\seealso{
  \code{\link{pScepticalMu}}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% \examples{
% scepticalCI(thetao = 0.4, thetar = 0.2, se_thetao = 0.1, se_thetar = 0.1, 
% alphas = thresholdSceptical(level = 0.05, type = "controlled", alternative = "two.sided"))
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
