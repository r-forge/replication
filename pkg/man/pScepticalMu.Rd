\name{pScepticalMu}
\alias{pScepticalMu}
\title{Sceptical p-value for a general point null hypothesis}
\description{
The sceptical p-value is computed based on the effect estimates of the original and replication studies and their respective standard errors for the general point null hypothesis \eqn{H_0: \theta = \mu}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pScepticalMu(thetao, thetar, se_thetao, se_thetar, mu, alternative = "two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thetao}{
    A vector of effect estimates from original studies.
  }
  \item{thetar}{
    A vector of effect estimates from replication studies.
  }
  
   \item{se_thetao}{
    Standard errors of \code{thetao}
  }
  \item{se_thetar}{
   Standard errors of \code{thetar}
  }
  \item{mu}{A vector of values for the null hypothesis}
  \item{alternative}{either \code{"one.sided"} or \code{"two.sided"}.  Defaults to \code{"two.sided"}. If \code{"one.sided"} then the
    sceptical p-value is based on a one-sided assessment of replication
    success in the direction of the original effect estimate. If \code{"two.sided"} then the
    sceptical p-value is based on a two-sided assessment of replication
    success regardless of the direction of the original and replication effect estimate.
  }
}
\value{
The sceptical p-value. 
}
\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 183(2):431 - 448.
  \url{https://doi.org/10.1111/rssa.12493}
  
  Held, L. (2020). The harmonic mean X2-test to substantiate scientific findings. 
  \emph{Journal of the Royal Statistical Society Series C}. To appear.
}
\author{
Charlotte Micheloud
}
\seealso{
  \code{\link{pSceptical}}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pScepticalMu(thetao = 0.4, thetar = 0.2, se_thetao = 0.1, se_thetar = 0.1, 
             mu = c(0, 0.2, 0.4, 0.6), alternative = "two.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
