\name{BFrep}
\alias{BFrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replication Bayes factor}
%%  ~~function to do ... ~~
\description{
  Computes the replication Bayes factor for the normal-normal model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BFrep(po = NULL, to = p2t(po, alternative = alternative), 
      pr = NULL, tr = p2t(pr, alternative = alternative), 
      c = 1, alternative = "two.sided")
}
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. 
    Is calculated from \code{po}, if necessary. 
  }
  \item{pr}{
    A vector of P-values from replication studies, default is \code{NULL}
  }
  \item{tr}{
    A vector of test statistics from replication studies. 
    Is calculated from \code{pr}, if necessary. 
  }
  \item{c}{
    A vector of ratios of the variance of the original and replication effect
    estimates. This is usually the ratio of the sample size of the
    replication study to the sample size of the original study. 
  }
  \item{alternative}{
    Alternative hypothesis of P-value. Default is \code{"two.sided"}. Is
    used to compute test-statistic from P-values.
  }
}
\details{
This function computes the Bayes factor comparing the null hypothesis;
the effect is zero, to the alternative; the effect is distributed according 
to the posterior distribution of the effect conditional on the original study
(assuming a normal model with fixed variances and a flat prior for the effect). 
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
The replication Bayes factor.
}
\references{
  
  Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result of a replication attempt. Journal of Experimental Psychology: General, 143:1457 - 1475. \url{https://dx.doi.org/10.1037/a0036731}

  Ly, A., Etz, A., Marsman, M., and Wagenmakers, E. J. (2018). Replication Bayes factors from evidence updating. Behavior Research Methods. \url{https://dx.doi.org/10.3758/s13428-018-1092-x}
  
}
\author{
  Samuel Pawel
}
\examples{
BFrep(to = c(2, 4), tr = c(3, 3.5), c = c(1, 1.5))
  
}