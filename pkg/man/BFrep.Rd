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
BFrep(zo, zr, c = 1)
}
\arguments{
  \item{zo}{
    A vector of z-values from original studies. 
  }
  \item{zr}{
    A vector of z-values from replication studies. 
  }
  \item{c}{
    A vector of ratios of the variance of the original and replication effect estimates. This is usually the ratio of the sample size of the
    replication study to the sample size of the original study. 
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
BFrep(zo = c(2, 4), zr = c(3, 3.5), c = c(1, 1.5))
  
}