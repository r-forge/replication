\name{CI.to.t}
\alias{CI.to.t}
\title{
    Transforms confidence intervals to test statistics. 
%%  ~~function to do ... ~~
}
\description{
    Transforms confidence intervals to the corresponding test statistics.
}
\usage{
CI.to.t(lower, upper, conf.level=0.95, ratio=FALSE)
}
\arguments{
  \item{lower}{
   a vector of lower confidence interval bounds
}
  \item{upper}{
   a vector of upper confidence interval bounds
}
    \item{conf.level}{
      The confidence level of the confidence intervals. Default is 0.95.
    }
   \item{ratio}{
     Indicates whether the confidence interval is for a ratio, e.g. an
     odds ratio, relative risk or hazard ratio. If \code{TRUE} the
     z-value of the log ratio is computed. 
    Defaults to \code{FALSE}.
}
}
\value{A numeric vector of t-values.
}
\examples{
CI.to.t(lower=1, upper=3)
CI.to.t(lower=1, upper=3, ratio=TRUE)
CI.to.t(lower=1, upper=3, conf.level=0.9)
}

