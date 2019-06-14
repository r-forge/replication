\name{CI2SE}
\alias{CI2SE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transforms confidence intervals to standard errors.
%%  ~~function to do ... ~~
}
\description{
    Transforms confidence intervals to the corresponding standard errors. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CI2SE(lower, upper, conf.level=0.95, ratio=FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
     odds ratio, relative risk or hazard ratio.
     If \code{TRUE} the standard error of the log ratio is computed. 
    Defaults to \code{FALSE}.
}
}
\value{A numeric vector of standard errors. 
}
\examples{
CI2SE(lower=1, upper=3)
CI2SE(lower=1, upper=3, ratio=TRUE)
CI2SE(lower=1, upper=3, conf.level=0.9)
}

