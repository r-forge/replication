\name{CI2p}
\alias{CI2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transforms confidence intervals to p-values.
%%  ~~function to do ... ~~
}
\description{
    Transforms confidence intervals to the corresponding p-values.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CI2p(lower, upper, conf.level=0.95, ratio=FALSE, alternative="two.sided")
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
    Defaults to \code{FALSE}.
}
   \item{alternative}{
    type of p-value. Either "two.sided" or "one.sided". 
    Defaults to "two.sided". Specifies if the P-values are one-sided or
    two-sided. 
}
}
\value{A numeric vector of p-values
}
\examples{
CI2p(1, 3)
CI2p(1, 3, alternative="one.sided")
}

