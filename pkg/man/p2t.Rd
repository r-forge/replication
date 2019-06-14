\name{p2t}
\alias{p2t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transforms P-values to t-values
%%  ~~function to do ... ~~
}
\description{
    Transforms one- or two-sided p-value(s) to the corresponding t-value(s).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p2t(p, alternative="two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
   a vector of p-values
%%     ~~Describe \code{bf} here~~
}
   \item{alternative}{
    type of p-value. Either "two.sided" or "one.sided". 
    Defaults to "two.sided". Specifies if the P-values are one-sided or
two-sided. 
%%     ~~Describe \code{prior.odds} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A numeric vector of t-values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
p2t(p=c(0.005, 0.01, 0.05))
p2t(p=c(0.005, 0.01, 0.05), alternative="one.sided")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
