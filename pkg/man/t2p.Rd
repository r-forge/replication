\name{t2p}
\alias{t2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transforms t-values to p-values
%%  ~~function to do ... ~~
}
\description{
    Transforms the t-value(s) to the corresponding one- or two-sided p-value(s).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
t2p(t, alternative="two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
   a vector of t-values
%%     ~~Describe \code{bf} here~~
}
   \item{alternative}{Underlying
    Specifies direction of the alternative of p-value. 
    Either \code{"less"}, \code{"greater"}, \code{"one.sided"}, \code{"two.sided"}.  
    Defaults to \code{"two.sided"}.  
%%     ~~Describe \code{prior.odds} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A numeric vector of p-values
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
t2p(t=c(1,2,5))
t2p(t=c(1,2,5), alternative="less")
t2p(t=c(1,2,5), alternative="greater")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
