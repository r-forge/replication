\name{thresholdSceptical}
\alias{thresholdSceptical}
\title{Computes the threshold for the sceptical p-value}
\description{
  The threshold for the sceptical p-value is computed.}
\usage{
thresholdSceptical(level, alternative = "one.sided", type = "nominal")
}
\arguments{
\item{level}{
The required significance level. This can be a vector.
}
\item{alternative}{
either \code{"one.sided"} or \code{"two.sided"}.
Defaults to \code{"one.sided"}.  Specifies if the threshold is for
one-sided or two-sided sceptical p-values.
}
\item{type}{
Type of threshold. Default is \code{"nominal"} as in 
  Held (2020). This ensures that the Type-I error is always smaller than
  \code{level}^2. Significance
  of both the original and replication study at \code{level} is then a necessary but not
  sufficient requirement for replication success.  If \code{"type"} is
  \code{"liberal"} then significance of both studies is a sufficient
  requirement for replication success if the variance ratio is equal to 1.
  If \code{"type"} is \code{"controlled"} then the Type-I error is equal
  to \code{level}^2 if the variance ratio is equal to 1. 
}
}
\value{
The threshold for the sceptical p-value.
}
\references{
Held, L. (2020). A new standard for the analysis and design of replication studies (with discussion). 
  \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}. 
  \url{https://doi.org/10.1111/rssa.12493}

Held, L. (2019). The harmonic mean chi squared test to substantiate scientific findings. 
  \url{https://arxiv.org/pdf/1911.10633.pdf}. Preprint.
}
\author{
Leonhard Held
}
\examples{
thresholdSceptical(level = 0.025, type = "nominal")
thresholdSceptical(level = 0.025, type = "liberal")
thresholdSceptical(level = 0.025, type = "controlled")
}
