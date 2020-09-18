\name{sampleSizeRES}
\alias{sampleSizeRES}
\title{
  Computes the required relative sample size to achieve replication success based on the relative effect size
  }
\description{
  The relative sample size to achieve replication success is computed based on the z-value of the original study, the level for replication success and the specification of the minimum relative effect size.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSizeRES(zo, 
              d = 1,
              level = thresholdSceptical(level = 0.025, 
                                         alternative = "one.sided", 
                                         type = "golden"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zo}{
    A vector of z-values from original studies. 
  }
  \item{d}{
    The relative effect size (ratio of the effect estimate 
    from the replication study to the effect estimate 
    from the original study)
  }
  \item{level}{
    One-sided replication success level. Default is the golden level (0.0617).
  }
}
\value{
  The relative sample size to achieve replication success with the
  specified relative effect size (or larger) at the specified level. 
  If the relative effect size is smaller than the lower bound on the relative 
  effect size at the specified level, then NA is returned. 
  If the p-value from the original study (calculated from the corresponding z-value)
  is larger than the level, then NA is returned. 
}

\author{
  Leonhard Held, 
  Charlotte Micheloud
}
\seealso{
  \code{\link{sampleSizeSignificance}}, 
  \code{\link{sampleSizeReplicationSuccess}}
}
\examples{
thresNominal <- thresholdSceptical(level = 0.025, 
                                  alternative = "one.sided", 
                                  type = "nominal")
thresGolden <- thresholdSceptical(level = 0.025, 
                                  alternative = "one.sided", 
                                  type = "golden")
                                  
po <- c(0.001, 0.002, 0.01, 0.02, 0.025)
zo <- p2z(po, alternative = "one.sided")

sampleSizeRES(zo = zo, 
              d = 1, 
              level = thresGolden)
              
sampleSizeRES(zo = zo, 
              d = 1, 
              level = thresNominal)
              
              
sampleSizeRES(zo = zo, 
              d = 0.9, 
              level = thresGolden)
              
sampleSizeRES(zo = zo, 
              d = 0.9, 
              level = thresNominal)
              
              

}
