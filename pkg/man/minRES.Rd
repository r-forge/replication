\name{minRES}
\alias{minRES}
\title{
  Computes the minimal relative effect size to achieve replication success
  }
\description{
  The minimal relative effect size to achieve replication success is computed based on the z-value of the original study, the level for replication success and the relative sample size.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
minRES(zo, 
       c,
       level = thresholdSceptical(level = 0.025, 
                                         alternative = "one.sided", 
                                         type = "golden"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zo}{
    A vector of z-values from original studies. 
  }
  \item{c}{
    The ratio of the sample size of the replication study to the sample size of the original study
  }
  \item{level}{
    One-sided replication success level. Default is the golden level (0.0617).
  }
}
\value{
  The relative effect size to achieve replication success with the specified relative sample size and at the specified level. 
  If the p-value from the original study (calculated from the corresponding z-value)
  is smaller than the level, then NA is returned. 
}

\author{
  Leonhard Held, 
  Charlotte Micheloud
}
\seealso{
  \code{\link{sampleSizeRES}}

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

minRES(zo = zo, 
       c = 1, 
       level = thresGolden)
              
minRES(zo = zo, 
       c = 1, 
       level = thresNominal)
              
              
minRES(zo = zo, 
       c = 2, 
       level = thresGolden)
       
minRES(zo = zo, 
       c = 2, 
       level = thresNominal)
              
              

}
