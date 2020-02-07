\name{powerSignificanceInterim}
\alias{powerSignificanceInterim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power of an adaptative replication study}
%%  ~~function to do ... ~~
\description{
Computes the probability that a replication study yields a significant effect estimate in the specified direction at the final analysis given the data so far (accumulated up to the interim analysis).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
powerSignificanceInterim(zo, zi, c = 1, f = 1/2,
                         level = 0.025, 
                         designPrior = "conditional",
                         analysisPrior = "flat", 
                         alternative="greater", 
                         shrinkage = 0)
}
\arguments{

  \item{zo}{
    A vector of z-values from original studies. 
  }
   \item{zi}{
    A vector of z-values from interim analyses of replication studies. 
  }
  \item{c}{
    The ratio of the sample size of the
    replication study up to the interim analysis to the sample size of the original study. 
  }
  \item{f}{
  Fraction of the replication study already completed.
  }
  \item{level}{
    The significance level.
  }
  \item{designPrior}{
    Either \code{"conditional"}, \code{"predictive"} or \code{"flat"}.  
    Defaults to \code{"conditional"}.
  }
    \item{analysisPrior}{
    Either \code{"flat"} or \code{"original"}.  
    Defaults to \code{"flat"}.
  }
  \item{alternative}{
    Either \code{"less"}, \code{"greater"}, or \code{"two.sided"}.  
    Defaults to \code{"greater"}.  
    Specifies direction of the alternative. 
    }
  \item{shrinkage}{
  A number in [0,1]. 
  Defaults to \code{0}. 
  Specifies how much the original effect estimate is shrunken towards zero (e.g. the effect is shrunken by a factor of 25\% for \code{shrinkage = 0.25}). 
  }
}
  \details{This is an extension of \code{PowerSignificance()} and adapts the `interim power' from section 6.6.3 of Spiegelhalter et al. (2004) to the setting of replication studies.
    
    %%  ~~ If necessary, more details than the description above ~~
  }

  \references{
  Spiegelhalter, D. J., Abrams, K. R., and Myles, J. P. (2004). Bayesian Approaches to Clinical Trials    and Health-Care
  Evaluation, volume 13. John Wiley & Sons
   
  }
  \author{
    Charlotte Micheloud
  }
\seealso{
  \code{\link{sampleSizeSignificance}}, \code{\link{powerSignificance}}
}
\examples{
powerSignificanceInterim(zo = 2, zi = 2, c = 1, f = 1/2, 
                         designPrior = "conditional", 
                         analysisPrior = "flat")

powerSignificanceInterim(zo = 2, zi = 2, c = 1, f = 1/2, 
                         designPrior = "predictive", 
                         analysisPrior = "flat")

powerSignificanceInterim(zo = 2, zi = 2, c = 1, f = 1/2,
                        designPrior = "flat", 
                        analysisPrior = "flat")

powerSignificanceInterim(zo = 2, zi = -2, c = 1, f = 1/2, 
                         designPrior = "conditional", 
                         analysisPrior = "flat")

powerSignificanceInterim(zo = 2, zi = 2, c = 1, f = 1/2, 
                         designPrior = "conditional", 
                         analysisPrior = "flat", 
                         shrinkage = 0.25)

}