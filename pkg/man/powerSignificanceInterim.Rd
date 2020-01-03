\name{powerSignificanceInterim}
\alias{powerSignificanceInterim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power of an Adaptative Replication Study}
%%  ~~function to do ... ~~
\description{
Computes the probability that a replication study yields a significant effect estimate in the specified direction at the final analysis given the data so far (accumulated up to the interim analysis).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
powerSignificanceInterim(to, ti, c=1, f = 1/2,
level=0.05, designPrior = "conditional", analysisPrior = "flat",  alternative="two.sided",
shrinkage=1)
}
\arguments{

  \item{to}{
    A vector of test statistics from original studies. 
  }
   \item{ti}{
    A vector of test statistics from interim analyses of replication studies. 
  }
  \item{c}{
    The ratio of the sample size of the
    replication study up to the interim analysis to the sample size of the original study. 
  }
  \item{f}{
  Fraction of the replicaton study already completed.
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
    Defaults to \code{"two.sided"}.  
    Specifies direction of the alternative. 
    If the alternative is one-sided, then sample size calculations are based on a one-sided assessment of significance.
  }
  \item{shrinkage}{
  A number in [0,1]. 
  Defaults to \code{1}. 
  Specifies how much the original effect estimate is shrunken towards zero (e.g. the effect is halfed for \code{shrinkage = 0.5}). 
  }
}
  \details{This is an extension of \code{PowerSignificance()} and adapts the `interim power' from Spiegelhalter et al. (2004) to the setting of replication studies.
    
    %%  ~~ If necessary, more details than the description above ~~
  }

  \references{
  Spiegelhalter, D. J., Abrams, K. R., and Myles, J. P. (2004). Bayesian Approaches to Clinical Trials    and Health-Care
  Evaluation, volume 13. John Wiley & Sons
   
  }
  \author{
    Charlotte Micheloud,
    Leonhard Held,
    Samuel Pawel 
  }
\seealso{
  \code{\link{sampleSizeSignificance}}, \code{\link{powerSignificance}}
}
  \examples{
  powerSignificanceInterim(to = 2, ti = 2, c = 1, f = 1/2, 
  designPrior = "conditional", analysisPrior = "flat")
  
  powerSignificanceInterim(to = 2, ti = 2, c = 1, f = 1/2, 
  designPrior = "predictive", analysisPrior = "flat")
  
  powerSignificanceInterim(to = 2, ti = 2, c = 1, f = 1/2, 
  designPrior = "flat", analysisPrior = "flat")
  
  powerSignificanceInterim(to = 2, ti = -2, c = 1, f = 1/2, 
  designPrior = "conditional", analysisPrior = "flat")
  
  powerSignificanceInterim(to = 2, ti = 2, c = 1, f = 1/2, 
  designPrior = "conditional", analysisPrior = "flat", shrinkage = 0.75)
  
    
}