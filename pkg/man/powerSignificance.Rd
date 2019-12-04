\name{powerSignificance}
\alias{powerSignificance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replication probabilities for P-values}
%%  ~~function to do ... ~~
\description{
  Computes the probability that a replication study yields a significant effect
  estimate in the specified direction. %in the same direction as in the original study. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
powerSignificance(po=NULL, to=p2t(po, alternative=alternative), c=1,
level=0.05, designPrior="conditional", alternative="two.sided", d=0,
shrinkage=1)
}
\arguments{
  \item{po}{
    A vector of P-values from original studies, default is \code{NULL}
  }
  \item{to}{
    A vector of test statistics from original studies. 
    Is calculated from \code{po}, if necessary. 
  }
  \item{c}{
    The ratio of the variances of the original and replication effect
    estimates. This is usually the ratio of the sample size of the
    replication study to the sample size of the original study. 
  }
  \item{level}{
    Significance level.
  }
  \item{designPrior}{
    Either \code{"conditional"}, \code{"predictive"}, or \code{"EB"}.  
    Defaults to \code{"conditional"}. If \code{"EB"}, the power is computed under a predictive distribution where the contribution of the original study is shrunken towards zero based on the evidence in the original study (with an empirical Bayes shrinkage estimator).
  }
  \item{alternative}{
    Either \code{"less"}, \code{"greater"}, or \code{"two.sided"}.  
    Defaults to \code{"two.sided"}.  
    Specifies direction of the alternative. 
    If the alternative is one-sided, then sample size calculations are based on a one-sided assessment of significance.
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  Is only taken into account when \code{designPrior = "predictive"} or \code{designPrior = "EB"}.
  } 
  \item{shrinkage}{
  A number in [0,1]. 
  Defaults to \code{1}. 
  Specifies how much the original effect estimate is shrunken towards zero (e.g. the effect is halfed for \code{shrinkage = 0.5}). 
  Recommended to set to 1 when designPrior is \code{"EB"}.
  }
}
  \details{
    This extends the 'replication probability', the probability of repeating a statistically significant result by Goodman (1992), to the case of possibly unequal sample sizes (see Senn (2002)), possible between study heterogeneity and shrinkage (see Pawel and Held (2019)). 
    
    %%  ~~ If necessary, more details than the description above ~~
  }
  \value{
    The probability that a replication study yields a significant effect estimate in the specified direction. 
  }
  \references{
    Goodman, S. N. (1992). A comment on replication, p-values and evidence, \emph{Statistics in Medicine}, \bold{11}, 875--879. \url{https://dx.doi.org/10.1002/sim.4780110705}
    
    Senn, S. (2002). Letter to the Editor, \emph{Statistics in Medicine}, \bold{21}, 2437--2444. 

    Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  To appear in \emph{Journal of the Royal Statistical Society, Series A}. 
  \url{https://arxiv.org/abs/1811.10287}
  
  Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. \url{https://dx.doi.org/10.31234/osf.io/fhwb7}

  }
  \author{
    Leonhard Held,
    Samuel Pawel (shrinkage and heterogeneity options)
  }
  \examples{
    powerSignificance(po=0.005, c=2)
    powerSignificance(po=0.005, c=2, designPrior="predictive")
    powerSignificance(po=0.005, c=2, alternative="greater")
    powerSignificance(to=-3, c=2, designPrior="predictive", alternative="less")
    powerSignificance(po=0.005, c=1/2)
    powerSignificance(po=0.005, c=1/2, designPrior="predictive")
    powerSignificance(po=0.005, c=1/2, alternative="greater")
    powerSignificance(po=0.005, c=1/2, designPrior="predictive", alternative="greater")
    powerSignificance(po=0.005, c=1/2, designPrior="predictive", alternative="greater", 
                      d = 0.5, shrinkage=0.5)
    powerSignificance(po=0.005, c=1/2, designPrior="EB", alternative="greater", 
                      d = 0.5)
}