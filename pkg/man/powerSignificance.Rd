\name{powerSignificance}
\alias{powerSignificance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replication probabilities for P-values}
%%  ~~function to do ... ~~
\description{
  Computes the probability that a replication study yields a significant effect
  estimate in the same direction as in the original study. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
powerSignificance(po=NULL, to=p2t(po, alternative=alternative), c=1,
level=0.05, designPrior="conditional", alternative="two.sided", 
shrinkage=FALSE, d=0)
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
    Either \code{"conditional"} or \code{"predictive"}.  
    Defaults to \code{"conditional"}.
  }
  \item{alternative}{
    Either \code{"one.sided"} or \code{"two.sided"}.  
    Defaults to \code{"two.sided"}.  
    Specifies if the P-value is one-sided or two-sided. 
    If the P-value is one-sided, then sample size calculations are based on a one-sided assessment of significance.
  }
  \item{shrinkage}{
  Either \code{TRUE} or \code{FALSE}. 
  Defaults to \code{FALSE}. If \code{TRUE} and \code{designPrior="predictive"}, the probability of significance is computed under a predictive distribution where the original effect is shrunken towards zero based on the evidence in the original study (with an empirical Bayes shrinkage estimator).
  }
  \item{d}{
  The relative between-study heterogeneity, i.e. the ratio of the heterogeneity variance to the variance of the original effect estimate.
  }  
}
  \details{
    This extends the 'replication probability', the probability of repeating a statistically significant result by Goodman (1992), to the case of possibly unequal sample sizes and possible between study heterogeneity, see also Senn (2002). 
    
    %%  ~~ If necessary, more details than the description above ~~
  }
  \value{
    The probability that a replication study yields a significant effect estimate in the same direction as in the original study. 
  }
  \references{
    Goodman, S. N. (1992). A comment on replication, p-values and evidence, \emph{Statistics in Medicine}, \bold{11}, 875--879. 
    
    Senn, S. (2002). Letter to the Editor, \emph{Statistics in Medicine}, \bold{21}, 2437--2444. 

    Held, L. (2019). A new standard for the analysis and design of replication studies (with discussion). 
  To appear in \emph{Journal of the Royal Statistical Society, Series A}. 
  \url{https://arxiv.org/abs/1811.10287}

  }
  \author{
    Leonhard Held,
    Samuel Pawel (shrinkage and heterogeneity options)
  }
  \examples{
    powerSignificance(po=0.005, c=2)
    powerSignificance(po=0.005, c=2, designPrior="predictive")
    powerSignificance(po=0.005, c=2, alternative="one.sided")
    powerSignificance(po=0.005, c=2, designPrior="predictive", alternative="one.sided")
    
    powerSignificance(po=0.005, c=1/2)
    powerSignificance(po=0.005, c=1/2, designPrior="predictive")
    powerSignificance(po=0.005, c=1/2, alternative="one.sided")
    powerSignificance(po=0.005, c=1/2, designPrior="predictive", alternative="one.sided")
    powerSignificance(po=0.005, c=1/2, designPrior="predictive", alternative="one.sided", shrinkage=TRUE, d = 1)

}