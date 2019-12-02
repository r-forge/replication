\name{SSRP}
\alias{SSRP}
\docType{data}
\title{
  Data from the Social Sciences Replication Project (SSRP)
%%   ~~ data name/kind ... ~~
}
\description{
  Data from the \emph{Social Sciences Replication Project} (SSRP) including the details of the interim analysis. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("SSRP")}
\format{
  A data frame with 21 observations on the following 18 variables:
  \describe{
    \item{\code{study}}{
    Study identifier, usually names of authors from original study
    }
    \item{\code{r_O}}{
    Effect estimate of original study transformed to correlation scale
    }
    \item{\code{r_I}}{
    Effect estimate of replication study at the interim analysis transformed to correlation scale
    }
    \item{\code{r_R}}{
    Effect estimate of replication study at the final analysis transformed to correlation scale
    }
    \item{\code{z_O}}{
    Effect estimate of original study transformed to Fisher-z scale
    }
    \item{\code{z_I}}{
    Effect estimate of replication study at the interim analysis transformed to Fisher-z scale
    }
    \item{\code{z_R}}{
    Effect estimate of replication study at the final analysis transformed to Fisher-z scale
    }
    \item{\code{z_se_O}}{
    Standard error of Fisher-z transformed effect estimate of original study
    }
    \item{\code{z_se_I}}{
    Standard error of Fisher-z transformed effect estimate of replication study at the interim analysis
    }
    \item{\code{z_se_R}}{
    Standard error of Fisher-z transformed effect estimate of replication study at the final analysis
    }
    \item{\code{n_O}}{
    Sample size in original study
    }
    \item{\code{n_I}}{
    Sample size in replication study at the interim analysis
    }
    \item{\code{n_R}}{
    Sample size in replication study at the final analysis
    }
    \item{\code{pval_O}}{
    Two-sided P-value from significance test of effect estimate from original study
    }
    \item{\code{pval_I}}{
    Two-sided P-value from significance test of effect estimate from replication study at the interim analysis
    }
    \item{\code{pval_R}}{
    Two-sided P-value from significance test of effect estimate from replication study at the final analysis 
    }
    \item{\code{n75}}{
   Sample size calculated to have 90\% power in replication study to detect 75\% of the original effect size (expressed as the correlation coefficient r)
    }
    \item{\code{n50}}{
      Sample size calculated to have 90\% power in replication study to detect 50\% of the original effect size (expressed as the correlation coefficient r)
    }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  Two-sided P-values were calculated assuming normality of Fisher-z transformed effect estimates.
  A two-stage procedure was used for the replications. In stage 1, the authors had 90\% power to detect 75\% of the original effect size at the 5\% significance level in a two-sided test. If the original result replicated in stage 1 (two-sided P-value < 0.05 and effect in the same direction as in the original study), the data collection was stopped. If not, a second data collection was carried out in stage 2 to have 90\% power to detect 50\% of the original effect size for the first and the second data collections pooled. \code{n75} and \code{n50} are the planned sample sizes calculated to reach 90\% power in stage 1 and 2, respectively. They sometimes differ from the sample sizes that were actually collected (\code{n_I} and \code{n_R}, respectively). See supplementary information of Camerer et al. (2018) for details.
 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., ... Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science  between 2010 and 2015. \emph{Nature Human Behaviour}, 2, 637–644. doi:10.1038/s41562-018-0399-z
  
Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. doi:10.31234/osf.io/fhwb7
}

  
\examples{
data("SSRP", package = "ReplicationSuccess")

# plot of the sample sizes

plot(n_I~n_O, data = SSRP, ylim = c(0,2500), 
     xlim = c(0,400), xlab = expression(n[o]), 
     ylab = expression(n[i]))
abline(a = 0, b = 1, col = "grey")


plot(n_R~n_O, data = SSRP, ylim = c(0,2500), 
     xlim = c(0,400), xlab = expression(n[o]), 
     ylab = expression(n[r]))
abline(a = 0, b = 1, col = "grey")
}
\keyword{datasets}