\name{ReplicationProjects}
\alias{ReplicationProjects}
\docType{data}
\title{
  Data from four large-scale replication projects.
%%   ~~ data name/kind ... ~~
}
\description{
  Data from \emph{Reproduciblity Project Psychology} (RPP), \emph{Experimental Economics Replication Project} (EERP), \emph{Social Sciences Replication Project} (SSRP), \emph{Experimental Philosophy Replicability Project} (EPRP).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("ReplicationProjects")}
\format{
  A data frame with 143 observations on the following 13 variables:
  \describe{
    \item{\code{study}}{
    Study identifier, usually names of authors from original study
    }
    \item{\code{project}}{
    Name of replication project
    }
    \item{\code{r_O}}{
    Effect estimate of original study transformed to correlation scale
    }
    \item{\code{r_R}}{
    Effect estimate of replication study transformed to correlation scale
    }
    \item{\code{z_O}}{
    Effect estimate of original study transformed to Fisher-z scale
    }
    \item{\code{z_R}}{
    Effect estimate of replication study transformed to Fisher-z scale
    }
    \item{\code{z_se_O}}{
    Standard error of Fisher-z transformed effect estimate of original study
    }
    \item{\code{z_se_R}}{
    Standard error of Fisher-z transformed effect estimate of replication study
    }
    \item{\code{n_O}}{
    Sample size in original study
    }
    \item{\code{n_R}}{
    Sample size in replication study
    }
    \item{\code{pval_O}}{
    Two-sided P-value from significance test of effect estimate from original study
    }
    \item{\code{pval_R}}{
    Two-sided P-value from significance test of effect estimate from replication study
    }
    \item{\code{pm_belief}}{
    Peer belief about whether replication effect estimate will achieve statistical significance elicited through prediction market (only available for EERP and SSRP)
    }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  Two-sided P-values were calculated assuming normality of Fisher-z transformed effect estimates.
  From the RPP only the \emph{meta-analytic subset} is included, which consists of 73 out of 100 study pairs for which the standard error of the z-transformed correlation coeffient can be computed. 
  From the EPRP only 31 out of 40 study pairs are included where effective sample size for original and replication study are available simultaneously.
  For details about how the the data was preprocessed see supplement S1 of Pawel and Held (2019).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  Camerer, C. F., Dreber, A., Forsell, E., Ho, T.-H., Huber, J., Johannesson, M., ... Hang, W. (2016). Evaluating replicability of laboratory experiments in economics. \emph{Science}, 351, 1433–1436. doi:10.1126/science.aaf0918
  
  Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., ... Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science  between 2010 and 2015. \emph{Nature Human Behaviour}, 2, 637–644. doi:10.1038/s41562-018-0399-z
  
Cova,  F.,  Strickland,  B.,  Abatista,  A.,  Allard,  A.,  Andow,  J., Attie, M., ... Zhou, X. (2018). Estimating the reproducibility of experimental philosophy. \emph{Review of Philosophy and Psychology}. doi:10.1007/s13164-018-0400-9
 
  Open Science Collaboration. (2015). Estimating the reproducibility of psychological science. \emph{Science}, 349, aac4716. doi:10.1126/science.aac4716
  
  Pawel, S., Held, L. (2019). Probabilistic forecasting of replication studies. doi:10.31234/osf.io/fhwb7
}

  
\examples{
data("ReplicationProjects", package = "ReplicationSuccess")

## Plots of effect estimates
par(mfrow = c(2, 2), las = 1)
for (p in unique(ReplicationProjects$project)) {
  data_project <- subset(ReplicationProjects, project == p)
  plot(r_R ~ r_O, data = data_project, ylim = c(-0.5, 1), 
       xlim = c(-0.5, 1), main = p, xlab = expression(italic(r)[o]),
       ylab = expression(italic(r)[r]))
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
}

## Plots of peer beliefs
ReplicationProjects$significant <- factor(ReplicationProjects$pval_R < 0.05,
                                          levels = c(FALSE, TRUE),
                                          labels = c("no", "yes"))
par(mfrow = c(1, 2), las = 1)
for (p in c("Experimental Economics", "Social Sciences")) {
  data_project <- subset(ReplicationProjects, project == p)
  boxplot(pm_belief ~ significant, data = data_project, ylim = c(0, 1), 
          main = p, xlab = "Replication effect significant", ylab = "Peer belief")
  points(pm_belief ~ significant, data = data_project)
}


## Computing the sceptical p-value
ps <- with(ReplicationProjects,
           pSceptical(to = z_O/z_se_O, tr = z_R/z_se_R, c = z_se_O^2/z_se_R^2))
}
\keyword{datasets}