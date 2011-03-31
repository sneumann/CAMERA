\name{calcCiS-methods}
\docType{methods}
\alias{calcCiS}
\alias{calcCiS,xsAnnotate-method}
\title{Calculate peak distance matrix after EIC correlation}
\description{
  Processing an xsAnnotate object and correlates peak EIC curves from one pseudospectrum, using a precalculated EIC matrix (\code{\link{getAllPeakEICs}}). It return a weighted edge list as distance matrix between peaks 
  according to the correlation analysis. The edge value is the pearson correlation coefficent. The list can be used as input for \code{\link{calcPC}}.
}
\usage{                                            
    calcCiS(object, EIC=EIC, corval=0.75, pval=0.05, psg_list=NULL)
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{EIC}{EIC Matrix}
  \item{corval}{Correlation threshold for the EIC correlation}
  \item{pval}{pvalue for testing correlation of significance}
  \item{psg_list}{Vector of pseudospectra indices. The correlation analysis will be only done for those groups}
}
\details{
  The algorithm correlates the EIC of a every peak with all others, to find the peaks that belong to one substance. LC/MS data should grouped with groupFWHM first. 
  This step reduce the runtime a lot and increased the number of correct classifications. Only correlation with a higher value than the correlation threshold and significant p-values
  will be returned.
}
\value{
  A matrix with 4 columns:

  \item{x}{
    peak index
  }
  \item{y}{
    peak index
  }
  \item{cor}{
    correlation value
  }
  \item{ps}{
    pseudospectrum index, which contains x and y
  }
}
\seealso{
  \code{\link{calcCaS}}
  \code{\link{groupCorr}}
  \code{\link{getAllPeakEICs}}
  \code{\link{xsAnnotate-class}}
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}