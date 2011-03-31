\name{calcCaS-methods}
\docType{methods}
\alias{calcCaS}
\alias{calcCaS,xsAnnotate-method}
\title{EIC correlation grouping of LC/ESI-MS data}
\description{ Calculate the correlation across samples. Filtering correlation with specific parameters and returns a correlation matrix. 
}
\usage{
  calcCaS(object,corval=0.75, pval=0.05, intval="into")
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{corval}{Correlation threshold for positive hits}
  \item{pval}{P-Value threshold for significance level of correlation}
  \item{intval}{Selection of the intensity values that should be used in the correlation analysis. Can be into, maxo or intb.}
}
\details{
  Calculate pearson correlation between the peak intensites over all samples. Afterwards use cor.test for returning only significant correlation. Returns only those correlation,
  which are above both threshold. Set corval and pval to 0 to get the unfiltered correlation matrix.
  If the object is pregrouped with groupFWHM, then the correlation is only calculated between peaks within a pseudospectrum. Otherwise between all peaks.
}
\examples{
 library(CAMERA)
 #Multiple sample 
 library(faahKO)
 xs.grp       <- group(faahko)
 #create xsAnnotate object 
 xsa          <- xsAnnotate(xs.grp)
 #generate pseudospectra
 xsa.group    <- groupFWHM(xsa)
 #calculate correlation
 correlationMatrix <- calcCaS(xsa.group)
}

\value{
  A matrix with 4 columns:
  \item{x}{
   peak index according to peaktable
  }
  \item{y}{
   peak index according to peaktable
  }
  \item{cor}{
    correlation value between peak x and peak y
  }
  \item{ps}{
    pseudospektrum index for both peaks
  }
}

\seealso{
  \code{\link{calcCiS}}
  \code{\link{groupCorr}}
  \code{\link{xsAnnotate-class}}
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}