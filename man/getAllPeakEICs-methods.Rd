\name{getAllPeakEICs}
\docType{methods}
\alias{getAllPeakEICs}
\alias{getAllPeakEICs,xsAnnotate-method}
\title{Generate EIC information from raw data}
\description{
 Generate EIC data out of the raw data, according to the peak peaker information.
}
\usage{
  getAllPeakEICs(object, index, maxscans)
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{index}{Sample index vector, with the same length as the number of peaks. 
  Encoding from with sample the peak should be extracted. If all peaks should be generated from the
  same sample set index = rep(sample index, peak count)}
  \item{maxscans}{Maximum number of scantime of all samples. If not available, it can be determined with CAMERA:::getMaxScans(object).}
}
\details{
  The function extract from the raw data the EIC curves. Therefore all .netcdf, .mzdata etc. files must be acessable.
  It returns a list with two item.
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
 #get maxScans
 maxscans <- CAMERA:::getMaxScans(xsa.group)
 #calculate correlation
 tmp <- getAllPeakEICs(xsa.group,index=rep(1,nrow(xsa.group@groupInfo)),maxscans=maxscans)
 #extract EIC matrix
 EIC.matrix <- tmp$EIC;
}

\value{
  A list with items:
  \item{EIC}{
   EIC Matrix with rows = number of peaks and columns = maxscans. It contains mostly NA values and only in that part, where
  a peak had been found, the intensity information.
  }
  \item{scantimes}{
    Not used anymore, ignore!
  }
}

\seealso{
  \code{\link{xsAnnotate-class}}
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}