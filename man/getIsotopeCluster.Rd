\name{getIsotopeCluster}
\alias{getIsotopeCluster}
\title{Retrieve the annotatad isotopes}
\description{
 Extract all annotated isotope cluster. Returns a list with one element per cluster.
 A element contains the charge of the molecule and a peakmatrix with mz and intensity value.
}
\usage{
  getIsotopeCluster(object, number=NULL, value="maxo", sampleIndex=NULL)
}
\arguments{
  \item{object}{ xsAnnotate object}
  \item{number}{ Set to NULL extract all isotope cluster or to specific chosen ones }
  \item{value}{ Which intensity values should be extracted. Allowed values are: maxo, into, intb}
  \item{sampleIndex}{ Selection vector with indexes to select from which sample(s)
  the intensity values should be retrieved. If set to NULL the sample is selected, which
  has been chosen for the pseudospectra in the grouping step }
}
\details{
  This method extract the isotope annotation from a xsAnnotate object. The order of the resulting list is the same as the one in the peaklist, see \code{\link{getPeaklist}}.
}

\examples{
  #single sample
  library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  an   <- findIsotopes(an) 
  isolist <- getIsotopeCluster(an)
  isolist[[10]] #get IsotopeCluster 10

  #multiple sample
  library(faahKO)
  xs <- group(faahko)
  xs <- fillPeaks(xs)
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  an   <- findIsotopes(an) 
  isolist <- getIsotopeCluster(an)

  #Select from multiple samples
  
  isolist <- getIsotopeCluster(an, sampleIndex=c(1,2,5))
  
  ##Interaction with Rdisop
 \dontrun{
  library(Rdisop)
  isotopes.decomposed <- lapply(isolist,function(x) {
    decomposeIsotopes(x$peaks[,1],x$peaks[,2],z=x$charge);
  }) #decomposed isotope cluster, filter steps are recommended
 }
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
