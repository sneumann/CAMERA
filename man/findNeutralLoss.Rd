\name{findNeutralLoss}
\alias{findNeutralLoss}
\title{Find pseudospectra that contains a specific neutral loss}
\description{
 The method searches in every pseudospectra for a distance between two ions matching a provided mass difference. It returns a xcmsSet object containing the matching peaks.
}
\usage{
  findNeutralLoss(object, mzdiff=NULL, mzabs=0, mzppm=10)
}
\arguments{
  \item{object}{ xsAnnotate object}
  \item{mzdiff}{ neutral loss in Dalton}
  \item{mzabs}{ absolut allowed mass difference}
  \item{mzppm}{ relative allowed mass difference}
}
\details{
  The function needs a xsAnnotate object after groupCorr or groupFWHM. The resulting object is a artificial xcmsSet, where the peaks with the specific neutral loss are stored in xcmsSet@peaks.
}
\examples{
  library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  #Searches for Peaks with water loss
  xs.pseudo <- findNeutralLoss(an,mzdiff=18.01,mzabs=0.01) 
  xs.pseudo@peaks #show Hits
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}