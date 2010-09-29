\name{findNeutralLossSpecs}
\alias{findNeutralLossSpecs}
\title{Find pseudospectra that contains a specific neutral loss}
\description{
 The method searches in every pseudospectra for a distance between two ions matching a provided mass difference. It returns a boolean vector with the length equals to the number of pseudospectra,
 where a hit is marked with TRUE.
}
\usage{
  findNeutralLossSpecs(object, mzdiff=NULL, mzabs=0, mzppm=10)
}
\arguments{
  \item{object}{ xsAnnotate object}
  \item{mzdiff}{ neutral loss in Dalton}
  \item{mzabs}{ absolut allowed mass difference}
  \item{mzppm}{ relative allowed mass difference}
}
\details{
  The function needs a xsAnnotate object after groupCorr or groupFWHM.
}
\examples{
  library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  hits <- findNeutralLossSpecs(an,mzdiff=18.01,mzabs=0.01) #Searches for Pseudspecta with water loss
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}