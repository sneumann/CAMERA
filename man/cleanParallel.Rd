\name{cleanParallel}
\alias{cleanParallel}
\title{Cleans up with spawned slave processes after use}
\description{
 The spawned slaves processes, which are created within the parallel mode, are closed explicit. 
}
\usage{
  cleanParallel(object)
}
\arguments{
  \item{object}{ xsAnnotate object}
}
\details{
  The function needs a xsAnnotate object after groupCorr or groupFWHM. The resulting object is a artificial xcmsSet, where the peaks with the specific neutral loss are stored in xcmsSet@peaks.
}
\examples{
  library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs, polarity="positive", nSlaves=2)
  an   <- groupFWHM(an)
  an   <- findAdducts(an)
  cleanParallel(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}