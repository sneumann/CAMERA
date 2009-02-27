\name{getPeaklist}
\alias{getPeaklist}
\title{Retrieve the annotatad peaklist}
\description{
 Extract all groups from a xsAnnotate object. Returns a peaklist with annotated peaks.
}
\usage{
  getPeaklist(object)
}
\arguments{
  \item{object}{ xsAnnotate object}
}
\details{
  xsAnnotate groups LC/MS Peaklist after there EIC correlation and FWHM. These function extract one of these so called "pseudo spectra groups" with include the peaklist with there annotations.
}
\examples{
library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(c(file), method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  an   <- findIsotopes(an)
  an   <- findAdducts(an,polarity="positive")
  peaklist <- getPeaklist(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}