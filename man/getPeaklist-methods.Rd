\name{getPeaklist}
\docType{methods}
\alias{getPeaklist}
\alias{getPeaklist,xsAnnotate-method}
\title{Generate the annotatad peaklist}
\description{
 Extract all information from an xsAnnotate object. Returns a peaklist with annotated peaks.
}
\usage{
  getPeaklist(object, intval="into")
}
\arguments{
  \item{object}{ xsAnnotate object}
  \item{intval}{ Choose intensity values. Allowed values are into, maxo, intb, intf, maxf}
}
\details{
  This function extract the peaktable from an xsAnnotate object, containing three additional columns (isotopes, adducts, pseudospectrum) with represents the annotation results. For a grouped xcmsSet it returns the grouped peaktable.
}
\examples{
library(CAMERA)
file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findIsotopes(an)
an   <- findAdducts(an,polarity="positive")
peaklist <- getPeaklist(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
