\name{groupFWHM-methods}
\docType{methods}
\alias{groupFWHM}
\alias{groupFWHM,xsAnnotate-method}
\title{FWHM-Grouping of LC/ESI-MS data}
\description{
  Grouping the peaks after the FWHM of the retentiontimes into pseudospectragroups for a xsAnnotate object. Returns a xsAnnotate object with group information.
}
\usage{
 groupFWHM(object,sigma = 6 , perfwhm = 0.6)
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{sigma}{the multiplier of the standard deviation}
  \item{perfwhm}{percentage of the width of the FWHM}
}
\details{ 
  The FWHM (full width at half maximum) of a peak, will be used for the grouping. Every peak who eluated at the same time as a selected peak, will be part of the group. Same time is defined about the Rt-med +/-FWHM*perfwhm. FWHM is calculated as the product of sigma * standard deviation.
}
\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}