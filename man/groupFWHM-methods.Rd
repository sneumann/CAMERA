\name{groupFWHM-methods}
\docType{methods}
\alias{groupFWHM}
\alias{groupFWHM,xsAnnotate-method}
\title{FWHM-Grouping of LC/ESI-MS data}
\description{
  Group peaks of one xsAnnotate object after the FWHM of there retentiontimes into pseudospectra-groups. Returns a xsAnnotate object with group information.
}
\usage{
  groupFWHM(object, sigma = 6 , perfwhm = 0.6)
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{sigma}{the multiplier of the standard deviation}
  \item{perfwhm}{percentage of the width of the FWHM}
}
\details{
  Every peak who eluated at the same time as a selected peak, will be part of the group. Same time is defined about the Rt_med +/- FWHM * perfwhm.
  For single sample the selection of peaks starts at the most abundant and goes down the smaller ones. In a multiple sample set you can define a sample, 
  that will be used for the process or use the automatic selection. This select the most abundant peak as an representative for every peak group. 
  Peak group means here those ones from the xcms grouping.
  The FWHM (full width at half maximum) of a peak is estimated as FWHM = SD * 2.35. The calculation of the SD, the peak is assumpted as a normal distribution.
}
\examples{
 library(CAMERA)
 #Single sample 
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)

 #Multiple sample 
 library(faahKO)
 xs   <- group(faahko)
 #With selected sample
 xs.anno  <- xsAnnotate(xs, sample=1)
 xs.group <- groupFWHM(xs.anno)
 
 #With automatic selection
 xs.anno.auto  <- xsAnnotate(xs)
 xs.group.auto <- groupFWHM(xs.anno.auto)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}