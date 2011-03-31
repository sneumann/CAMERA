\name{groupFWHM-methods}
\docType{methods}
\alias{groupFWHM}
\alias{groupFWHM,xsAnnotate-method}
\title{FWHM-Grouping of LC/ESI-MS data}
\description{
  Group peaks of a xsAnnotate object according to there retention time
  into pseudospectra-groups. Uses the peak FWHMs as grouping borders.
  Returns xsAnnotate object with pseudospectra informations.
}
\usage{
  groupFWHM(object, sigma = 6 , perfwhm = 0.6, intval = "maxo")
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{sigma}{the multiplier of the standard deviation}
  \item{perfwhm}{percentage of the width of the FWHM}
  \item{intval}{intensity values for ordering. Allowed values are into, maxo,
intb }
}
\details{
  Every peak who eluate at the same time as a selected peak, will be part of the
  group. Same time is defined about the Rt_med +/- FWHM * perfwhm.
  For a single sample \code{xcmsSet} the selection of peaks starts at the
  most abundant and goes down to the smaller ones. With a multiple sample set
  the automatic selection uses the most abundant peak as an representative for
  every feature group, according to the xcms grouping. With the xsAnnotate
  sample parameter a sample selection can be defined to use only specific
  samples. See \code{\link{xsAnnotate-class}} for further information.
  The FWHM (full width at half maximum) of a peak is estimated as FWHM = SD *
  2.35. For the calculation of the SD, the peak is assumed as normal
  distributed.
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

 #With specific selected sample
 xs.anno  <- xsAnnotate(xs, sample=1)
 xs.group <- groupFWHM(xs.anno)
 
 #With automatic selection
 xs.anno.auto  <- xsAnnotate(xs)
 xs.group.auto <- groupFWHM(xs.anno.auto)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}