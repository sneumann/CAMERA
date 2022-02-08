\name{getpspectra}
\alias{getpspectra}
\title{Retrieve a peaklist of one or more pseudospectra}
\description{
 Extract group(s) from a xsAnnotate object. Returns a peaklist as matrix with annotated peaks.
}
\usage{
  getpspectra(object, grp)
}
\arguments{
  \item{object}{xsAnnotate object}
  \item{grp}{index of pseudo-spectra-group}
}
\details{
  xsAnnotate groups LC/MS Peaklist after there EIC correlation and FWHM. These function extract one or more of these so called "pseudo spectra groups" with include the peaklist with there annotations. The annotation depends on a before called findAdducts() ( and findIsotopes() ).
  Important: The indices for the isotopes, are those from the whole peaklist. See \code{getPeaklist()}.
}
\examples{
  library(CAMERA)
  file <- system.file('mzML/MM14.mzML', package = "CAMERA")
  xs   <- xcmsSet(c(file), method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  #For one group
  peaklist <- getpspectra(an, 1)
  #For two groups
  peaklist <- getpspectra(an, c(1,2))
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}