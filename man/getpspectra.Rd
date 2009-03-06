\name{getpspectra}
\alias{getpspectra}
\title{Retrieve a spectra-group peaklist}
\description{
 Extract one group from a xsAnnotate object. Returns a peaklist with annotated peaks.
}
\usage{
  getpspectra(object,grp)
}
\arguments{
  \item{object}{xsAnnotate object}
  \item{grp}{index of pseudo-spectra-group}
}
\details{
  xsAnnotate groups LC/MS Peaklist after there EIC correlation and FWHM. These function extract one of these so called "pseudo spectra groups" with include the peaklist with there annotations. The annotation depends on a before called findAdducts() ( and findIsotopes() ).
Attention: The indices for the isotopes, are those from the whole peaklist. See \code{getPeaklist()}.
}
\examples{
  library(CAMERA)
  file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
  xs   <- xcmsSet(c(file), method="centWave", ppm=30, peakwidth=c(5,10))
  an   <- xsAnnotate(xs)
  an   <- groupFWHM(an)
  peaklist <- getpspectra(an, 1)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}