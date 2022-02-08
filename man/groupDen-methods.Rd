\name{groupDen}
\docType{methods}
\alias{groupDen}
\alias{groupDen,xsAnnotate-method}
\title{Density-Grouping of LC/ESI-MS data}
\description{
  Group peaks of a xsAnnotate object according to peak
  distributions in chromatographic time 
  into pseudospectra-groups. Works analogous as the group.density
  method of xcms. Returns xsAnnotate object with pseudospectra informations.
}
\usage{
  groupDen(object, bw = 5 , ...)
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{bw}{bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram}
  \item{...}{Further Arguments, NYI}
}
\value{
  Returns a grouped xsAnnotate object. 
}
\examples{
 library(CAMERA)
 #Single sample 
 file <- system.file('mzML/MM14.mzML', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 xsa   <- xsAnnotate(xs)
 xsa.grp   <- groupDen(xsa, bw=0.5)

 #Multiple sample 
 library(faahKO)
 xs   <- group(faahko)

 #With specific selected sample
 xsa     <- xsAnnotate(xs, sample=1)
 xsa.grp <- groupDen(xsa)
 
 #With automatic selection
 xsa.auto     <- xsAnnotate(xs)
 xsa.grp.auto <- groupDen(xsa.auto)
}

\details{
  The grouping strongly depends on the bw parameter. 
  For an UPLC a good starting point is smaller or around 1.
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}