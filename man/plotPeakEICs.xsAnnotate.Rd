\name{plotPeakEICs-methods}
\docType{methods}
\alias{plotPeakEICs.xsAnnotate}
\alias{plotPeakEICs}
\alias{plotPeakEICs,xsAnnotate-method}
\title{Plot all extracted ion chromatograms of a compound spectrum}
\description{
  Plot all extracted ion chromatograms of one or more compound spectra 
  to the current graphics device. 
}
\usage{
  plotPeakEICs(object, pspec=1:length(object@pspectra),
    maxlabel=0, sleep=0,...)
}

\value{
  None.
}
\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 xsa   <- xsAnnotate(xs)
 xsa.grp   <- groupFWHM(xsa)
 #Plot all EICs of pseudospectrum 1
 plotPeakEICs(xsa.grp, pspec=1)
 }
 
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{pspec}{ pseudospectra indices that should be plotted }
  \item{maxlabel}{How many m/z labels to print}
  \item{sleep}{seconds to pause between plotting EICs of the different
  pseudospectra}
  \item{...}{other graphical parameters}
}
\author{Carsten Kuhl, \email{ckuhl@ipb-halle.de}}
\seealso{
  \code{\link{xsAnnotate-class}},
  \code{\link{png}},
  \code{\link{pdf}},
  \code{\link{postscript}},
}
\keyword{methods}
\keyword{hplot}
