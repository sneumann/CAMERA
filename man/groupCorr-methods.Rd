\name{groupCorr-methods}
\docType{methods}
\alias{groupCorr}
\alias{groupCorr,xsAnnotate-method}
\title{EIC correlation grouping of LC/ESI-MS data}
\description{
 Grouping the peaks after the correlation of the EICs into pseudospectragroups for a xsAnnotate object. Return a xsAnnotate object with group information.
}
\usage{
 groupCorr(object,cor_eic_th=0.75)
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{cor_eic_th}{correlation threshold for the EIC correlation}
}
\details{
  The algorithm correlates the EIC of a every peak with all others, to find the peaks that belong to one substance. LC/MS data should grouped with groupFWHM first. This step reduce the runtime a lot and increased the number of correct classifications.
}
\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(c(file), method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 an   <- groupCorr(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}