\name{findIsotopes}
\docType{methods}
\alias{findIsotopes}
\alias{findIsotopes,xsAnnotate-method}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Annotate isotope peaks for a xsAnnotate object. Returns a xsAnnotate object with annotated isotopes.
}
\usage{
  findIsotopes(object, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, intval="maxo", minfrac=0.5)
}

\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{maxcharge}{max. number of the isotope charge}
  \item{maxiso}{max. number of the isotope peaks}
  \item{ppm}{ppm error for the search}
  \item{mzabs}{allowed variance for the search}
  \item{intval}{choose intensity values for C12/C13 check. Allowed values are into, maxo, intb }
  \item{minfrac}{in case of multiple samples, percentaged value of samples, which have to contain the correct C12/C13 ratio}
}
\details{
  Isotope peaks are annotated for a xsAnnotate object according to given rules (maxcharge, maxiso). The algorithm benefits from a earlier grouping of the data, with groupFWHM. 
  Generates a list of all possible isotopes, which is stored in object@isotopes. Those isotope information will be used in the groupCorr funtion.
  The itensity of the C13 isotope peak is checked against the C12 of proper ratio. In the case of mulitiple sample, all samples will be tested. Minfrac describe the minimal percentaged
  of samples, which must passed the test. If peaks are NA, then this sample is skipped and the ratio is (found correct C12/C13 ratio) / (samples containing C12 and C13 peak).
  
}
\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 an   <- findIsotopes(an)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}