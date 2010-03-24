\name{findIsotopes-methods}
\docType{methods}
\alias{findIsotopes}
\alias{findIsotopes,xsAnnotate-method}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Annotate isotope peaks for a xsAnnotate object. Returns a xsAnnotate object with annotated isotopes.
}
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{findIsotopes(object, maxcharge=3, maxiso=3, ppm=5, mzabs=0.01)}
}
}}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{maxcharge}{max. number of the isotope charge}
  \item{maxiso}{max. number of the isotope peaks}
  \item{ppm}{ppm error for the search}
  \item{mzabs}{allowed variance for the search}
}
\details{
  Isotope peaks are annotated for a xsAnnotate object according to given rules (maxcharge, maxiso). The algorithm benefits from a earlier grouping of the data, with groupFWHM. Generates a list of all possible isotopes, which is stored in object@isotopes. Those isotope information will be used in the groupCorr funtion. 
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