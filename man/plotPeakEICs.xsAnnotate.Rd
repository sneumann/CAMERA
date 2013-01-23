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
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{
    plotEICs(object,
    pspec=1:length(object@pspectra),
    maxlabel=0, sleep=0)}
  }
}
}
\value{
  None.
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{maxlabel}{How many m/z labels to print}
  \item{sleep}{seconds to pause between plotting EICs}
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
