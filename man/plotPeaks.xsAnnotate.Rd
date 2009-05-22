\name{plotPeaks-methods}
\docType{methods}
\alias{plotPeaks.xsAnnotate}
\alias{plotPeaks}
\alias{plotPeaks,xsAnnotate-method}
\title{Plot a Pseudospectrum}
\description{
  Plot a pseudospectrum, with the most intense peaks labelled,
  to the current graphics device.
}
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{
    plotPeaks(object, pspec=NULL, log=FALSE,
                       value="into", maxlabel=0)
}
}}}
\value{
  None.
}
\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{pspec}{ID of the pseudospectrum to print}
  \item{log}{Boolean, whether the log(intensity) should be shown}
  \item{value}{Which of a peak's intensities should be used}
  \item{maxlabel}{How many m/z labels to print}
}
\author{Steffen Neumann, \email{sneumann@ipb-halle.de}}
\seealso{
  \code{\link{xsAnnotate-class}},
  \code{\link{png}},
  \code{\link{pdf}},
  \code{\link{postscript}},
}
\keyword{methods}
\keyword{hplot}
