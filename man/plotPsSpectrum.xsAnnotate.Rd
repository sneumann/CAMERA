\name{plotPsSpectrum-methods}
\docType{methods}
\alias{plotPsSpectrum.xsAnnotate}
\alias{plotPsSpectrum}
\alias{plotPsSpectrum,xsAnnotate-method}
\title{Plot a Pseudospectrum}
\description{
  Plot a pseudospectrum, with the most intense peaks labelled,
  to the current graphics device.
}
\usage{
  plotPsSpectrum(object, pspec=1:length(object@pspectra), log=FALSE, value="into", maxlabel=0, title=NULL,mzrange=numeric(), sleep=0, cexMulti = 1,...)
}

\section{Methods}{
\describe{
\item{\code{signature(object = "xsAnnotate")}}{
  object deriviving from class "xsAnnotate"
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
  \item{title}{Main title of the Plot}
  \item{mzrange}{Which m/z range should plotted}
  \item{sleep}{Time (in seconds) to wait between successive Spectra, if
    multiple \code{pspec} are requested.}
  \item{cexMulti}{Cex multiplier for peak labels}
  \item{...}{Additional parameter for function plot}
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
