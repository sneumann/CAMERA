\name{plotEICs-methods}
\docType{methods}
\alias{plotEICs.xsAnnotate}
\alias{plotEICs}
\alias{plotEICs,xsAnnotate-method}
\title{Plot extracted ion chromatograms from (multiple) Pseudospectra}
\description{
  Batch plot a list of extracted ion chromatograms to the current
  graphics device.
}
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{
    plotEICs(object,
    xraw,
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
  \item{xraw}{\code{xcmsRaw} object underlying the the xsAnnotate}
  \item{maxlabel}{How many m/z labels to print}
  \item{sleep}{seconds to pause between plotting EICs}
  \item{...}{other graphical parameters}
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
