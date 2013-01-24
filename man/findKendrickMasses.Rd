\name{findKendrickMasses}
\alias{findKendrickMasses}
\title{Find specfic mass defects using Kendrick mass scales}
\description{
 Todo
}
\usage{
  findKendrickMasses(object, masses=c(14, 14.01565),
  maxHomologue=4, error=0.002, time=60, intval="maxo",
  plot=FALSE)
}
\arguments{
  \item{object}{ xsAnnotate object}
  \item{masses}{ nominal mass and exact mass}
  \item{error}{ allowed mass difference in Da for matching Kendrick mass defect}
  \item{maxHomologue}{ max number of homologue}
  \item{time}{ allowed retention time difference between homologues}
  \item{intval}{ intensity value (allowed values: maxo,into or intb)}
  \item{plot}{ plot hits}
}
\details{
}
\examples{
  library(CAMERA)
  library(faahKO)
  xs   <- group(faahko)

  #With specific selected sample
  xsa     <- xsAnnotate(xs)
  #Screen for substance with CH2 differences
  findKendrickMasses(xsa, masses=c(14, 14.01565), plot=TRUE)
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}