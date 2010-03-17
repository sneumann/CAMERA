\name{groupCorr-methods}
\docType{methods}
\alias{groupCorr}
\alias{groupCorr,xsAnnotate-method}
\title{EIC correlation grouping of LC/ESI-MS data}
\description{
 Grouping the peaks after the correlation of the EICs into pseudospectragroups for a xsAnnotate object. Return a xsAnnotate object with group information.
}
\usage{
 groupCorr(object, cor_eic_th=0.75, psg_list = NULL, polarity = NA)
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{cor_eic_th}{Correlation threshold for the EIC correlation}
  \item{psg_list}{Vector of pseudospectra indices. The correlation analysis will be only done for those groups}
  \item{polarity}{Optional feature, that ensure that [M+H]+, [M+Na]+, [M+K]+ in pos. mode and [M-H]- [M+Cl]- ions  in neg. mode will not be seperated in different groups}
}
\details{
  The algorithm correlates the EIC of a every peak with all others, to find the peaks that belong to one substance. LC/MS data should grouped with groupFWHM first. This step reduce the runtime a lot and increased the number of correct classifications. If you perform the findIsotope step first, annotated isotopes will always be in the same group.
}
\examples{
 library(CAMERA)
 file  <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs    <- xcmsSet(c(file), method="centWave", ppm=30, peakwidth=c(5,10))
 an    <- xsAnnotate(xs)
 an    <- groupFWHM(an)
 an    <- findIsotopes(an) #optional step
 index <- c(1,4) #Only group one and four will be calculate
 an    <- groupCorr(an, psg_list=index, polarity="positive")
 #For csv output
 # write.csv(file="peaklist_with_isotopes.csv",getPeaklist(an))
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}