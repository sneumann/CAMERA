\name{findIsotopesWithValidation}
\docType{methods}
\alias{findIsotopesWithValidation}
\alias{findIsotopesWithValidation,xsAnnotate-method}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Annotate validated isotope clusters for a xsAnnotate object. Returns a xsAnnotate object with annotated isotopes. Validation of isotope clusters is based on statistics of the KEGG database implemented in S4 class object \code{compoundQuantiles}.
}
\usage{
  findIsotopesWithValidation(object, maxcharge=3, ppm=5, mzabs=0.01, intval=c("maxo","into","intb"), validateIsotopePatterns = TRUE)
}

\arguments{
  \item{object}{the \code{xsAnnotate} object}
  \item{maxcharge}{max. number of the isotope charge}
  \item{ppm}{ppm error for the search}
  \item{mzabs}{allowed variance for the search}
  \item{intval}{choose intensity values for C12/C13 check. Allowed values are into, maxo, intb }
  \item{validateIsotopePatterns}{logical, if TRUE putative isotope clusters are validated based on KEGG database statistics.}
}
\details{
  Isotope peaks are annotated for a xsAnnotate object according to given rules (maxcharge, maxiso). The algorithm benefits from a earlier grouping of the data, with groupFWHM. 
  Generates a list of all possible isotopes, which is stored in object@isotopes. Those isotope information will be used in the groupCorr funtion.
  The ratios between isotope peaks are checked against the mass--specific $99\%$ confidence interval based on statistics of the KEGG database. 
}
\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 an   <- findIsotopesWithValidation(an)
}
\author{Hendrik Treutler <hendrik.treutler@ipb-halle.de>}
\references{
  Hendrik Treutler and Steffen Neumann. 
  "Prediction, detection, and validation of isotope clusters in mass spectrometry data".
  Submitted to Metabolites 2016, Special Issue "Bioinformatics and Data Analysis".
}
\seealso{
  \code{\link{findIsotopes}}
}
\keyword{methods}
