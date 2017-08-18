\name{getReducedPeaklist}
\docType{methods}
\alias{getReducedPeaklist}
\alias{getReducedPeaklist,xsAnnotate-method}
\title{Generate reduced peaklist from the annotatad peaklist}
\description{
 Extract information from an xsAnnotate object. Returns a reduced peaklist with annotated peaks. For any putative compound in the pcgroup, all found adducts are pooled into one putative compound per group. Thus, the reduced peaklist only contains one annotated adduct per pcgroup.
}
\usage{
  getReducedPeaklist(object, method = "median", intval = "into", default.adduct.info = "first", mzrt.range = FALSE, npeaks.sum = FALSE, cleanup = FALSE)
}
\arguments{
  \item{object}{ xsAnnotate object. }
  \item{method}{ Choose reduction method. Allowed values are "sum", "median", "maxint", "pca". }
  \item{intval}{ Choose intensity values. Allowed values are "into", "maxo", "intb". }
  \item{default.adduct.info}{ Choose method to select adduct information. Allowed values are "first", "maxint", "maxpeaks" }
  \item{mzrt.range}{ If TRUE, max and min values of mz and rt values of all adducts winthin a pcgroup are saved (not recommended). }
  \item{npeaks.sum}{ If TRUE, the sum of all peaks of all adducts within a pcgroup is saved (not recommended). }
  \item{cleanup}{ If TRUE, NA values and negative abundances are being set to zero and constant features (rows) are being removed. }
}
\details{
  This function extracts a reduced peaktable from an xsAnnotate object. Normally, all adducts are grouped for any putative compounds and saved within the peaklist (see method \link{getPeaklist}). However, for statistical computation it is sometimes better to only work with putative compounds rather than with all of their adducts. Thus, this function pools all adducts for any putative compound into one putative compound per pcgroup. There are several methods to choose from how this is being done.
  Selection methods:
    "sum": The intensities of adducts are summed for each sample.
    "median" (default): The median intensities of adducts is calculated for each sample.
    "maxint": Only the adduct with the highest intensities throughout the samples is returned.
    "pca": A Principal Component Analysis is being performed for the adducts for the samples. and the PC1 values are taken as intensity information.
  Select mz / rt methods:
    "first" (default): The mz & rt information of the first adduct are taken.
    "maxint": The mz & rt information of the adduct that has highest intensities are taken.
    "maxpeaks": The mz & rt information of the adduct that has the most peaks are taken.
    In addition, when mzrt.range is TRUE, the min and max values of all mz and rt found in a group are stored within mzmin, mzmax and rtmin and rtmax (not recommended).
    In addition, when npeaks.sum is TRUE, all peaks within a pcgroup are summed (not recommended).
}
\examples{
library(CAMERA)
file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findIsotopes(an)
an   <- findAdducts(an,polarity="positive")
peaklist.reduced <- getReducedPeaklist(an)
}
\author{Kristian Peters <kpeters@ipb-halle.de>}
\keyword{methods}

