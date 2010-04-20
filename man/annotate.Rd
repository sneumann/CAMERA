\name{annotate}
\alias{annotate}
\alias{annotate,xsAnnotate-methods}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Annotate isotope peaks, adducts and fragments for a (grouped) xcmsSet \code{xs}. Returns a xsAnnotate object.
}
\usage{
annotate(object, sigma = 6, perfwhm = 0.6, cor_eic_th = 0.75, maxcharge = 3, maxiso = 4, ppm = 5, mzabs = 0.01, multiplier = 3, sample=NA, quick= FALSE, psg_list=NULL, polarity="positive", nSlaves = 1, max_peaks=100)
}
\arguments{
  \item{object}{ xcmsSet with peak group assignments}
  \item{sigma}{Isotopic peak relationship table}
  \item{perfwhm}{Adduct/Fragment peak relationship table}
  \item{cor_eic_th}{correlation threshold (0..1)}
  \item{maxcharge}{max charge of the ions}
  \item{maxiso}{max number of isotopes}
  \item{ppm}{ppm error}
  \item{mzabs}{absolut error}
  \item{multiplier}{max. number n of [nM+x] clusterions}
  \item{sample}{Index of which sample is used for the correlation, NA for auto-selection}
  \item{quick}{Use only groupFWHM and findIsotopes}
  \item{psg_list}{Calculation will only be done for the selected groups}
  \item{polarity}{Which polarity mode was used for measuring of the ms sample}
  \item{nSlaves}{Number of slaves for parallel calculation (Warning: Beta)}
  \item{max_peaks}{If run in parallel mode, value defines how much peaks will be calculated in every thread }
}
\details{
  Batch script for a annotation for a (grouped) xcmsSet \code{xs}. Generate intern a xsAnnotate object and calls the member function for the annotation step.
  Function list: 1: groupFWHM() , 2: findIsotopes() , 3: groupCorr(), 4: findAdducts()
  Return the xsAnnotate object which all the annotations.
  For more information see the funtion manpages.
}
\value{
\code{annotate}	returns an xsAnnotate object. For more information about see \link{xsAnnotate-class}
}

\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 xsa  <- annotate(xs)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
