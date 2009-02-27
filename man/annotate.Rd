\name{annotate}
\alias{annotate}
\alias{annotate,xsAnnotate-methods}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Annotate isotope peaks, adducts and fragments for a (grouped) xcmsSet \code{xs}. Returns a xsAnnotate object.
}
\usage{
annotate(xs, sigma = 6, perfwhm = 0.6, cor_eic_th = 0.75, maxcharge = 3, maxiso = 4, ppm = 5, mzabs = 0.01, multiplier = 3, sample=1, category = NA, polarity="positive")
}
\arguments{
  \item{xs}{ xcmsSet with peak group assignments}
  \item{sigma}{Isotopic peak relationship table}
  \item{perfwhm}{Adduct/Fragment peak relationship table}
  \item{cor_eic_th}{correlation threshold (0..1)}
  \item{maxcharge}{max charge of the ions}
  \item{maxiso}{max number of isotopes}
  \item{ppm}{ppm error}
  \item{mzabs}{absolut error}
  \item{multiplier}{max. number n of [nM+x] clusterions}
  \item{sample}{Index of which sample is used for the correlation}
  \item{category}{Which class label should be used}
  \item{polarity}{Which polarity mode was used for measuring of the ms sample}
}
\details{
  Batch script for a annotation for a a (grouped) xcmsSet \code{xs}. Generate intern a xsAnnotate object and calls the member function for the annotation step. Returns a peaklist with additonal columns for isotopes, adducts and an index for the pseudo spectra group. xsAnnotate groups Peaks into spectra groups, after there Retentiontime and EIC correlations.
  Attention: For the EIC correlation only one sample can be used, so if the xcmsSet contains more than one sample, one must be chosen.
}
\value{
\code{annotate}	returns an xsAnnotate object. For more information about see \link{xsAnnotate-class}
}

\examples{
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 annotated_xs<- annotate(xs)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
