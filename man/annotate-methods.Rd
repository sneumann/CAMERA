\name{annotate-methods}
\docType{methods}
\alias{annotate}
\alias{annotate,xcmsSet-method}
\title{Automatic deconvolution/annotation of LC/ESI-MS data}
\description{
Wrapper skript for automatic annotation of isotope peaks, adducts and fragments for a (grouped) xcmsSet \code{xs}. The function returns an xsAnnotate object.
}
\usage{
annotate(object, sample=NA, nSlaves=1, sigma=6, perfwhm=0.6,
  cor_eic_th=0.75, graphMethod="hcs", pval=0.05, calcCiS=TRUE,
  calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5,
  ppm=5, mzabs=0.015, quick=FALSE, psg_list=NULL, rules=NULL,
  polarity="positive", multiplier=3, max_peaks=100 ,intval="into")
}
\arguments{
  \item{object}{xcmsSet with peak group assignments}
  \item{sample}{xsAnnotate: Sample selection for grouped xcmsSet, see \link{xsAnnotate-class}}
  \item{nSlaves}{xsAnnotate: Use parallel CAMERA mode, require Rmpi}
  \item{sigma}{groupFWHM: multiplier of the standard deviation}
  \item{perfwhm}{groupFWHM: percentage of FWHM width}
  \item{cor_eic_th}{groupCorr: correlation threshold (0..1)}
  \item{graphMethod}{groupCorr: Method selection for grouping peaks after correlation analysis into pseudospectra}
  \item{pval}{groupCorr: significant correlation threshold}
  \item{calcCiS}{groupCorr: Use correlation inside samples for peak grouping}
  \item{calcIso}{groupCorr: Use isotopic relationship for peak grouping}
  \item{calcCaS}{groupCorr: Use correlation across samples for peak grouping}
  \item{maxcharge}{findIsotopes: max. ion charge}
  \item{maxiso}{findIsotopes: max. number of expected isotopes}
  \item{minfrac}{findIsotopes: The percentage number of samples, which must satisfy the C12/C13 rule for isotope annotation}
  \item{ppm}{General ppm error}
  \item{mzabs}{General absolut error in m/z}  
  \item{quick}{Use only groupFWHM and findIsotopes}
  \item{psg_list}{Calculation will only be done for the selected groups}
  \item{rules}{findAdducts: User defined ruleset}
  \item{polarity}{findAdducts: Which polarity mode was used for measuring of the ms sample}
  \item{multiplier}{findAdducts: If no ruleset is provided, calculate ruleset with max. number n of [nM+x] clusterions}
  \item{max_peaks}{How much peaks will be calculated in every thread using the parallel mode}
  \item{intval}{General used intensity value (into, maxo, intb)}
}
\details{
  Batch script for annotation of an (grouped) xcmsSet \code{xs}. Generates an xsAnnotate object by calling all involved functions for the annotation step.
  Function list: 1: groupFWHM() , 2: findIsotopes() , 3: groupCorr(), 4: findAdducts()
  Return the xsAnnotate object, which inherits all annotations.
  For more information about the parameters see the specific function manpages.
}
\value{
  \code{annotate} returns an xsAnnotate object. For more information about the xsAnnotate object see \link{xsAnnotate-class}.
}

\examples{
 library(CAMERA)
 file <- system.file('mzML/MM14.mzML', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 xsa  <- annotate(xs)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
