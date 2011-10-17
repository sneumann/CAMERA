\name{annotateDiffreport}
\alias{annotateDiffreport}
\alias{annotateDiffreport,xsAnnotate-methods}
\title{Automatic deconvolution/annotation of LC/ESI-MS data}
\description{
 Wrapper function for the xcms diffreport and the annotate function. Returns a diffreport within the annotation results.
}
\usage{
annotateDiffreport(object, sample=NA, nSlaves=1, sigma=6, perfwhm=0.6,
  cor_eic_th=0.75, graphMethod="hcs", pval=0.05, calcCiS=TRUE,
  calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5,
  ppm=5, mzabs=0.015, quick=FALSE, psg_list=NULL, rules=NULL,
  polarity="positive", multiplier=3, max_peaks=100, intval="into",
  pval_th = NULL, fc_th = NULL, sortpval=TRUE, ...)
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
  \item{pval_th}{pval threshold. Creates a new psg_list. A pseudospectra is selected if it contains peaks, with pval < pval_th }
  \item{fc_th}{Same as pval. Select those groups with contains peaks with fold-change > fc_th. Pval_th and fc_th can be combined}
  \item{sortpval}{Sort diffreport after pvalues}
  \item{...}{Diffreport parameters see \link{diffreport}}
}

\details{
  Batch script wrapper for combining the annotation and the diffreport for a (grouped) xcmsSet \code{xs}.
  Function list:  1: diffreport(), 2: groupFWHM(), 3: findIsotopes(), 4: groupCorr(), 5: findAdducts()
  For a speedup calculation users can create a quick run, with quick = TRUE to preselect pseudospectra of interest.
  The indices of those pseudospectra are set with psg_list in a second run.
  On the other hand, a automatic selection with pval_th and/or fc_th can be performed.
  Returns the normal xcms diffreport table, with the additional CAMERA slots
}
\value{
 \code{annotateDiffreport} returns an diffreport, see \link{diffreport}, within additional columns containing the annotation results. 
}

\examples{
 #Multiple sample 
 library(CAMERA)
 library(faahKO)
 xs.grp     <- group(faahko)
 xs.fill    <- fillPeaks(xs.grp)
 
 #fast preselection
# diffreport  <- annotateDiffreport(xs.fill,quick=TRUE)
# index <- c(1,18,35,45,56) #Make only for those grps a adduct annotation
# diffreport2 <- annotateDiffreport(xs.fill,psg_list=index,metlin = TRUE)

 #automatic selection for groups with peaks p-val < 0.05 and fold-change > 3
# diffreport <- annotateDiffreport(xs.fill,pval_th=0.05,fc=3)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
