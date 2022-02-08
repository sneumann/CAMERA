\name{groupCorr}
\docType{methods}
\alias{groupCorr}
\alias{groupCorr,xsAnnotate-method}
\title{EIC correlation grouping of LC/ESI-MS data}
\description{
 Peak grouping after correlation information into pseudospectrum groups for an xsAnnotate object. Return an xsAnnotate object with grouping information.
}
\usage{
  groupCorr(object,cor_eic_th=0.75, pval=0.05, graphMethod="hcs",
  calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE, psg_list=NULL, xraw=NULL, 
  cor_exp_th=0.75, intval="into", ...)
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{cor_eic_th}{Correlation threshold for EIC correlation}
  \item{pval}{p-value threshold for testing correlation of significance}
  \item{graphMethod}{Clustering method for resulting correlation graph. See \link{calcPC} for more details.}
  \item{calcIso}{Include isotope detection informationen for graph clustering}
  \item{calcCiS}{Calculate correlation inside samples}
  \item{calcCaS}{Calculate correlation accross samples}
  \item{psg_list}{Vector of pseudospectra indices. The correlation analysis will be only done for those groups}
  \item{xraw}{Optional xcmsRaw object, which should be used for raw data extraction}
  \item{cor_exp_th}{Threshold for intensity correlations across samples}
  \item{intval}{Selection of the intensity values (such as "into") that should be used in the correlation analysis.
  See \code{\link{getPeaklist}} for all allowed values.}
  \item{...}{Additional parameter}
}
\details{
  The algorithm calculates different informations for group peaks into so called pseudospectra. This pseudospectra contains peaks, with have a high correlation between each other.
  So far three different kind of information are available. Correlation of intensities across samples (need more than 3 samples), EIC correlation between peaks inside a sample and
  additional the informationen about recognized isotope cluster can be included. After calculation of all these informations, they are combined as edge value into a graph object.
  A following graph clustering algorithm separate the peaks (nodes in the graph) into the pseudospectra.
}
\examples{
 library(CAMERA)
 file        <- system.file('mzML/MM14.mzML', package = "CAMERA");
 xs          <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5, 10));
 an          <- xsAnnotate(xs);
 an.group    <- groupFWHM(an);
 an.iso      <- findIsotopes(an.group); #optional step for using isotope information
 an.grp.corr <- groupCorr(an.iso, calcIso=TRUE);
 
 #For csv output
 # write.csv(file="peaklist_with_isotopes.csv",getPeaklist(an))

 #Multiple sample 
 library(faahKO)
 xs.grp       <- group(faahko)
 
 #With selected sample
 xsa          <- xsAnnotate(xs.grp, sample=1)
 xsa.group    <- groupFWHM(xsa)
 xsa.iso      <- findIsotopes(xsa.group) #optional step
 xsa.grp.corr <- groupCorr(xsa.iso, calcIso=TRUE)

 #With automatic selection
 xsa.auto     <- xsAnnotate(xs.grp)
 xsa.grp      <- groupFWHM(xsa.auto)
 xsa.iso      <- findIsotopes(xsa.grp) #optional step
 index        <- c(1,4) #Only group one and four will be calculate
 #We use also correlation across sample
 xsa.grp.corr <- groupCorr(xsa.iso, psg_list=index, calcIso=TRUE, calcCaS=TRUE)
 #Note: Group 1 and 4 have no subgroups
}

\seealso{
  \code{\link{calcCiS}}
  \code{\link{calcCaS}}
  \code{\link{calcPC}}
  \code{\link{xsAnnotate-class}}
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
