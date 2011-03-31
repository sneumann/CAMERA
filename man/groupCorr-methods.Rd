\name{groupCorr-methods}
\docType{methods}
\alias{groupCorr}
\alias{groupCorr,xsAnnotate-method}
\title{EIC correlation grouping of LC/ESI-MS data}
\description{
 Grouping the peaks after the correlation of the EICs into pseudospectragroups for a xsAnnotate object. Return a xsAnnotate object with group information.
}
\usage{
  groupCorr(object,cor_eic_th=0.75, pval=0.05, graphMethod="hcs", calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE, psg_list=NULL, polarity=NA, ...)
}
\arguments{
  \item{object}{The \code{xsAnnotate} object}
  \item{cor_eic_th}{Correlation threshold for the EIC correlation}
  \item{pval}{pvalue for testing correlation of significance}
  \item{graphMethod}{Clustering method for correlation graph}
  \item{calcIso}{Use isotope detection for graph}
  \item{calcCiS}{Calculate correlation in samples for }
  \item{calcCaS}{Calculate correlation accross samples}
  \item{psg_list}{Vector of pseudospectra indices. The correlation analysis will be only done for those groups}
  \item{polarity}{Optional feature, that ensure that [M+H]+, [M+Na]+, [M+K]+ in pos. mode and [M-H]- [M+Cl]- ions  in neg. mode will not be seperated in different groups}
  \item{...}{Additional parameter}
}
\details{
  The algorithm correlates the EIC of a every peak with all others, to find the peaks that belong to one substance. LC/MS data should grouped with groupFWHM first. This step reduce the runtime a lot and increased the number of correct classifications. If you perform the findIsotope step first, annotated isotopes will always be in the same group.
}
\examples{
 library(CAMERA)
 file        <- system.file('mzdata/MM14.mzdata', package = "CAMERA");
 xs          <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5, 10));
 an          <- xsAnnotate(xs);
 an.group    <- groupFWHM(an);
 an.iso      <- findIsotopes(an.group); #optional step
 index       <- c(1,4); #Only group one and four will be calculate
 an.grp.corr <- groupCorr(an.iso, psg_list=index, polarity="positive");
 
 #For csv output
 # write.csv(file="peaklist_with_isotopes.csv",getPeaklist(an))

 #Multiple sample 
 library(faahKO)
 xs.grp       <- group(faahko)
 
 #With selected sample
 xsa          <- xsAnnotate(xs.grp, sample=1)
 xsa.group    <- groupFWHM(xsa)
 xsa.iso      <- findIsotopes(xsa.group) #optional step
 index        <- c(1,4) #Only group one and four will be calculate
 xsa.grp.corr <- groupCorr(xsa.iso, psg_list=index, polarity="positive")

 #With automatic selection
 xsa.auto     <- xsAnnotate(xs.grp)
 xsa.grp      <- groupFWHM(xsa.auto)
 xsa.iso      <- findIsotopes(xsa.grp) #optional step
 index        <- c(1,4) #Only group one and four will be calculate
 xsa.grp.corr <- groupCorr(xsa.iso, psg_list=index, polarity="positive")
 #Note: Group 1 and 4 have no subgroups
}

\seealso{
  \code{\link{calcCiS}}
  \code{\link{calcCaS}}
  \code{\link{xsAnnotate-class}}
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}