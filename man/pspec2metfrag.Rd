\name{pspec2metfrag}
\alias{pspec2metfrag}
\alias{pspec2metfusion}
\title{Export the putative fragments as MetFrag query files}
\description{
  MetFrag is an in-silico metabolite identification system, which aims
  to putatively identify compounds from fragmentation MS data,
  expecially from tandem-MS, but also in-source fragments might
  give additional hints on top of the accurate mass of the precursor
  alone.

}
\usage{
  pspec2metfrag(object, pspecidx=NULL, filedir=NULL) 
  pspec2metfusion(object, pspecidx=NULL, filedir=NULL) 
}
\arguments{
  \item{object}{an xsAnnotate object}
  \item{pspecidx}{Index of pspectra to export, if NULL then all are exported.}
  \item{filedir}{Directory for placement of batch query files}
}

\details{
  For each spectrum in pspecidx (or all in the xsAnnotate object),
  for each [M] mass hypothesis, remove all non-fragment peaks (isotopes,
  clusters, adducts) and pass them to MetFrag and MetFusion batch query files.
}

\value{
  Returns a list 
}

\examples{
library(CAMERA)
file        <- system.file('mzML/MM14.mzML', package = "CAMERA");
xs          <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5, 10));
an          <- xsAnnotate(xs);
an          <- groupFWHM(an);
an          <- findIsotopes(an); #optional step
an          <- findAdducts(an, polarity="positive")

pspec2metfrag(an, pspecidx=c(1))

}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
