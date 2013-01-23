\name{pspec2metfrag}
\alias{pspec2metfrag}
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
}
\arguments{
  \item{object}{an xsAnnotate object}
  \item{pspecidx}{Index of pspectra to export, if NULL then all are exported.}
  \item{filedir}{Directory for placement of batch query files}
}

\details{
  For each spectrum,

  for each mass hypothesis,

  remove all non-fragments
  
}

\value{
  Returns XXX
}

\examples{
 \dontrun{
 }
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
