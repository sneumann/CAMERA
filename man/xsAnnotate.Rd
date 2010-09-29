\name{xsAnnotate}
\alias{xsAnnotate}
\title{Constructor for a xsAnnotate object which one assigned xcmsSet object}
\description{
  This function handles the construction of xsAnnotate object. It extract the peaktable from a given xcmsSet, which is used for all further analysis.
  The xcmsSet can be a single sample or multiple sample experiment. Since some functions needs to go back into the raw data, the selection algorithm must be choosen for a multiple sample.
  CAMERA includes two different strategies: A fixed selection, which sample = index of sample, or the default automatic solution (sample = NA). 
  The automatic solution chooses the best sample for a specifc groups called pseudospectrum, see groupFWHM and groupCorr.
}
\usage{
  xsAnnotate(xs = NULL,sample=NA, nSlaves = 1)
}
\arguments{
  \item{xs}{a \code{xcmsSet} object}
  \item{sample}{
    Index of the group xcmsSet sample, that is used for the EIC corelations step. For automatic selection don''t set a value.
  }
  \item{nSlaves}{
    For parallel mode set nSlaves higher than 1, but not higher than the number of cpu cores.
  }
}
\examples{
 library(faahKO)
 xs <- group(faahko)
 xs.anno <- xsAnnotate(xs, sample=1)
 #With automatic selection
 xs.anno.autoselect <- xsAnnotate(xs)
}
\value{
  A \code{xsAnnotate} object.
}
\seealso{
  \code{\link{xsAnnotate-class}}
}
\author{Carsten Kuhl,
 \email{ckuhl@ipb-halle.de}}
\keyword{file}
