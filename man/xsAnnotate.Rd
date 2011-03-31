\name{xsAnnotate}
\alias{xsAnnotate}
\title{xsAnnotate constructor for an provided xcmsSet object}
\description{
  This function deals with the construction of an xsAnnotate object. It extracts the peaktable from a provided xcmsSet, which is used for all further analysis.
  The xcmsSet can be a single sample or multiple sample experiment. Since some functions needs the raw data a selection algorithm must be choosen in the case of a multiple sample.
  CAMERA includes two different strategies: A defined selection of samples (sample = indices of samples) or the default automatic solution (sample = NA). 
  The automatic solution chooses the best sample for a specifc groups called pseudospectrum, see \link{groupFWHM} and \link{groupCorr}.
  It returns a xsAnnotate object, see \link{xsAnnotate-class}.
}
\usage{
  xsAnnotate(xs = NULL, sample=NA, nSlaves = 1)
}
\arguments{
  \item{xs}{a \code{xcmsSet} object}
  \item{sample}{
    Indices of the group xcmsSet sample, that are used for the EIC correlation step. For automatic selection don't set a value.
    For use all samples simply define sample = c(1:n), with n = number of samples.
  }
  \item{nSlaves}{
    For parallel mode set nSlaves higher than 1, but not higher than the number of cpu cores.
  }
}
\examples{
 library(faahKO)
 xs <- group(faahko)
 xsa <- xsAnnotate(xs, sample=c(1:12))

 #With automatic selection
 xsa.autoselect <- xsAnnotate(xs)
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
