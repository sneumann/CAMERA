\name{xsAnnotate}
\alias{xsAnnotate}
\title{Constructor for a xsAnnotate object which one assigned xcmsSet object}
\description{
  This function handles the construction of xsAnnotate object. It
  is generate from a xcmsSet object adapt its peaktable.
}
\usage{
  xsAnnotate(xs = NULL,sample=NA,category=NA)
}
\arguments{
  \item{xs}{a \code{xcmsSet} object}
  \item{sample}{
    Index of the group xcmsSet sample, that is used for the EIC corelations step.
  }
  \item{category}{
    Class label, for the class label of a grouped xcmsSet with different classes.
  }
}
\examples{
 library(faahKO)
 xs <- group(faahko)
 xs_anno <- xsAnnotate(xs, sample=1, category="WT")
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
