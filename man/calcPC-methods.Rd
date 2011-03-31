\name{calcPC-methods}
\docType{methods}
\alias{calcPC,xsAnnotate-method}
\alias{calcPC}
\alias{calcPC,calcPC-method}
\title{Peakclustering into pseudospectra according to a distance matrix }
\description{
  A number of clustering methods exist in CAMERA. \code{calcPC}
  is the generic method.
}

\usage{
  calcPC(object, method, ...)
}

\arguments{
  \item{object}{\code{\link{xsAnnotate-class}} object}
  \item{method}{Method to use for clustering. See details.}
  \item{...}{Optional arguments to be passed along}
}
\details{
  This algorithms cluster peaks from a xsAnnotate object into pseudospectra
  according to a provided distance matrix. Therefore all peaks are transformend
  into a graph, with peaks as nodes and the value from the distance matrix as edges.
  Afterwards a graph separation algorithm is applied, which searches in the graph for
  clusters. See the manpages of the specific clustering algorithms for more information.

  If the xsAnnotate is pregrouped, for example \link{groupFWHM},
  only the already existing groups will be further processed.

  The different algorithms that can be used by specifying them with the
  \code{method} argument. For example to use the highly connected 
  subgraphs approach by E. Hartuv, R. Shamir, (1999), one would use:
  \code{calcPC(object, method="hcs")}. This is also
  the default, see \code{\link{calcPC.hcs}}.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}, which are most likely ajc.
  The parameter ajc is the peak distance matrix.
  
  \code{getOption("BioC")$CAMERA$findPeaks.methods} returns
  a character vector of \emph{nicknames} for the 
  algorithms available. 
  
  The function returns a xsAnnotate object with grouping information,
  as list of peak indices. They are stored as object@pspectra.
}

\seealso{
  \code{\link{calcPC.lpc}}
  \code{\link{calcPC.hcs}}
  \code{\link{xsAnnotate-class}}
}
\keyword{methods}
