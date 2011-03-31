\name{calcPC.lpc}
\docType{methods}
\alias{calcPC.lpc}
\alias{calcPC.lpc,xsAnnotate-method}
\title{Peakclustering into pseudospectra with the label-propagation-community algorithm}
\description{
  Cluster peaks from an xsAnnotate object into pseudospectra
}
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{
    calcPC.lpc(object, ajc=NULL, psg_list=NULL)
    }
}
}}
\details{
  In some cases, is the peak grouping after retentiontime with \code{\link{groupFWHM}} not enough to
  separate co-elution compounds. Therefore \code{\link{groupCorr}} use additional correlation analysis to achieve a separation.
  \code{calcPC} is part of this approach, which takes the calculated weighted edge list and performs the graph clustering.
  It returns an xsAnnotate object with further separated pseudospectra. 
}
\arguments{
  \item{object}{\code{xsAnnotate} object}
  \item{ajc}{Weighted symbolic edge list as four column matrix ("x","y","cor","ps"). Columns x,y are peak indices, cor the edge value
  and ps the pseudospectrum index, where both peaks occur.
  }
  \item{psg_list}{
    additional vector ps pseudospectra indices, which are used in the clustering. If set to NULL all pseudospectra will be processed.
  }
}
\author{Carsten Kuhl, \email{ckuhl@ipb-halle.de}}

\seealso{
  \code{\link{calcPC}}
  \code{\link{groupCorr}}
  \code{\link{xsAnnotate-class}}
  \code{\link{label.propagation.community}}
}
\keyword{methods}
