\name{calcIsotopes-methods}
\docType{methods}
\alias{calcIsotopes}
\alias{calcIsotopes,xsAnnotate-method}
\title{Calculate isotope distance matrix from xsAnnotate object}
\description{
  Processing an xsAnnotate object with annotated isotopes (\link{findIsotopes}). It return a weighted edge list as distance matrix between peaks 
  according to the isotope annotation. The edge value for recognized isotopes is 1 for all cases. The list can be used as input for \code{\link{calcPC}}.
}
\section{Methods}{
\describe{
\item{object = "xsAnnotate"}{
  \code{
    calcIsotopes(object)
    }
}
}}

\arguments{
  \item{object}{\code{xsAnnotate} object}
}
\value{
  A matrix with 4 columns:

  \item{x}{
    peak index
  }
  \item{y}{
    peak index
  }
  \item{cor}{
    edge value, always 1
  }
  \item{ps}{
    pseudospectrum index, which contains x and y
  }
}

\author{Carsten Kuhl, \email{ckuhl@ipb-halle.de}}
\seealso{
  \code{\link{calcPC}}
  \code{\link{xsAnnotate-class}}
}
\keyword{methods}
