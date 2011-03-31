\name{psDist-methods}
\docType{methods}
\alias{psDist-methods}
\alias{psDist}
\alias{psDist,CAMERA-method}
\title{Distance methods for xsAnnotate}
\description{
  The package xcms contains several methods for calculating a distance between two sets of peaks. the CAMERA method \code{psDist} is the generic wrapper to use these methods for processing two pseudospectra from two different xsAnnotate objects.
 }
\section{Methods}{
	\describe{
	\item{object1 = "xsAnnotate"}{
		\code{
		specDist(object1, object2, pspectrum1, pspectrum2, 
        method,...)
		}
		}
	}
}
\arguments{
  \item{object1}{a xsAnnotate object with pseudospectra}
  \item{object2}{a xsAnnotate object with pseudospectra}
  \item{PSpec1}{index of pseudospectrum in object1}
  \item{PSpec2}{index of pseudospectrum in object2}
  \item{method}{method to use for distance calculation. See details.}
  \item{...}{mzabs, mzppm and parameters for the distance function.}
}
\details{
  Different algorithms can be used by specifying them with the
  \code{method} argument. For example to use the "meanMZmatch"
  approach one would use:
  \code{specDist(object1, object2, pspectrum1, pspectrum2, method="meanMZmatch")}. This is also
  the default.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}. 
  
  A character vector of \emph{nicknames} for all the 
  algorithms which are available is returned by 
  \code{getOption("BioC")$xcms$specDist.methods}.
  If the nickname of a method is called "meanMZmatch",
  the help page for that specific method can
  be accessed with \code{?specDist.meanMZmatch}. 
}
\value{
\item{mzabs}{
    maximum absolute deviation for two matching peaks
  }
  \item{mzppm}{
    relative deviations in ppm for two matching peaks
  }
  \item{symmetric}{
    use symmetric pairwise m/z-matches only, or each match
    }
}

\author{Joachim Kutzera, \email{jkutzer@ipb-halle.de}}
\keyword{methods}