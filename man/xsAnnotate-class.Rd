\name{xsAnnotate-class}
\docType{class}
\alias{xsAnnotate-class}
\alias{show,xsAnnotate-method}
\title{Class xsAnnotate, a class for annotate peak data}
\description{
  This class transforms a set of peaks from multiple LC/MS or GC/MS
  samples into a matrix of preprocessed data. It groups the peaks
  and does nonlinear retention time correction without internal
  standards. It fills in missing peak values from raw data. Lastly,
  it generates extracted ion chromatograms for ions of interest.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{xsAnnotate}} constructor
  which gathers peaks from a set NetCDF files. Objects can also be
  created by calls of the form \code{new("xsAnnotate", ...)}.
}

%Slots for class 'xsAnnotate'
%  Code: annoGrp annoID category derivativeIons formula grp_info isoID
%        isotopes peaks polarity pspectra ruleset sample xcmsSet


\section{Slots}{
  \describe{
    \item{\code{annoGrp}:}{Assignment of mass hypotheses to correlation groups}
    \item{\code{annoID}:}{The assignemnt of peaks to the mass difference
  rule used.}
    \item{\code{category}:}{A single string category}
    \item{\code{grp_info}:}{grp_info}
    \item{\code{isoID}:}{isoID}
    \item{\code{polarity}:}{A single string with the polarity mode of
  the peaks}
    \item{\code{ruleset}:}{A dataframe describing the
  mass difference rules used for the annotion}
    \item{\code{xcmsSet}:}{The embedded xcmsSet, from which the
  annotation was created}
    \item{\code{peaks}:}{
      matrix containing peak data
    }
\item{\code{pspectra}:}{
      list containing peak data
    }
\item{\code{isotopes}:}{
      matrix containing peak data
    }
\item{\code{derivativeIons}:}{
      matrix containing peak data
    }
\item{\code{formula}:}{
      matrix containing peak data
    }
\item{\code{sample}:}{
      the number of the used xcmsSet sample
    }
  }
}
\section{Methods}{
  \describe{
    \item{\link{groupFWHM}}{
      \code{signature(object = "xsAnnotate")}: group the peak data after the FWHM of the retention time
    }
    \item{\link{groupCorr}}{
      \code{signature(object = "xsAnnotate")}: group the peak data after the Correlation of the EIC
    }
    \item{\link{findIsotopes}}{
      \code{signature(object = "xsAnnotate")}: search for possible isotopes in the spectra
    }
    \item{\link{findAdducts}}{
      \code{signature(object = "xsAnnotate")}: search for possible
      adducts in the spectra
    }
    \item{\link{plotEICs}}{
      \code{signature(object = "xsAnnotate")}: plot EICs of pseudospectra 
    }
  }
}
\author{Carsten Kuhl, \email{ckuhl@ipb-halle.de}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{xsAnnotate}}
}
\keyword{classes}
