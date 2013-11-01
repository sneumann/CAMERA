\name{xsAnnotate-class}
\docType{class}
\alias{xsAnnotate-class}
\alias{show,xsAnnotate-method}
\title{Class xsAnnotate, a class for annotated peak data}
\description{
  This class transforms a \code{\link{xcmsSet}} object with peaks from multiple LC/MS or GC/MS
  samples into a set of annotation results. It contains searching algorihms for isotopes and adducts,
  peak grouping algorithms to find connected peak, which originate from the same molecule.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{xsAnnotate}} constructor
  which include the peaktable from a provided \code{\link{xcmsSet}}. Objects can also be
  created by calls of the form \code{new("xsAnnotate", ...)}.
}

%Slots for class 'xsAnnotate'
%  Code: annoGrp annoID derivativeIons formula isoID groupInfo
%        isotopes polarity pspectra psSamples ruleset runParallel sample xcmsSet

\section{Slots}{
  \describe{
    \item{\code{annoGrp}:}{Assignment of mass hypotheses to correlation groups}
    \item{\code{annoID}:}{The assignemnt of peaks to the mass difference rule used}
    \item{\code{derivativeIons}:}{List with annotation result for every peak}
    \item{\code{formula}:}{Matrix containing putative sum formula (intended for future use)}
    \item{\code{isoID}:}{Matrix containing IDs and additional of all annotated isotope peaks}
    \item{\code{groupInfo}:}{(grouped) Peaktable with "into" values}
    \item{\code{isotopes}:}{List with annotated isotopid results for every peak}
    \item{\code{polarity}:}{A single string with the polarity mode of the peaks}
    \item{\code{pspectra}:}{List contains all pseudospectra with there peak IDs}
    \item{\code{psSamples}:}{List containing information with sample was sample was selecteted as representative (automatic selection)}
    \item{\code{ruleset}:}{A dataframe describing the mass difference rules used for the annotion}
    \item{\code{runParallel}:}{Flag if CAMERA runs in serial or parallel mode}
    \item{\code{sample}:}{Number of the used xcmsSet sample (beforehand sample selection)}
    \item{\code{xcmsSet}:}{The embedded xcmsSet}
  }
}

\section{Methods}{
  \describe{
    \item{\link{groupFWHM}}{
      \code{signature(object = "xsAnnotate")}: group the peak data after the FWHM of the retention time
    }
    \item{\link{groupCorr}}{
      \code{signature(object = "xsAnnotate")}: group the peak data after the correlation of the EICs
    }
    \item{\link{findIsotopes}}{
      \code{signature(object = "xsAnnotate")}: search for possible isotopes in the spectra
    }
    \item{\link{findAdducts}}{
      \code{signature(object = "xsAnnotate")}: search for possible adducts in the spectra
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
