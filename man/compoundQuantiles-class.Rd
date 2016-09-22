\name{compoundQuantiles-class}
\docType{class}
\alias{compoundQuantiles-class}
\title{Class compoundQuantiles encapsulates compound statistics from different databases.}
\description{
  The user is able to get the expected number of atoms of element e (C, N, ...) 
  for a compound of mass m for a q-quantile.
  I.e. getAtomCount(object = compoundQuantiles(), element = e, mass = m, quantile = q) returns the number of atoms
  of element e in a compound of mass m in the lowest-(q*100)% of a sorted set of compounds
  (sorted ascending by the possible number of atoms of element e for compounds of such mass).
  
  The user is able to get the expected proportion between the intensities of two isotope peaks
  for a compound of mass m for a q-quantile.
  I.e. getIsotopeProportion(object = compoundQuantiles(), isotope1 = i1, isotope2 = i2, mass = m, quantile = q) returns the
  isotope proportion i1 / i2 for a compound of mass m in the lowest-(q*100)% of a sorted set of compounds
  (sorted ascending by the possible isotope proportions for compounds of such mass).
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{compoundQuantiles}} constructor.
}

\section{Slots}{
  \describe{
    \item{\code{compoundLibrary}:}{The compound library to rely on (kegg, chebi, ...)}
    \item{\code{massWindowSize}:}{The mass window size of the compound statistics (25, 100, ...)}
    \item{\code{minCompoundMass}:}{Minimum compound mass for which there are statistics}
    \item{\code{maxCompoundMass}:}{Maximum compound mass for which there are statistics}
    \item{\code{numberOfMassWindows}:}{Number of mass windows}
    \item{\code{numberOfIsotopes}:}{Number of isotopes for which there are isotope ratio quantiles}
    \item{\code{isotopeSet}:}{The set of isotopes for which there are isotope ratio quantiles}
    \item{\code{elementSet}:}{The set of elements for which there are element count statistics}
    \item{\code{quantileSet}:}{The set of quantiles for which there are isotope ratio statistics}
    \item{\code{eleCounters_e_q_mw}:}{Three dimensional array containing the element count statistics (element, quantile, mass window index)}
    \item{\code{proportions_i_q_mw}:}{Three dimensional array containing the isotope ratio quantiles relative to the monoisotopic peak (isotope index, quantile, mass window index)}
  }
}

\section{Methods}{
  \describe{
    \item{\link{getAtomCount}}{
      \code{signature(object = "xsAnnotate")}: returns the number of atoms of the specified element for the given quantile and mass window index
    }
    \item{\link{getIsotopeProportion,compoundQuantiles-method}}{
      \code{signature(object = "xsAnnotate")}: returns the isotope ratio of the specified isotope for the given quantile and mass window index relative to the monoisotopic peak
    }
  }
}
\author{Hendrik Treutler, \email{hendrik.treutler@ipb-halle.de}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{compoundQuantiles}}
  \code{\link{getAtomCount}}
  \code{\link{getIsotopeProportion}}
}
\keyword{classes}
