\name{annotateDiffreport}
\alias{annotateDiffreport}
\alias{annotateDiffreport,xsAnnotate-methods}
\title{Deconvolute/Annotate LC/ESI-MS data}
\description{
 Wrapper function for the xcms diffreport and the annotate function. Returns a diffreport with the result from the annotation progress.
}
\usage{
annotateDiffreport(object, sample=NA, sigma = 6, perfwhm = 0.6, cor_eic_th = 0.75, maxcharge = 3, maxiso = 4, ppm = 5, mzabs = 0.01,
                 multiplier = 3, polarity = "positive", nSlaves = 1, psg_list = NULL, pval_th = NULL, fc_th = NULL, quick = FALSE,
                 class1 = levels(sampclass(object))[1], class2 = levels(sampclass(object))[2], filebase = character(),
                 eicmax = 0, eicwidth = 200, sortpval = TRUE, classeic = c(class1, class2), value = c("into", "maxo", "intb"),
                  metlin = FALSE, h = 480, w = 640, ...)
}
\arguments{
  \item{object}{ xcmsSet with peak group assignments}
  \item{sample}{Index of which sample is used for the correlation, NA for auto-selection}
  \item{sigma}{Isotopic peak relationship table}
  \item{perfwhm}{Adduct/Fragment peak relationship table}
  \item{cor_eic_th}{correlation threshold (0..1)}
  \item{maxcharge}{max charge of the ions}
  \item{maxiso}{max number of isotopes}
  \item{ppm}{ppm error}
  \item{mzabs}{absolut error}
  \item{multiplier}{max. number n of [nM+x] clusterions}
  \item{polarity}{Which polarity mode was used for measuring of the ms sample}
  \item{nSlaves}{Number of slaves for parallel calculation (Warning: Beta)}
  \item{psg_list}{Index-list of pseudospectra. Only for those ones will the annotation be calculated. Useful as a speed-up}
  \item{pval_th}{pval threshold. Creates a new psg_list. A pseudospectra is selected if it contains peaks, with pval < pval_th }
  \item{fc_th}{Same as pval. Select those groups with contains peaks with fold-change > fc_th. Pval_th and fc_th can be combined}
  \item{quick}{If quick is TRUE, only groupFWHM and findIsotopes will be calculated. Useful for preselecting groups with psg_list}
  \item{class1}{
    character vector with the first set of sample classes to be
    compared
  }
  \item{class2}{
    character vector with the second set of sample classes to be
    compared
  }
  \item{filebase}{
    base file name to save report, \code{.tsv} file and \code{_eic}
    will be appended to this name for the tabular report and EIC
    directory, respectively. if blank nothing will be saved
  }
  \item{eicmax}{
    number of the most significantly different analytes to create
    EICs for
  }
  \item{eicwidth}{
    width (in seconds) of EICs produced
  }
  \item{sortpval}{
    logical indicating whether the reports should be sorted by
    p-value
  }
  \item{classeic}{
    character vector with the sample classes to include in the EICs
  }
  \item{value}{
    intensity values to be used for the diffreport. \cr
    If \code{value="into"}, integrated peak intensities are used. \cr
    If \code{value="maxo"}, maximum peak intensities are used. \cr
    If \code{value="intb"}, baseline corrected integrated peak intensities are used (only available if peak detection was done by \code{\link{findPeaks.centWave}}).
  }
  \item{metlin}{
    mass uncertainty to use for generating link to Metlin metabolite
    database. the sign of the uncertainty indicates negative or
    positive mode data for M+H or M-H calculation. a value of FALSE
    or 0 removes the column
  }
  \item{h}{
    Numeric variable for the height of the eic and boxplots that are printed out.
  }
  \item{w}{
    Numeric variable for the width of the eic and boxplots print out made.
  }
  \item{...}{
    optional arguments to be passed to \code{\link[multtest]{mt.teststat}} 
  }
}
\details{
  Batch script wrapper for combining the annotation and the diffreport for a (grouped) xcmsSet \code{xs}.
  Function list:  1: diffreport(), 2: groupFWHM(), 3: findIsotopes(), 4: groupCorr(), 5: findAdducts()
  For a speedup calculation users can create a quick run, with quick = TRUE to preselect pseudospectra of interest.
  The indices of those pseudospectra are set with psg_list in a second run.
  On the other hand, a automatic selection with pval_th and/or fc_th can be performed.
  Returns the normal xcms diffreport table, with the additional CAMERA slots
}
\value{
\code{annotate}	returns an xsAnnotate object. For more information about see \link{xsAnnotate-class}
}

\examples{
 #Multiple sample 
 library(CAMERA)
 library(faahKO)
 xs.grp     <- group(faahko)
 xs.fill    <- fillPeaks(xs.grp)
 
 #fast preselection
 diffreport  <- annotateDiffreport(xs.fill,quick=TRUE)
 index <- c(1,18,35,45,56) #Make only for those grps a adduct annotation
 diffreport2 <- annotateDiffreport(xs.fill,psg_list=index)

 #automatic selection, all groups with peaks p-val < 0.05 and fold-change > 3
 diffreport <- annotateDiffreport(xs.fill,pval_th=0.05,fc=3)
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
