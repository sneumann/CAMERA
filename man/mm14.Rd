\name{mm14}
\docType{data}
\alias{mm14}
\title{Extract of marker mixture 14 LC/MS data}
\description{
 xcmsSet object containing quantitated LC/MS peaks from a marker mixture.
 The data is a centroided subset from 117-650 m/z and 271-302 seconds with 134 peaks.
 Positive ionization mode data in mzData file format.
}
\usage{data(mm14)}
\format{
The format is: \preformatted{
Formal class 'xcmsSet' [package "xcms"] with 8 slots
    @ peaks    : num [1:83, 1:11] 117 117 118 119 136
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt"
  ..@ groups   : logi[0 , 0 ]
  ..@ groupidx : list()
  ..@ phenoData:'data.frame':	1 obs. of  1 variable:
  .. ..$ class: Factor w/ 1 level "mzdata": 1
  ..@ rt       :List of 2
  .. ..$ raw      :List of 1
  .. .. ..$ : num [1:112] 270 271 271 271 272 ...
  .. ..$ corrected:List of 1
  .. .. ..$ : num [1:112] 270 271 271 271 272 ...
  ..@ filepaths: chr "mzdata/MM14.mzdata"
  ..@ profinfo :List of 2
  .. ..$ method: chr "bin"
  .. ..$ step  : num 0.1
  ..@ polarity : chr(0)
 }}
\details{
 The corresponding raw mzData files are located in the mzData subdirectory of this package.
}
\source{http://doi:10.1186/1471-2105-9-504}
\references{
 Data originally reported in "Highly sensitive feature detection for high resolution LC/MS"
 BMC Bioinformatics; 2008; 9:504.
}
\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{datasets}
