## single sample
test.groupFWHM_CORR <- function() {
    file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
    an   <- xsAnnotate(xs)
    anF   <- groupFWHM(an)
    checkEqualsNumeric(nrow(xs@peaks),126)
    checkEqualsNumeric(length(anF@pspectra),14)
    checkEqualsNumeric(anF@pspectra[[5]][5],86)
    anC <- groupCorr(an) ## groupCORR without groupFWHM
    checkEqualsNumeric(length(anC@pspectra),55)
    anFC <- groupCorr(anF) ## groupCORR with groupFWHM
    checkEqualsNumeric(length(anFC@pspectra),43)

 ## groupCorr with findIsotopes before
    anI <- findIsotopes(an)
    anIC <- groupCorr(anI)

    }

test.groupFWHM_CORR_multi <- function() {
    library(faahKO)         
    filepath <- system.file("cdf", package = "faahKO")
    xsg <- group(faahko)
  ##  groupCorr after groupFWHM
    ## manual selection
    xsa <- xsAnnotate(xsg, sample=1)
    xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
    xsaC <- groupCorr(xsaF)
    checkEqualsNumeric(length(xsaC@pspectra),275)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
    xsaC <- groupCorr(xsaF)
    checkEqualsNumeric(length(xsaC@pspectra),344)
 ## groupCorr with polarity = "negative"
    xsaCN <- groupCorr(xsaF,polarity="negative")
    checkEqualsNumeric(length(xsaCN@pspectra),344)

  ##  groupCorr without groupFWHM
    ## manual selection
    xsa <- xsAnnotate(xsg, sample=1)
    xsaC <- groupCorr(xsa)
    checkEqualsNumeric(length(xsaC@pspectra),353)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaC <- groupCorr(xsa)
    ## checkEqualsNumeric(length(xsaC@pspectra),344) --> erzeugt blos ein pspectrum .. wieso?
    }