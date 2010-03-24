## single sample
test.groupFWHM_CORR <- function() {
    file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
    an   <- xsAnnotate(xs)
    anF   <- groupFWHM(an)
    checkEqualsNumeric(nrow(xs@peaks),126)
    checkEqualsNumeric(length(anF@pspectra),14)
    checkEqualsNumeric(anF@pspectra[[5]][5],86)
 ## groupCORR without groupFWHM
    anC <- groupCorr(an) 
    checkEqualsNumeric(length(anC@pspectra),55)
 ## groupCORR with groupFWHM
    anFC <- groupCorr(anF) 
    checkEqualsNumeric(length(anFC@pspectra),43)
 ## groupCorr with findIsotopes before
    anI <- findIsotopes(anF)
    anIC <- groupCorr(anI)
    checkEqualsNumeric(length(anIC@pspectra),35)
 ## groupCorr with polarity = "negative"
    anCN <- groupCorr(anF,polarity="negative")
    checkEqualsNumeric(length(anCN@pspectra),41)
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
    checkEqualsNumeric(length(xsaC@pspectra),222)
  ##  groupCorr without groupFWHM
    ## manual selection
    xsa <- xsAnnotate(xsg, sample=1)
    xsaC <- groupCorr(xsa)
    checkEqualsNumeric(length(xsaC@pspectra),353)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaC <- groupCorr(xsa)
    checkEqualsNumeric(length(xsaC@pspectra),363)
    }