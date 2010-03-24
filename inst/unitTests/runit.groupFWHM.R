## single sample
test.groupFWHM_CORR <- function() {
 library(faahKO)         
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 anF   <- groupFWHM(an)
 checkEqualsNumeric(nrow(xs@peaks),126)
 checkEqualsNumeric(length(an@pspectra),14)
 checkEqualsNumeric(an@pspectra[[5]][5],86)
 anC <- groupCorr(an) ## check to be continued tomorrow here
}


library(faahKO)         
filepath <- system.file("cdf", package = "faahKO")
xsg <- group(faahko)
## manual selection
xsa <- xsAnnotate(xsg, sample=1)
xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
xsaC <- groupCorr(xsaF)
## highestPeak-selection
xsa <- xsAnnotate(xsg, sample=NA)
xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
xsaC <- groupCorr(xsaF)