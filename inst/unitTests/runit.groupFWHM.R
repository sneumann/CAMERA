## single sample
test.anno_single <- function() {
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
    ## groupCorr with  psg_list
         anFCp <- groupCorr(anF, psg_list=c(5,6,7,8,9,10,11,12)) 
        
 ## groupCorr with findIsotopes before
    anI <- findIsotopes(anF)
    anIC <- groupCorr(anI)
    checkEqualsNumeric(length(anIC@pspectra),35)
 ## groupCorr with polarity = "negative"
    anCN <- groupCorr(anF,polarity="negative")
    checkEqualsNumeric(length(anCN@pspectra),41)
 ## findIsotopes and findAdducts
    
    anFI <- findIsotopes(anFC)
    checkEqualsNumeric(nrow(anFI@isoID),27)
    anFA <- findAdducts(anFI, polarity="positive")
    checkEqualsNumeric(length(unique(anFA@annoID[,1])),39)
    ## findAdducts with psg_list
    anFAc <- findAdducts(anFI, polarity="positive", psg_list=c(5,6,7,8,9,10,11,12))
    checkEqualsNumeric(length(unique(anFAc@annoID[,1])),9)
    }

test.anno_multi <- function() {
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
    checkEqualsNumeric(length(xsaC@pspectra),213)
  ##  groupCorr without groupFWHM
    ## manual selection
    xsa <- xsAnnotate(xsg, sample=1)
    xsaC <- groupCorr(xsa)
    checkEqualsNumeric(length(xsaC@pspectra),353)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaC <- groupCorr(xsa)
    checkEqualsNumeric(length(xsaC@pspectra),353)
 ## findIsotopes and findAdducts
    xsaFI <- findIsotopes(xsaC)
    checkEqualsNumeric(nrow(xsaFI@isoID),14)
    xsaFA <- findAdducts(xsaFI, polarity="positive")
    checkEqualsNumeric(length(unique(xsaFA@annoID[,1])),27)
    }