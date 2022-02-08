## single sample
test.anno_single <- function() {
    file <- system.file('mzML/MM14.mzML', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))

    an   <- xsAnnotate(xs)

    anF   <- groupFWHM(an)

    checkEqualsNumeric(nrow(xs@peaks),126)

    checkEqualsNumeric(length(anF@pspectra),14)

    checkEqualsNumeric(anF@pspectra[[5]][5],86)

    ## groupCORR without groupFWHM
    anC <- groupCorr(an) 
    checkEqualsNumeric(length(anC@pspectra),48)

    ## groupCORR with groupFWHM
    anFC <- groupCorr(anF) 
    checkEqualsNumeric(length(anFC@pspectra),48)

    ## groupCorr with  psg_list
    anFCp <- groupCorr(anF, psg_list=c(5,6,7,8,9,10,11,12)) 
        
    ## groupCorr with findIsotopes before
    anI <- findIsotopes(anF)
    anIC <- groupCorr(anI)
    checkEqualsNumeric(length(anIC@pspectra),48)

        ## findIsotopes without group before
    anI2 <- findIsotopes(an)
    checkEqualsNumeric(nrow(anI2@isoID),32)

    ## findAdducts without anything before
    file  <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
    rules <- read.csv(file)
    anA <- findAdducts(an, polarity="positive",rules=rules)
    checkEqualsNumeric(length(unique(anA@annoID[,1])),55)

    ## findIsotopes and findAdducts 
    anFI <- findIsotopes(anFC)
    checkEqualsNumeric(nrow(anFI@isoID),23)
    anFA <- findAdducts(anFI, polarity="positive")
    checkEqualsNumeric(length(unique(anFA@annoID[,1])),38)

    ## findAdducts with psg_list
    anFAc <- findAdducts(anFI, polarity="positive", psg_list=c(5,6,7,8,9,10,11,12))
    checkEqualsNumeric(length(unique(anFAc@annoID[,1])),7)
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
    #Change from 171 to 159 due to applying correlation accross samples
    checkEqualsNumeric(length(xsaC@pspectra),333)
    #checkEqualsNumeric(length(xsaC@pspectra),171)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
    xsaC <- groupCorr(xsaF)
    #Change from 211 to 236 due to applying correlation accross samples
    checkEqualsNumeric(length(xsaC@pspectra),329)
    #checkEqualsNumeric(length(xsaC@pspectra),211)
  ##  groupCorr without groupFWHM
    ## manual selection
    xsa <- xsAnnotate(xsg, sample=1)
    xsaC <- groupCorr(xsa)
    #Change from 316 to 8 due to applying correlation accros samples
    checkEqualsNumeric(length(xsaC@pspectra),316)
#    checkEqualsNumeric(length(xsaC@pspectra),316)
    ## highestPeak-selection
    xsa <- xsAnnotate(xsg, sample=NA)
    xsaC <- groupCorr(xsa)
    #Change from 211 to 236 due to applying correlation accros samples
    checkEqualsNumeric(length(xsaC@pspectra),316)
    #checkEqualsNumeric(length(xsaC@pspectra),316)
 ## findIsotopes and findAdducts
    xsaFI <- findIsotopes(xsaC)
    #Change from 20 to 109 due to applying correlation accros samples
    #checkEqualsNumeric(nrow(xsaFI@isoID),109)
    #checkEqualsNumeric(nrow(xsaFI@isoID),20)
    file  <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
    rules <- read.csv(file)
    xsaFA <- findAdducts(xsaFI, polarity="positive",rules=rules)
    #Change from 41 to 154 due to applying correlation accros samples
    checkEqualsNumeric(length(unique(xsaFA@annoID[,1])),20)
    }
