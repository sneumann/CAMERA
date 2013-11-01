## single sample MPI
test.anno_single <- function() {
    if(require("Rmpi", quietly=TRUE)){
    file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
    an   <- xsAnnotate(xs,nSlaves=2)
    anF  <- groupFWHM(an)
    anI  <- findIsotopes(anF)
    file  <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
    rules <- read.csv(file)
    anFA <- findAdducts(anI,polarity="positive",rules=rules)
    checkEqualsNumeric(length(unique(anFA@annoID[,1])),30)    
  }
}

test.anno_multi <- function() {
    library(faahKO)
    if(require("Rmpi", quietly=TRUE)){         
      filepath <- system.file("cdf", package = "faahKO")
      xsg <- group(faahko)

      file  <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
      rules <- read.csv(file)

      xsa <- xsAnnotate(xsg, sample=1,nSlaves=2)
      xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
      xsaC <- groupCorr(xsaF)
      xsaFI <- findIsotopes(xsaC)
      xsaFA <- findAdducts(xsaFI, polarity="positive",rules=rules)
      checkEqualsNumeric(length(unique(xsaFA@annoID[,1])),20)

      xsa <- xsAnnotate(xsg, sample=c(1:8),nSlaves=2)
      xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
      xsaC <- groupCorr(xsaF)
      xsaFI <- findIsotopes(xsaC)
      xsaFA <- findAdducts(xsaFI, polarity="positive",rules=rules)
      checkEqualsNumeric(length(unique(xsaFA@annoID[,1])),28)

      xsa <- xsAnnotate(xsg, sample=NA,nSlaves=2)
      xsaF <- groupFWHM(xsa, sigma=6, perfwhm=0.6)
      xsaC <- groupCorr(xsaF)
      xsaFI <- findIsotopes(xsaC)
      xsaFA <- findAdducts(xsaFI, polarity="positive",rules=rules)
      checkEqualsNumeric(length(unique(xsaFA@annoID[,1])),16)
    }
}
