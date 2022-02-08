library(CAMERA)

test.grp.lpc <- function() {
    file <- system.file('mzML/MM14.mzML', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))

    an   <- xsAnnotate(xs)

    anF   <- groupFWHM(an)

    checkEqualsNumeric(nrow(xs@peaks),126,"xset contains 126 Peaks")

    checkEqualsNumeric(length(anF@pspectra),14,"xsa has 14 compound spectra")

    checkEqualsNumeric(anF@pspectra[[5]][5],86,
                       "5th peak in 5th compound spectra is number 86")

    ## groupCORR without groupFWHM
    #anC <- groupCorr(an,graphMethod="lpc") 
    #checkEqualsNumeric(length(anC@pspectra),18,
    #                   " lpc created 18 compound spectra")

    ## groupCORR with groupFWHM
    
    anFC <- groupCorr(anF,graphMethod="lpc") 
    checkEqualsNumeric(length(anFC@pspectra),35,
                       " lpc created 35 compound spectra")
    
    }

test.grpDen <- function() {
  file <- system.file('mzML/MM14.mzML', package = "CAMERA")
  xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
  
  xsa   <- xsAnnotate(xs)
  
  xsa.grp   <- groupDen(xsa, bw=0.4)
  
  checkEqualsNumeric(nrow(xs@peaks),126)
  
  checkEqualsNumeric(length(xsa.grp@pspectra),11)
  
  checkEqualsNumeric(xsa.grp@pspectra[[5]][5],100)
  
}

