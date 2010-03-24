test.psDist <- function() {
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 checkEqualsNumeric(CAMERA:::psDist(an,an,PSpec1=10,PSpec2=10, method="cosine", nPmin=1),0.75)
 checkEqualsNumeric(CAMERA:::psDist(an,an,PSpec1=10,PSpec2=10, method="meanMZmatch"),1)
 }