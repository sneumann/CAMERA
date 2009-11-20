test.groupFWHM <- function() {
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 checkEqualsNumeric(nrow(xs@peaks),126)
 checkEqualsNumeric(length(an@spectra),14)
 checkEqualsNumeric(an@spectra[[5]][5],86)
}