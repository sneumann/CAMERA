test.psp3c2metfrag <- function() {

file        <- system.file('mzML/MM14.mzML', package = "CAMERA");
xs          <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5, 10));
an          <- xsAnnotate(xs);
an          <- groupFWHM(an);
an          <- findIsotopes(an); #optional step
an          <- findAdducts(an, polarity="positive")

##
## Example call for three pspec
## TODO: 1) guess precursor
##       2) add mode + ppm etc.
##

pspec2metfrag(an, pspecidx=c(3))

pspec2metfrag(an, pspecidx=c(7))

pspec2metfusion(an, pspecidx=c(1,3))

pspec2metfrag(an, pspecidx=c(2))

pspec2metfrag(an)

}

