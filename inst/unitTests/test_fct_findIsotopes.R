test.findIsotopePspec <- function() {
    ## Test whether the old and new implementation yield the same result.
    library(CAMERA)
    file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
    xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
    an   <- xsAnnotate(xs)
    an   <- groupFWHM(an)
    
    intval <- "maxo"
    object <- an
    ppm <- 5
    maxcharge <- 3
    maxiso <- 4
    mzabs <- 0.01
    minfrac <- 0.5
    filter <- TRUE
    
    isotopeMatrix <- calcIsotopeMatrix(maxiso=maxiso)
    npspectra <- length(object@pspectra)
    
    devppm <- ppm / 1000000;
    params <- list(maxiso=maxiso, maxcharge=maxcharge, devppm=devppm,
                   mzabs=mzabs, IM=isotopeMatrix, minfrac=minfrac, filter=filter)
    ncl <- sum(sapply(object@pspectra, length));
    
    imz  <- object@groupInfo[, "mz", drop=FALSE];
    irt  <- object@groupInfo[, "rt", drop=FALSE];
    mint <- object@groupInfo[, intval, drop=FALSE];      
    
    isotope   <- vector("list", length(imz));
    
    isomatrix <- matrix(ncol=5, nrow=0);
    colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
    
    for(i in seq(along = object@pspectra)){
        ipeak <- object@pspectra[[i]];
        if(length(ipeak) > 1){
            mz  <- imz[ipeak];
            int <- mint[ipeak, , drop=FALSE];
            isomatrix <-  CAMERA:::findIsotopesPspec(isomatrix, mz, ipeak,
                                                     int, params)         
        }
    }
    
    isomatrix_old <- matrix(ncol=5, nrow=0)
    colnames(isomatrix_old) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
    
    for(i in seq(along = object@pspectra)){
        ipeak <- object@pspectra[[i]];
        if(length(ipeak) > 1){
            mz  <- imz[ipeak];
            int <- mint[ipeak, , drop=FALSE];
            isomatrix_old <-  CAMERA:::findIsotopesPspec_orig(
                                           isomatrix_old, mz, ipeak,
                                           int, params)         
        }
    }
    checkEquals(isomatrix, isomatrix_old)
}
