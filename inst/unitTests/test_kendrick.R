test.grp.lpc <- function() {

  library(RUnit)
  library(faahKO)
  library(CAMERA)

  xs   <- group(faahko)
  
  ## With specific selected sample
  xsa     <- xsAnnotate(xs)

  ## Screen for substance with CH2 mass gain and RT increase
  k1 <- findKendrickMasses(xsa, masses=c(14, 14.01565),
                           maxHomologue=1, error=0.05, time=60,
                           plot=FALSE)

  checkEqualsNumeric(nrow(k1), 22)

  k2 <- findKendrickMasses(xsa, masses=c(14, 14.01565),
                           maxHomologue=1, error=0.05, time=-60,
                           plot=TRUE)

  checkEqualsNumeric(nrow(k2), 22)
  
}
