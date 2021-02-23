
test.annotateDiffreport.calcIso <- function() {

  library(faahKO)

  xs.group <- group(faahko)
  xs.fill <- fillPeaks(xs.group)

  calcIsoFALSE <- annotateDiffreport(xs.fill, calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE)
  calcIsoTRUE  <- annotateDiffreport(xs.fill, calcIso = TRUE,  calcCiS = TRUE, calcCaS = FALSE)

  checkTrue(
    !identical(calcIsoTRUE[,"pcgroup"], calcIsoFALSE[,"pcgroup"]),
    "The iso was not successfully computed with calcIso = TRUE flag."
  )

}