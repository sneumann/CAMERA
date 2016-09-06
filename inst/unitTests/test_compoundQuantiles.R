test.compoundQuantiles <- function() {
  #library(TODO)
  
  ## libraries
  checkTrue(length(compoundLibraries()) > 0, "There is at least one compound library supported")
  checkTrue(length(massWindowSizes(compoundLibraries()[[1]])) > 0, "There is at least one mass window supported")
  
  ## instantiate
  o <- compoundQuantiles()
  
  ## check general properties
  checkTrue(length(o@elementSet) > 0, "There is at least one element supported")
  checkTrue(length(o@isotopeSet) > 0, "There is at least one isotope supported")
  checkTrue((o@maxCompoundMass - o@minCompoundMass) > 0, "There is a mass interval supported")
  checkTrue(length(o@quantileSet) > 0, "There is at least one quantile supported")
  
  ## examples
  compoundMass <- 503
  quantileLow  <- 0.05
  quantileHigh <- 0.95
  
  ## example for element count
  element <- "C"
  countLow  <- getAtomCount(object = o, element = element, mass = compoundMass, quantile = quantileLow)
  countHigh <- getAtomCount(object = o, element = element, mass = compoundMass, quantile = quantileHigh)
  
  checkTrue(countLow > 0 & countHigh > 0, "There are positive results")
  
  ## example for isotope proportions
  isotope1 <- 0
  isotope2 <- 1
  propLow  <- getIsotopeProportion(object = o, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileLow)
  propHigh <- getIsotopeProportion(object = o, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileHigh)
  
  checkTrue(propLow > 0 & propHigh > 0, "There are positive results")
}