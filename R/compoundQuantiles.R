#########################################################################################
## Class compoundQuantiles
## 
## The user is able to get the expected number of atoms of element e (C, N, ...) 
## for a compound of mass m for a q-quantile.
## I.e. getAtomCount(object = compoundQuantiles(), element = e, mass = m, quantile = q) returns the number of atoms
## of element e in a compound of mass m in the lowest-(q*100)% of a sorted set of compounds
## (sorted ascending by the possible number of atoms of element e for compounds of such mass).
## 
## The user is able to get the expected proportion between the intensities of two isotope peaks
## for a compound of mass m for a q-quantile.
## I.e. getIsotopeProportion(object = compoundQuantiles(), isotope1 = i1, isotope2 = i2, mass = m, quantile = q) returns the
## isotope proportion i1 / i2 for a compound of mass m in the lowest-(q*100)% of a sorted set of compounds
## (sorted ascending by the possible isotope proportions for compounds of such mass).
## 

#' definition of S4 class "compoundQuantiles"
#' @export
compoundQuantiles <- setClass(
  # Set the name for the class
  "compoundQuantiles",
  
  # Define the slots
  slots = c(
    compoundLibrary     = "character",
    massWindowSize      = "numeric",
    minCompoundMass     = "numeric",
    maxCompoundMass     = "numeric",
    numberOfMassWindows = "integer",
    numberOfIsotopes    = "integer",
    isotopeSet          = "numeric",
    elementSet          = "character",
    quantileSet         = "numeric",
    eleCounters_e_q_mw  = "array",
    proportions_i_q_mw  = "array"
  ),
  
  # Set the default values for the slots. (optional)
  prototype=list(
    compoundLibrary = "NA",
    massWindowSize  = -1,
    minCompoundMass = -1,
    maxCompoundMass = -1,
    numberOfIsotopes    = as.integer(0),
    numberOfMassWindows = as.integer(0),
    isotopeSet  = c(),
    elementSet  = c(),
    quantileSet = c(),
    eleCounters_e_q_mw = array(dim = c(0, 0, 0)),
    proportions_i_q_mw = array(dim = c(0, 0, 0))
  ),
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    if( (object@maxCompoundMass - object@minCompoundMass) / object@massWindowSize != object@numberOfMassWindows ){
      return("Mass bounds variables are inconsistent")
    }
    if(
      length(object@elementSet)   != dim(object@eleCounters_e_q_mw)[[1]] ||
      length(object@quantileSet)  != dim(object@eleCounters_e_q_mw)[[2]] ||
      object@numberOfMassWindows  != dim(object@eleCounters_e_q_mw)[[3]] ||
      object@numberOfIsotopes     != dim(object@proportions_i_q_mw)[[1]] ||
      length(object@quantileSet)  != dim(object@proportions_i_q_mw)[[2]] ||
      object@numberOfMassWindows  != dim(object@proportions_i_q_mw)[[3]]
    ){
      return("Array dimensions are inconsistent.")
    }
    return(TRUE)
  }
)

##' constructor of class compoundQuantiles
##' @title compoundQuantiles constructor
##' @param compoundLibrary the database; see compoundLibraries() for a list of supported databases
##' @param massWindowSize the mass window size for grouping compounds; see massWindowSizes(compoundLibrary = "kegg") for a list of supported databases for e.g. the database kegg
##' @return the compoundQuantiles object
##' @export
##' @exportClass compoundQuantiles
##' @author Hendrik Treutler
##' @examples
##' library(compoundQuantiles)
##' cpObj <- compoundQuantiles(compoundLibrary = "kegg")
compoundQuantiles <- function(compoundLibrary = "kegg", massWindowSize = 50) {
  ######################################################
  ## create new object
  object <- new("compoundQuantiles")
  
  ######################################################
  ## sanity check
  if(is.na(match(compoundLibrary, compoundLibraries())))
    stop(paste("Library '", compoundLibrary, "' is not present. See compoundLibraries() for a full list of supported libraries."))
  if(is.na(match(massWindowSize, massWindowSizes(libraryName = compoundLibrary))))
    stop(paste("Mass-window-size '", massWindowSize, "' is not present. See massWindowSizes() for a full list of supported window sizes."))
  
  ######################################################
  ## fetch data file paths
  lib.loc=.libPaths()
  folder  <- system.file("data", package = "compoundQuantiles", lib.loc=lib.loc)
  
  ## lib path
  regExPattern <- paste("^library_", compoundLibrary, "__maxDa_[0-9]+$", sep = "")
  path <- list.files(path = folder, pattern = regExPattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)[[1]]
  folderName <- tail(x = strsplit(x = path, split = c("/"))[[1]], n = 1)
  
  ## parse max mass
  maxMassTagged  <- strsplit(x = folderName, split = c("__"))[[1]][[2]]
  maxMass        <- strsplit(x = maxMassTagged, split = c("_"))[[1]][[2]]
  
  ## window size path
  regExPattern <- paste("^windowSize_", massWindowSize, "$", sep = "")
  path <- list.files(path = path, pattern = regExPattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)[[1]]
  
  ## files
  regExPatternElementCounts <- "^[a-zA-Z]+_quantiles.tsv$"
  filePathsElementCounts    <- list.files(path = path, pattern = regExPatternElementCounts, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  regExPatternIsotopes <- "^IsotopeQuantiles__iso_[0-9]+.tsv$"
  filePathsIsotopes    <- list.files(path = path, pattern = regExPatternIsotopes, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  ######################################################
  ## parse data files
  
  #########################
  ## element counts
  numberOfFiles <- length(filePathsElementCounts)
  elements_e    <- vector(length = length(filePathsElementCounts), mode = "character")
  dataFrames_e  <- vector(length = length(filePathsElementCounts), mode = "list")
  
  for(elementIdx in 1:numberOfFiles){
    filePath <- filePathsElementCounts[[elementIdx]]
    
    ## get element
    fileName <- tail(strsplit(x = filePath, split = c("/"))[[1]], n = 1)
    element <- strsplit(x = fileName, split = c("_"))[[1]][[1]]
    
    ## get data
    dataFrame <- read.table(filePath, header=FALSE, sep = "\t", as.is=TRUE)
    dataFrame <- dataFrame[, 1:(ncol(dataFrame) - 1)]
    
    ## extract
    rows <- c(2,nrow(dataFrame))
    cols <- c(2,ncol(dataFrame))
    
    quantiles <- as.numeric(dataFrame[rows[[1]]:rows[[2]], 1])
    masses    <- as.numeric(dataFrame[1, cols[[1]]:cols[[2]]])
    data      <- dataFrame[rows[1]:rows[2], cols[1]:cols[2]]
    
    ## annotate
    rownames(data) <- quantiles
    colnames(data) <- masses
    
    ## java to R
    data[data == "Infinity"] <- "Inf"
    data[data == "NaN"] <- "Na"
    
    ## add
    elements_e[[elementIdx]]   <- element
    dataFrames_e[[elementIdx]] <- data
  }
  
  numberOfElements  <- length(elements_e)
  numberOfQuantiles <- length(quantiles)
  numberOfMasses    <- length(masses)
  
  #########################
  ## isotope proportions
  
  ## TODO
  numberOfFiles <- length(filePathsIsotopes)
  isotopes_i    <- vector(length = length(filePathsIsotopes), mode = "numeric")
  dataFrames_i  <- vector(length = length(filePathsIsotopes), mode = "list")
  
  for(isotopeIdx in 1:numberOfFiles){
    filePath <- filePathsIsotopes[[isotopeIdx]]
    
    ## get element
    fileName      <- tail(strsplit(x = filePath, split = c("/"))[[1]], n = 1)
    fileName      <- strsplit(x = fileName, split = c("\\."))[[1]][[1]]
    isotopeTag    <- strsplit(x = fileName, split = c("__"))[[1]][[2]]
    isotopeNumber <- as.numeric(strsplit(x = isotopeTag, split = c("_"))[[1]][[2]])
    
    ## get data
    dataFrame <- read.table(filePath, header=FALSE, sep = "\t", as.is=TRUE)
    dataFrame <- dataFrame[, 1:(ncol(dataFrame) - 1)]
    
    ## extract
    rows <- c(2,nrow(dataFrame))
    cols <- c(2,ncol(dataFrame))
    
    quantiles <- as.numeric(dataFrame[rows[[1]]:rows[[2]], 1])
    masses    <- as.numeric(dataFrame[1, cols[[1]]:cols[[2]]])
    data      <- dataFrame[rows[1]:rows[2], cols[1]:cols[2]]
    
    ## annotate
    rownames(data) <- quantiles
    colnames(data) <- masses
    
    ## add
    isotopes_i[[isotopeIdx]]   <- isotopeNumber
    dataFrames_i[[isotopeIdx]] <- data
  }
  
  numberOfIsotopes <- length(isotopes_i)
  
  ######################################################
  ## initialise
  object@compoundLibrary  = compoundLibrary
  object@massWindowSize   = massWindowSize
  object@minCompoundMass  = masses[[1]]
  object@maxCompoundMass  = masses[[length(masses)]] + object@massWindowSize
  object@numberOfMassWindows = as.integer((object@maxCompoundMass - object@minCompoundMass) / object@massWindowSize)
  object@numberOfIsotopes = numberOfIsotopes
  object@isotopeSet    = isotopes_i
  object@elementSet    = elements_e
  object@quantileSet   = as.numeric(quantiles)
  object@eleCounters_e_q_mw  = array(
    dim = c(numberOfElements, numberOfQuantiles, numberOfMasses), 
    dimnames = list(
      "element"  = as.vector(elements_e, mode = "character"), 
      "quantile" = as.vector(quantiles,  mode = "numeric"), 
      "mass"     = as.vector(masses,     mode = "numeric")
    )
  )
  object@proportions_i_q_mw  = array(
    dim = c(numberOfIsotopes, numberOfQuantiles, numberOfMasses), 
    dimnames = list(
      "isotope"  = as.vector(isotopes_i, mode = "numeric"), 
      "quantile" = as.vector(quantiles,  mode = "numeric"), 
      "mass"     = as.vector(masses,     mode = "numeric")
    )
  )
  
  ## fill data in array
  for(elementIdx in 1:length(dataFrames_e))
    object@eleCounters_e_q_mw[elementIdx, , ] <- data.matrix(dataFrames_e[[elementIdx]])[ , ]
  for(isotopeIdx in 1:length(dataFrames_i))
    object@proportions_i_q_mw[isotopeIdx, , ] <- data.matrix(dataFrames_i[[isotopeIdx]])[ , ]
  
  return(object)
}

##' Returns a set of supported compound databases
##' @title The supported compound databases
##' @return Vector of supported compound databases
##' @export
##' @author Hendrik Treutler
##' @examples
##' library(compoundQuantiles)
##' compoundLibraries()
compoundLibraries <- function() {
  ## get parent folder
  ## TODO
  lib.loc <- .libPaths()
  folder  <- system.file("data", package = "compoundQuantiles", lib.loc=lib.loc)
  #folder <- "/home/htreutle/Data/MsStatistics/Compounds/"
  
  ## get library folders
  regExPattern <- "^library_[a-zA-Z]+__maxDa_[0-9]+$"
  folderNames <- list.files(path = folder, pattern = regExPattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  ## get library names
  firstSplitOfFolderNames  <- strsplit(x = folderNames, split = c("__"))
  librariesTagged <- sapply(firstSplitOfFolderNames, function (x) x[[1]])
  secondSplitOfFolderNames <- strsplit(x = librariesTagged, split = c("_"))
  libraries       <- sapply(secondSplitOfFolderNames, function (x) x[[2]])
  
  return(libraries)
}
##' Returns the set of supported mass window sizes for the given compound database
##' @title The supported mass window sizes
##' @param libraryName The compound database
##' @return Vector of supported mass window sizes
##' @export
##' @author Hendrik Treutler
##' @examples
##' library(compoundQuantiles)
##' massWindowSizes(libraryName = "kegg")
massWindowSizes <- function(libraryName = "kegg") {
  ## get parent folder
  ## TODO
  lib.loc <- .libPaths()
  folder  <- system.file("data", package = "compoundQuantiles", lib.loc=lib.loc)
  #folder <- "/home/htreutle/Data/MsStatistics/Compounds/"
  
  ## get library folders
  regExPattern <- paste("^library_", libraryName, "__maxDa_[0-9]+$", sep = "")
  paths <- list.files(path = folder, pattern = regExPattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  if(length(paths) != 1)
    stop("There is not exactly one folder for library '", libraryName, "'")
  
  path <- paths[[1]]
  regExPattern <- "^windowSize_[0-9]+$"
  folderNames <- list.files(path = path, pattern = regExPattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  ## get library names
  splitOfFolderNames  <- strsplit(x = folderNames, split = c("_"))
  windowSizes         <- sapply(splitOfFolderNames, function (x) x[[2]])
  windowSizes <- sort(as.numeric(windowSizes))
  
  return(windowSizes)
}

setGeneric(
  name="getAtomCount", 
  def=function(object, element, mass, quantile){
    standardGeneric("getAtomCount")
  }
)

##' Returns the number of atoms the specified element in a compound of the specified mass for the specified quantile level
##'
##' @title The number of atoms of the given element
##' @param object A compoundQuantiles object
##' @param element The element of interest specified by element symbol
##' @param mass The mass of the compound specified in atomic units (=dalton)
##' @param quantile The quantile level for the number of atoms
##' @return The number of atoms
##' @export
##' @author Hendrik Treutler
##' @examples
##' library(compoundQuantiles)
##' cpObj <- compoundQuantiles(compoundLibrary = "kegg")
##' 
##' compoundMass <- 503
##' quantileLow   <- 0.05
##' quantileHigh  <- 0.95
##' element <- "C"
##' countLow  <- getAtomCount(object = cpObj, element = element, mass = compoundMass, quantile = quantileLow)
##' countHigh <- getAtomCount(object = cpObj, element = element, mass = compoundMass, quantile = quantileHigh)
##' 
##' print(paste("The ", (quantileHigh - quantileLow) * 100, "% confidence interval for the number of atoms of element ", element, " in a compound with mass ", compoundMass, " is [", countLow, ", ", countHigh, "]", sep = ""))
setMethod(f="getAtomCount", signature="compoundQuantiles", definition=function(object, element, mass, quantile){
  ###################################################################################################################
  ## sanity checks of the given parameters
  if(is.na(match(x = element, table = object@elementSet)))
    stop(paste("Element '", element, "' is not supported! Please refer to 'object@elementSet' for the complete set of supported elements.", sep = ""))
  if(is.na(match(x = quantile, table = object@quantileSet)))
    stop(paste("The quantile level '", quantile, "' is not supported! Please refer to 'object@quantileSet' for the complete set of supported quantile levels.", sep = ""))
  if(mass < object@minCompoundMass || mass > object@maxCompoundMass)
    stop(paste("Mass '", mass, "' is out of range! Please refer to 'object@minCompoundMass' and 'object@maxCompoundMass' for the supported mass range.", sep = ""))
  
  quantile <- as.character(quantile)
  
  ###################################################################################################################
  ## get and return expected atom count
  roundedMass <- object@massWindowSize * floor(mass / object@massWindowSize)
  roundedMass <- as.character(roundedMass)
  
  expectedAtomCount <- object@eleCounters_e_q_mw[
    "element"   = element, 
    "quantile"  = quantile, 
    "mass"      = roundedMass
  ]
  
  return(expectedAtomCount)
})

setGeneric(
  name="getIsotopeProportion", 
  def=function(object, isotope1, isotope2, mass, quantile){
    standardGeneric("getIsotopeProportion")
  }
)

##' Returns the proportion of the intensities of isotope1 versus isotope2 for a compound of the given mass for the given quantile level
##'
##' @title The proportion of the intensities of two isotope peaks
##' @param object A compoundQuantiles object
##' @param isotope1 The divident isotope ranging from 0 (the monoisotopic peak) to 5
##' @param isotope2 The divisor isotope ranging from 0 (the monoisotopic peak) to 5
##' @param mass The mass of the compound specified in atomic units (=dalton)
##' @param quantile The quantile level for the isotope proportion
##' @return The isotope proportion
##' @export
##' @author Hendrik Treutler
##' @examples
##' library(compoundQuantiles)
##' cpObj <- compoundQuantiles(compoundLibrary = "kegg")
##' 
##' compoundMass <- 503
##' isotope1 <- 0
##' isotope2 <- 1
##' quantileLow   <- 0.05
##' quantileHigh  <- 0.95
##' 
##' propLow  <- getIsotopeProportion(object = cpObj, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileLow)
##' propHigh <- getIsotopeProportion(object = cpObj, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileHigh)
##' print(paste("The ", (quantileHigh - quantileLow) * 100, "% confidence interval for the proportion of isotopes ", isotope1, " / ", isotope2, " in a compound with mass ", compoundMass, " is [", propLow, ", ", propHigh, "]", sep = ""))
setMethod(f="getIsotopeProportion", signature="compoundQuantiles", definition=function(object, isotope1, isotope2, mass, quantile){
  ###################################################################################################################
  ## sanity checks of the given parameters
  if(isotope1 != 0 & is.na(match(x = isotope1, table = object@isotopeSet)))
    stop(paste("Isotope1 '", isotope1, "' is not supported! Please refer to 'object@isotopeSet' for the complete set of supported isotopes.", sep = ""))
  if(isotope2 != 0 & is.na(match(x = isotope2, table = object@isotopeSet)))
    stop(paste("Isotope2 '", isotope2, "' is not supported! Please refer to 'object@isotopeSet' for the complete set of supported isotopes.", sep = ""))
  #if(isotope1 == isotope2)
  #  stop(paste("Isotope1 '", isotope1, "' and isotope2 '", isotope2, "' are equal.", sep = ""))
  if(is.na(match(x = quantile, table = object@quantileSet)))
    stop(paste("The quantile level '", quantile, "' is not supported! Please refer to 'object@quantileSet' for the complete set of supported quantile levels.", sep = ""))
  if(mass < object@minCompoundMass || mass > object@maxCompoundMass)
    stop(paste("Mass '", mass, "' is out of range! Please refer to 'object@minCompoundMass' and 'object@maxCompoundMass' for the supported mass range.", sep = ""))
  
  isotope1 <- as.character(isotope1)
  isotope2 <- as.character(isotope2)
  quantile <- as.character(quantile)
  
  ###################################################################################################################
  ## get and return expected isotope proportion
  if(isotope1 == isotope2)
    return(1)
  
  roundedMass <- object@massWindowSize * floor(mass / object@massWindowSize)
  roundedMass <- as.character(roundedMass)
  
  if(isotope1 == 0){
    expectedIsotopeProportion <- object@proportions_i_q_mw[
      "isotope"   = isotope2, 
      "quantile"  = quantile, 
      "mass"      = roundedMass
    ]
    
    return(expectedIsotopeProportion)
  }
  if(isotope2 == 0){
    expectedIsotopeProportion <- object@proportions_i_q_mw[
      "isotope"   = isotope1, 
      "quantile"  = quantile, 
      "mass"      = roundedMass
      ]
    expectedIsotopeProportion <- 1 / expectedIsotopeProportion
    
    return(expectedIsotopeProportion)
  }
  
  ## isotope1 > 0 vs isotope2 > 0
  expectedIsotopeProportion1 <- object@proportions_i_q_mw[
    "isotope"   = isotope1, 
    "quantile"  = quantile, 
    "mass"      = roundedMass
  ]
  expectedIsotopeProportion2 <- object@proportions_i_q_mw[
    "isotope"   = isotope2, 
    "quantile"  = quantile, 
    "mass"      = roundedMass
  ]
  
  expectedIsotopeProportion <- expectedIsotopeProportion2 / expectedIsotopeProportion1
  
  return(expectedIsotopeProportion)
})

##' Runs an example for the usage of package compoundQuantiles
##'
##' @title Runs an example
##' @return void
##' @export
##' @author Hendrik Treutler
##' @examples
##' example()
example <- function() {
  ## attach
  library(compoundQuantiles)
  ## instantiate
  cpObj <- compoundQuantiles(compoundLibrary = "kegg")
  
  ## meta information
  print(paste("Available libraries = {", paste(compoundLibraries(), collapse = ", "), "}", sep = ""))
  print(paste("Compound library = ", cpObj@compoundLibrary, sep = ""))
  print(paste("Available mass window sizes = {", paste(massWindowSizes(cpObj@compoundLibrary), collapse = ", "), "}", sep = ""))
  print(paste("Mass window size = ", cpObj@massWindowSize, sep = ""))
  print(paste("Elements = {", paste(cpObj@elementSet, collapse = ", "), "}", sep = ""))
  print(paste("Isotopes = {", paste(cpObj@isotopeSet, collapse = ", "), "}", sep = ""))
  print(paste("Mass interval = [", cpObj@minCompoundMass, ", ", cpObj@maxCompoundMass, "] Da; #mass windows = ", cpObj@numberOfMassWindows, " a ", cpObj@massWindowSize, "Da", sep = ""))
  print(paste("Quantile levels = {", paste(cpObj@quantileSet, collapse = ", "), "}", sep = ""))
  
  ## examples
  compoundMass <- 503
  quantileLow   <- 0.05
  quantileHigh  <- 0.95
  
  ## example for element count
  element <- "C"
  countLow  <- getAtomCount(object = cpObj, element = element, mass = compoundMass, quantile = quantileLow)
  countHigh <- getAtomCount(object = cpObj, element = element, mass = compoundMass, quantile = quantileHigh)
  
  print(paste("The ", (quantileHigh - quantileLow) * 100, "% confidence interval for the number of atoms of element ", element, " in a compound with mass ", compoundMass, " is [", countLow, ", ", countHigh, "]", sep = ""))
  
  ## example for isotope proportion
  isotope1 <- 0
  isotope2 <- 1
  propLow  <- getIsotopeProportion(object = cpObj, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileLow)
  propHigh <- getIsotopeProportion(object = cpObj, isotope1 = isotope1, isotope2 = isotope2, mass = compoundMass, quantile = quantileHigh)
  
  print(paste("The ", (quantileHigh - quantileLow) * 100, "% confidence interval for the proportion of isotopes ", isotope1, " / ", isotope2, " in a compound with mass ", compoundMass, " is [", propLow, ", ", propHigh, "]", sep = ""))
}
