
setGeneric("findIsotopesWithValidation", function(object, maxcharge=3, ppm=5, mzabs=0.01, intval=c("maxo","into","intb"), validateIsotopePatterns = TRUE) standardGeneric("findIsotopesWithValidation"));

setMethod("findIsotopesWithValidation", "xsAnnotate", function(object, maxcharge=3, ppm=5, mzabs=0.01, intval=c("maxo","into","intb"), validateIsotopePatterns = TRUE)
{
  #searches in every pseudospectrum mass differences, which match isotope distances
  
  ################################################################################################
  ## sanity checks
  
  ## test maxcharge  
  if(!is.wholenumber(maxcharge) || maxcharge < 1)
    stop("Invalid argument 'maxcharge'. Must be integer and > 0.\n")
  ## test ppm
  if(!is.numeric(ppm) || ppm < 0)
    stop("Invalid argument 'ppm'. Must be numeric and not negative.\n")
  ## test mzabs
  if(!is.numeric(mzabs) || mzabs < 0)
    stop("Invalid argument 'mzabs'. Must be numeric and not negative.\n")
  ## test intval
  intval <- match.arg(intval)
  
  ################################################################################################
  ## init  
  numberOfPS <- length(object@pspectra)
  
  ## scaling
  devppm <- ppm / 1000000
  
  ## Check if object have been preprocessed with groupFWHM
  if(numberOfPS < 1) {
    cat("xsAnnotate contains no pseudospectrum. Regroup all peaks into one!\n")
    numberOfPS <- 1
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo))
    object@psSamples  <- 1
  }
  
  ## number of peaks in pseudospectrum
  numberOfPeaks <- sum(sapply(object@pspectra, length))
  
  ## get mz, rt, and intensity values from peaktable
  
  ## "mz"     "mzmin"  "mzmax"  "rt"     "rtmin"  "rtmax"  "into"   "intb"   "maxo"   "sn"
  ## "egauss" "mu"     "sigma"  "h"      "f"      "dppm"   "scale"  "scpos"  "scmin"  "scmax"  "lmin"   "lmax"   "sample"
  cat("Generating peak matrix!\n")
  mzValues  <- object@groupInfo[, "mz", drop=FALSE]
  rtValues  <- object@groupInfo[, "rt", drop=FALSE]
  
  if(nrow(object@xcmsSet@groups) > 0){
    ## multiple sample or grouped single sample
    if(is.na(object@sample[1])){
      index <- 1:length(object@xcmsSet@filepaths)
    }else{
      index <- object@sample
    }
    #intValues <- groupval(object@xcmsSet,value=intval)[,index,drop=FALSE]
    intValues <- apply(X = groupval(object@xcmsSet,value=intval)[,index,drop=FALSE], MARGIN = 1, FUN = function(x){median(x = x, na.rm = TRUE)})
    snValues  <- unlist(lapply(X = 1:nrow(object@groupInfo), FUN = function(x){median(object@xcmsSet@peaks[object@xcmsSet@groupidx[[x]], "sn"])}))
  }else{
    ## one sample case
    intValues <- object@groupInfo[, intval, drop=FALSE]
    snValues  <- object@groupInfo[, "sn", drop=FALSE]
  }
  
  ## isotope matrix
  isoMatrix <- matrix(ncol=4, nrow=0)
  colnames(isoMatrix) <- c("mpeak", "isopeak", "iso", "charge")
  
  ################################################################################################
  ## find isotopes
  cat("Run isotope peak annotation\n % finished: ")
  
  ## look for isotopes in every pseudospectrum
  numberOfPS <- length(object@pspectra)
  isotopeClusterCounter <- 0
  
  for(psIdx in 1:numberOfPS){
    ## progress
    if(numberOfPS >= 10)      if((psIdx %% round(numberOfPS / 10)) == 0)   cat((psIdx / round(numberOfPS / 10) * 10), ' ')
    
    ## get peak indizes for psIdx-th pseudospectrum
    peakIndeces <- object@pspectra[[psIdx]]
    if(length(peakIndeces) <= 1)
      ## Pseudospectrum has only one peak
      next
    
    ## calculate isotopes
    isoMatrixForPS.list <- findIsotopesForPS(
      peakIndeces, mzValues[peakIndeces], intValues[peakIndeces], snValues[peakIndeces],
      maxcharge=maxcharge, devppm=devppm, mzabs=mzabs, validateIsotopePatterns=validateIsotopePatterns
    )
    
    if(length(isoMatrixForPS.list) == 0)
      ## no isotope cluster found
      next
    
    ## add isotope cluster to isoMatrix
    for(isotopeCluster in 1:length(isoMatrixForPS.list)){
      isoClusterMatrix <- isoMatrixForPS.list[[isotopeCluster]]
      numberOfClusterPeaks <- nrow(isoClusterMatrix)
      isotopeClusterCounter <- isotopeClusterCounter + 1
      
      ## assign cluster index and box
      isoCluster2 <- cbind(
        isoClusterMatrix[, "peak index"], 
        rep(x = isotopeClusterCounter, times = numberOfClusterPeaks), 
        isoClusterMatrix[, "isotope"], 
        isoClusterMatrix[, "charge"]
      )
      isoMatrix <- rbind(isoMatrix, isoCluster2)
    }
  }
  
  ################################################################################################
  ## create isotope matrix for peak list
  isotopeMatrix <- vector(mode="list", length=numberOfPeaks)
  
  for(peakIdx in 1:numberOfPeaks){
    row <- match(peakIdx, isoMatrix[, "mpeak"])
    if(is.na(row))
      ## no isotope found
      next
    
    ## register isotope data
    isotopeMatrix[[peakIdx]] <- list(
      "y"       = isoMatrix[row[[1]], "isopeak"],
      "iso"     = isoMatrix[row[[1]], "iso"] - 1,
      "charge"  = isoMatrix[row[[1]], "charge"],
      "val"     = 0
    )
  }
  
  ## assign isotope matrix to object
  object@isoID <- isoMatrix
  object@isotopes <- isotopeMatrix
  
  ## out
  cat("\nNumber of isotopes:", nrow(object@isoID), "\n")
  return(object)
})

findIsotopesForPS <- function(peakIndeces, mzValues, intValues, snValues, maxcharge, devppm, mzabs, validateIsotopePatterns = TRUE){
  ## peakIndeces  - peak indeces
  ## mzValues     - m/z vector, contains all m/z values from specific pseudospectrum
  ## intValues    - int vector, see above
  ## snValues     - signal-to-noise vector, see above
  ## maxcharge    - maximum allowed charge
  ## devppm       - scaled ppm error
  ## mzabs        - absolut error in m/z
  
  ###################################################################################################
  ## algorithm
  ## 
  ## 1) compute triangular matrix of peak differences regarding m/z
  ## 2) compute for all charges possible isotope chains M, M+1, M+2, ...
  ##    (in case of multiple possible successors take the successor with minimum m/z distance from expected m/z)
  ## 3) select the longest chain and proceed with the remaining peaks
  ## 4) check isotope cluster validity and split chains to chain segments
  ## 5) box results
  ## 
  
  ###################################################################################################
  ## init
  
  ## matrix with all important informationen
  snValues[snValues==0] <- 1
  noiseEstimate <- intValues / snValues
  intensityMin  <- intValues - noiseEstimate
  intensityMax  <- intValues + noiseEstimate
  
  spectrum <- cbind(peakIndeces, mzValues, intValues, intensityMin, intensityMax)
  colnames(spectrum) <- c("peak index", "mz", "intensity", "intensityMin", "intensityMax")
  
  ## order peaks by m/z
  spectrum <- spectrum[order(spectrum[, "mz"]), ]
  numberOfPeaksHere <- nrow(spectrum)
  
  if(numberOfPeaksHere <= 1)
    return(list())
  
  ## calculate allowed m/z errors for the given masses; at least mzabs
  mzErrors <- devppm * spectrum[, "mz"]
  mzErrors[mzErrors < mzabs] <- mzabs
  
  ###################################################################################################
  ## compute m/z difference and possible isotopic connection for every peak pair in pseudospectrum
  isotopeDifference <- 1.0033548378
  expectedDistances <- isotopeDifference / 1:maxcharge
  hitsForCharge_c_p1_p2 <- array(
    dim = c(maxcharge, numberOfPeaksHere - 1, numberOfPeaksHere - 1),
    dimnames = c("charge", "peak1", "peak2")
  )
  
  for(peakIdx in 1:(numberOfPeaksHere - 1)){
    ## create distance matrix
    mzDiff <- spectrum[(peakIdx + 1):numberOfPeaksHere, "mz"] - spectrum[peakIdx, "mz"]
    ## compare distances to expected distances
    mzDiffExp <- outer(-expectedDistances, mzDiff, FUN="+")
    ## create hit matrix
    hits <- abs(mzDiffExp) <= mzErrors[[peakIdx]]
    ## add hit matrix to big picture
    for(chargeIdx in 1:maxcharge)
      hitsForCharge_c_p1_p2[chargeIdx, peakIdx, 1:(numberOfPeaksHere - peakIdx)] <- hits[chargeIdx, ]
  }#end for peakIdx
  
  ###################################################################################################
  ## find and select isotope chains
  goOnSearching <- TRUE
  peakIsProcessed <- rep(x = FALSE, times = numberOfPeaksHere)
  resultChains <- list()
  chainCharges <- list()
  
  while(goOnSearching){
    candidateChains <- list()
    candidateCharge <- list()
    
    ## check all charges and potential starting peaks to built all possible isotope chains
    for(chargeIdx in 1:maxcharge){
      ## get potential monoisotopic peaks to built all possible isotope chains
      potentialStartPeaks <- which(!peakIsProcessed[1:(numberOfPeaksHere - 1)])
      for(peakIdx in potentialStartPeaks){
        ## start new chain with this peak
        peakIdx1 <- peakIdx
        candidateChain <- c(peakIdx1)
        
        ## assemble matching m/z distances to a chain of peaks
        while(
          ## peak is not in any prior chain
          (!peakIsProcessed[peakIdx1]) && 
          ## peak is not the last peak in m/z dimension in the current pseudo spectrum
          (peakIdx1 <= (numberOfPeaksHere - 1)) && 
          ## there are peaks in the right m/z distance
          any(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ], na.rm = TRUE)){
          ## get matching peak index
          matchingIdx <- which(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ])
          
          peakIdx2 <- matchingIdx + peakIdx1
          
          if(length(peakIdx2) > 1){
            ## more than one peak is candidate as successor: take the best matching one regarding m/z
            
            ## get mass of current peak and successor peaks
            mass1   <- spectrum[peakIdx1, "mz"]
            masses2 <- spectrum[peakIdx2, "mz"]
            
            ## take peak with minimum deviation from the expected value
            masses2 <- masses2 - mass1 - expectedDistances[chargeIdx]
            masses2 <- abs(masses2)
            tempIdx <- which.min(masses2)
            peakIdx2 <- peakIdx2[[tempIdx]]
          }
          
          if(peakIsProcessed[peakIdx2])
            ## next peak is already part of a isotope chain --> do not elongate chain
            break
          
          ## elongate chain
          candidateChain <- c(candidateChain, peakIdx2)
          ## set peak for next iteration
          peakIdx1 <- peakIdx2
        }#end of while loop for current peak chain
        if(length(candidateChain) == 1)
          ## a single peak
          next
        
        ## add candidate chain
        candidateChains[[length(candidateChains) + 1]] <- candidateChain
        candidateCharge[[length(candidateCharge) + 1]] <- chargeIdx
      }#end for peakIdx
    }#end for chargeIdx
    
    ## select the longest chain
    if(length(candidateChains) == 0){
      ## no more chains left -> stop searching
      goOnSearching <- FALSE
    } else {
      ## select the longest chain
      maxChainIdx <- which.max(sapply(candidateChains, function(x) length(x)))
      maxChain    <- candidateChains[[maxChainIdx]]
      maxCharge   <- candidateCharge[[maxChainIdx]]
      ## add chain
      resultChains[[length(resultChains) + 1]] <- maxChain
      chainCharges[[length(chainCharges) + 1]] <- maxCharge
      ## mark comprised peaks
      peakIsProcessed[maxChain] <- TRUE
    }
  }#end of while loop for peak chains
  
  if(length(resultChains) == 0)
    ## no isotope cluster found
    return(list())
  
  ###################################################################################################
  ## validate chains
  validatedResultChains <- list()
  validatedChainCharges <- list()
  ## for validation
  cpObj <- compoundQuantiles(compoundLibrary = "kegg")
  maximumIsotopeNumber <- max(cpObj@isotopeSet)
  
  ## available quantiles:
  ## 0.000005 0.999995 0.000010 0.999990 0.000050 0.999950 0.000100 0.999900 0.000500 0.999500 0.001000 0.999000
  ## 0.005000 0.995000 0.010000 0.990000 0.025000 0.975000 0.050000 0.950000 0.100000 0.900000 0.500000
  #quantileLow  <- 0.025
  #quantileHigh <- 0.975
  quantileLow  <- 0.01
  quantileHigh <- 0.99
  
  for(chainIdx in seq_len(length.out = length(resultChains))){
    ## get data
    charge  <- chainCharges[[chainIdx]]
    chain   <- resultChains[[chainIdx]]
    numberOfIsotopes <- length(chain)
    
    mass <- spectrum[[chain[[1]], "mz"]] * charge
    compoundMassInRange <- mass < cpObj@maxCompoundMass
    
    if(validateIsotopePatterns & compoundMassInRange){
      ## validate and decompose chain to segments
      monoisotopicPeakIdx <- 1
      isotopePatternMembership <- integer(length = length(chain))
      isotopePatternLabel <- 1
      isotopePatternMembership[[1]] <- isotopePatternLabel
      for(peakIdx in 2:length(chain)){
        ## check proportion
        #proportion <- spectrum[chain[[monoisotopicPeakIdx]], "intensity"] / spectrum[chain[[peakIdx]], "intensity"]
        proportionObservedMin <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMin"]] / spectrum[[chain[[peakIdx]], "intensityMax"]]
        proportionObservedMax <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMax"]] / spectrum[[chain[[peakIdx]], "intensityMin"]]
        
        isotopeNumber <- peakIdx - monoisotopicPeakIdx
        
        if(isotopeNumber > maximumIsotopeNumber){
          ## pattern is too long to be checkable
          proportionExpectedMin <- 0
          proportionExpectedMax <- Inf
        } else {
          ## fetch expected proportion interval
          proportionExpectedMin <- getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileLow)
          proportionExpectedMax <- getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileHigh)
        }
        
        centerObserved = (proportionObservedMin + proportionObservedMax) / 2;
        centerExpected  = (proportionExpectedMin + proportionExpectedMax) / 2;
        radiusObserved  = (proportionObservedMax - proportionObservedMin) / 2;
        radiusExpected  = (proportionExpectedMax - proportionExpectedMin) / 2;
        isotopeProportionFits <- abs(centerObserved - centerExpected) <= (radiusObserved + radiusExpected)
        #isotopeProportionFits <- proportion >= proportionInterval[[1]] & proportion <= proportionInterval[[2]]
        
        if(!isotopeProportionFits){
          ## isotope proportion does not fit --> start new pattern
          isotopePatternLabel <- isotopePatternLabel + 1
          monoisotopicPeakIdx <- peakIdx
        }
        ## assign pattern membership
        isotopePatternMembership[[peakIdx]] <- isotopePatternLabel
      }## end for peakIdx
      
      if(FALSE & length(unique(isotopePatternMembership)) > 1){
        if(exists("tableMM48"))
          for(labelIdx in unique(isotopePatternMembership)){
            chainPos <- which(isotopePatternMembership == labelIdx)[[1]]
            massHere <- spectrum[[chain[[chainPos]], "mz"]] * charge
            mm48 <- any(abs(tableMM48$Isotope.peak.0.Exact.mass - massHere) <= max(devppm * massHere, mzabs))
            cat(paste(", ", mm48))
          }
        cat(" - ")
        cat(paste(
          chainIdx, 
          numberOfIsotopes, isotopePatternLabel, mass, 
          "[", paste(isotopePatternMembership, collapse = ";"), "]", 
          "[", paste(spectrum[chain, "intensity"], collapse = ";"), "]",
          ifelse(test = isotopePatternLabel > 1, yes = 
                   paste(spectrum[chain[[1]], "intensity"] / spectrum[chain[[(which(diff(isotopePatternMembership) == 1)[[1]]+1)]], "intensity"], "not in",
                         "[", getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = which(diff(isotopePatternMembership) == 1)[[1]], mass = mass, quantile = quantileLow),
                         getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = which(diff(isotopePatternMembership) == 1)[[1]], mass = mass, quantile = quantileHigh), "]"), no = "n/a")
        ), "\n")
      }
      
      ## extract chain segments
      chainSegments <- list()
      for(chainSegmentLabel in 1:isotopePatternLabel){
        chainSegment <- chain[which(isotopePatternMembership == chainSegmentLabel)]
        if(length(chainSegment) <= 1)
          next
        
        chainSegments[[chainSegmentLabel]] <- chainSegment
      }## end for chainSegment
    } else {
      chainSegments <- list()
      chainSegments[[1]] <- chain
    }
    
    ## box chain segments
    for(chainSegment in chainSegments){
      validatedResultChains[[length(validatedResultChains) + 1]] <- chainSegment
      validatedChainCharges[[length(validatedChainCharges) + 1]] <- charge
    }## end for chainSegment
  }## end for chain
  
  ###################################################################################################
  ## assemble list of isotope clusters as result
  isoMatrixForPS.list <- list()
  for(chainIdx in seq_len(length.out = length(validatedResultChains))){
    ## get data
    chain  <- validatedResultChains[[chainIdx]]
    charge <- validatedChainCharges[[chainIdx]]
    
    ## create and fill isotope cluster
    numberOfIsotopes <- length(chain)
    
    isoClusterMatrix <- matrix(nrow = numberOfIsotopes, ncol = 3)
    colnames(isoClusterMatrix) <- c("peak index", "isotope", "charge")
    for(isoIdx in 1:numberOfIsotopes){
      ## translate local peak index to global peak index
      peakIdx <- spectrum[chain[[isoIdx]], "peak index"]
      ## add isotope entry
      isoClusterMatrix[isoIdx, ] <- c(peakIdx, isoIdx, charge)
    }
    ## add isotope cluster to result list
    isoMatrixForPS.list[[length(isoMatrixForPS.list) + 1]] <- isoClusterMatrix
  }## end for chainSegment
  
  return(isoMatrixForPS.list)
}
