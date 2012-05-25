##Functions for findIsotopes
calcIsotopeMatrix <- function(maxiso=4){
  
  if(!is.numeric(maxiso)){
    stop("Parameter maxiso is not numeric!\n")  
  } else if(maxiso < 1 | maxiso > 8){
    stop(paste("Parameter maxiso must between 1 and 8. ",
          "Otherwise use your own IsotopeMatrix.\n"),sep="")
  }
  
  isotopeMatrix <- matrix(NA, 8, 4);
  colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
  
  isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
  isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
  isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
  isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
  isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
  isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
  isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
  isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)  
  
  return(isotopeMatrix[1:maxiso, , drop=FALSE])

}

findIsotopesPspec <- function(isomatrix, mz, ipeak, int, params){
  #isomatrix - isotope annotations (5 column matrix)
  #mz - m/z vector, contains all m/z values from specific pseudospectrum
  #int - int vector, see above
  #maxiso - how many isotopic peaks are allowed
  #maxcharge - maximum allowed charge
  #devppm - scaled ppm error
  #mzabs - absolut error in m/z
  
  #matrix with all important informationen
  spectra <- matrix(c(mz, ipeak), ncol=2)
  int     <- int[order(spectra[, 1]), , drop=FALSE]
  spectra <- spectra[order(spectra[, 1]), ];    
  cnt     <- nrow(spectra);
  #isomatrix <- matrix(NA, ncol=5, nrow=0)
  #colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
  
  #calculate error
  error.ppm <- params$devppm * mz;
  #error.abs <- ),1, function(x) x + params$mzabs*rbind(1,2,3)));
  
  #for every peak in pseudospectrum
  for ( j in 1:(length(mz) - 1)){
    #create distance matrix
    MI <- spectra[j:cnt, 1] - spectra[j, 1];
    #Sum up all possible/allowed isotope distances + error(ppm of peak mz and mzabs)
    max.index <- max(which(MI < (sum(params$IM[1:params$maxiso, "mzmax"]) + error.ppm[j] + params$mzabs )))
    #check if one peaks falls into isotope window
    if(max.index == 1){
      #no promising candidate found, move on
      next;
    }
    
    #IM - isotope matrix (column diffs(min,max) per charge, row num. isotope)
    IM <- t(sapply(1:params$maxcharge,function(x){
      mzmin <- (params$IM[, "mzmin"] - error.ppm[j]-params$mzabs*x) / x
      mzmax <- (params$IM[, "mzmax"] + error.ppm[j]+params$mzabs*x) / x
      res   <- c(0,0);
      for(k in 1:length(mzmin)){
        res <- c(res, mzmin[k]+res[2*k-1], mzmax[k]+res[2*k])
      }
      return (res[-c(1:2)])
    } ))
    
    #find peaks, which m/z value is in isotope interval
    hits <- t(apply(IM, 1, function(x){ findInterval(MI[1:max.index], x)}))
    rownames(hits) <- c(1:nrow(hits))
    colnames(hits) <- c(1:ncol(hits))
    
    
    #checking first isotopic peak
    hit <- apply(hits, 1, function(x) any(x==1))
    hits <- hits[hit, , drop=FALSE]
    
    if(nrow(hits) == 0){
      next;
    }
    
    #getting max. isotope cluster length
    #TODO: unique or not????
    #isolength <- apply(hits, 1, function(x) length(which(unique(x) %% 2 !=0)))
    isohits   <- lapply(1:nrow(hits), function(x) which(hits[x, ] %% 2 !=0))
    isolength <- sapply(isohits, length)

    #Check if any result is found
    if(all(isolength==0)){
      next;
    }
    
    #itensity checks
    #candidate.matrix
    #first column - how often succeded the isotope intensity test
    #second column - how often could a isotope int test be performed
    candidate.matrix <- matrix(0, nrow=length(isohits), ncol=max(isolength)*2);
    
    for(iso in 1:length(isohits)){
      for(candidate in 1:length(isohits[[iso]])){
        for(sample.index in c(1:ncol(int))){
          #Test if C12 Peak is NA
          if(!is.na(int[j, sample.index])){              
            #candidate.matrix[maxIso, 1] <- candidate.matrix[maxIso, 1] + 1
          }
          charge <- as.numeric(row.names(hits)[iso])
          int.c12 <- int[j, sample.index]
          isotopePeak <- hits[iso,isohits[[iso]][candidate]]%/%2 + 1;
          if(isotopePeak == 1){
            #first isotopic peak, check C13 rule
            int.c13 <- int[isohits[[iso]][candidate]+j-1, sample.index];
            int.available <- all(!is.na(c(int.c12, int.c13)))
            if (int.available){
              theo.mass <- spectra[j, 1] * charge; #theoretical mass
              numC      <- round(theo.mass / 12); #max. number of C in molecule
              inten.max <- int.c12 * numC * 0.011; #highest possible intensity
              inten.min <- int.c12 * 1    * 0.011; #lowest possible intensity
              if((int.c13 < inten.max && int.c13 > inten.min) || !params$filter){
                candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }else{
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }
            } else {
              #todo
            } 
          } else {
            #x isotopic peak
            int.cx <- int[isohits[[iso]][candidate]+j-1, sample.index];
            int.available <- all(!is.na(c(int.c12, int.cx)))
            if (int.available) {
              intrange <- c((int.c12 * params$IM[isotopePeak,"intmin"]/100),
                            (int.c12 * params$IM[isotopePeak,"intmax"]/100))
              #filter Cx isotopic peaks muss be smaller than c12
              if(int.cx < intrange[2] && int.cx > intrange[1]){
                candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1                        
              }else{
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }
            } else {
              candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
            }#end int.available
          }#end if first isotopic peak
        }#for loop samples
      }#for loop candidate
    }#for loop isohits
    
    #calculate ratios
    candidate.ratio <- candidate.matrix[, seq(from=1, to=ncol(candidate.matrix),
                       by=2)] / candidate.matrix[, seq(from=2, 
                      to=ncol(candidate.matrix), by=2)];
    if(is.null(dim(candidate.ratio))){
      candidate.ratio <- matrix(candidate.ratio, nrow=1)
    }
    if(any(is.nan(candidate.ratio))){
      candidate.ratio[which(is.nan(candidate.ratio))] <- 0;
    }
    
    #decision between multiple charges or peaks
    for(charge in 1:nrow(candidate.matrix)){
      if(any(duplicated(hits[charge, isohits[[charge]]]))){
        #One isotope peaks has more than one candidate
        ##check if problem is still consistent
        for(iso in unique(hits[charge, isohits[[charge]]])){
          if(length(index <- which(hits[charge, isohits[[charge]]]==iso))== 1){
            #now duplicates next
            next;
          }else{
            #find best
            index2 <- which.max(candidate.ratio[charge, index]);
            save.ratio <- candidate.ratio[charge, index[index2]]
            candidate.ratio[charge,index] <- 0
            candidate.ratio[charge,iso] <- save.ratio
            index <- index[-index2]
            isohits[[charge]] <- isohits[[charge]][-index]
          }
        }
      }
      
      for(isotope in 1:ncol(candidate.ratio)){
        if(candidate.ratio[charge, isotope] >= params$minfrac){
          isomatrix <- rbind(isomatrix, 
                             c(spectra[j, 2],
                               spectra[isohits[[charge]][isotope]+j-1, 2], 
                               isotope, as.numeric(row.names(hits)[charge]), 0))
        } else{
          break;
        }
      }
    }#end for charge
  }#end for j
      
  return(isomatrix)
}

getderivativeIons <- function(annoID, annoGrp, rules, npeaks){
  #generate Vector length npeaks
  derivativeIons <- vector("list", npeaks);
  #intrinsic charge
  #TODO: Not working at the moment
  charge <- 0;
  
  #check if we have annotations
  if(nrow(annoID) < 1){
    return(derivativeIons);
  }
  
  for(i in 1:nrow(annoID)){
    
    peakid  <-  annoID[i, 1];
    grpid   <-  annoID[i, 2];
    ruleid  <-  annoID[i, 3];
    
    if(is.null(derivativeIons[[peakid]])){
      #Peak has no annotations so far
      if(charge == 0 | rules[ruleid, "charge"] == charge){
        derivativeIons[[peakid]][[1]] <- list( rule_id = ruleid, 
                                           charge = rules[ruleid, "charge"], 
                                           nmol = rules[ruleid, "nmol"], 
                                           name = paste(rules[ruleid, "name"]),
                                           mass = annoGrp[grpid, 2])
      }
    } else {
      #Peak has already an annotation
      if(charge == 0 | rules[ruleid, "charge"] == charge){
        derivativeIons[[peakid]][[(length(
          derivativeIons[[peakid]])+1)]] <- list( rule_id = ruleid, 
                                              charge = rules[ruleid, "charge"],
                                              nmol = rules[ruleid, "nmol"],
                                              name=paste(rules[ruleid, "name"]),
                                              mass=annoGrp[grpid, 2])
      }
    }
    
    charge <- 0;
  }
  return(derivativeIons);
}

getIsotopeCluster <- function(object, number=NULL, value="maxo", 
                              sampleIndex=NULL){
 
  #check values
  if(is.null(object)) { 
    stop("No xsa argument was given.\n"); 
  }else if(!class(object)=="xsAnnotate"){
    stop("Object parameter is no xsAnnotate object.\n");
  }
  
  value <- match.arg(value, c("maxo", "into", "intb"), several.ok=FALSE)

  if(!is.null(number) & !is.numeric(number)){
    stop("Number must be NULL or numeric");
  }

  if(!is.null(sampleIndex) & !all(is.numeric(sampleIndex))){
    stop("Parameter sampleIndex must be NULL or numeric");
  }
  
  if(is.null(sampleIndex)){
      nSamples <- 1;
  } else if( all(sampleIndex <= length(object@xcmsSet@filepaths) & sampleIndex > 0)){
      nSamples <- length(sampleIndex);
  } else {
      stop("All values in parameter sampleIndex must be lower equal 
         the number of samples and greater than 0.\n")
  }
  
  if(length(sampnames(object@xcmsSet)) > 1){  ## more than one sample
      gvals <- groupval(object@xcmsSet, value=value);
      groupmat <- object@groupInfo;
      iso.matrix <- matrix(0, ncol=nSamples, nrow=length(object@isotopes));
      if(is.null(sampleIndex)){
        for(i in 1:length(object@pspectra)){
          iso.matrix[object@pspectra[[i]],1] <- gvals[object@pspectra[[i]],object@psSamples[i]]; 
        }
      } else {
        for(i in 1:length(object@pspectra)){
          iso.matrix[object@pspectra[[i]], ] <- gvals[object@pspectra[[i]], sampleIndex]
        }
      }
      peakmat <- cbind(groupmat[, "mz"], iso.matrix );
      rownames(peakmat) <- NULL;
      if(is.null(sampleIndex)){
        colnames(peakmat) <- c("mz",value);
      }else{
        colnames(peakmat) <- c("mz", sampnames(object@xcmsSet)[sampleIndex]);
      }
      
      if(any(is.na(peakmat))){
        cat("Warning: peak table contains NA values. To remove apply fillpeaks on xcmsSet.\n");
      }
      
   } else if(length(sampnames(object@xcmsSet)) == 1){  ## only one sample was 
      peakmat <- object@groupInfo[, c("mz", value)];
   } else { 
     stop("sampnames could not extracted from the xcmsSet.\n"); 
   }

  #collect isotopes

  index <- which(!sapply(object@isotopes, is.null));

  tmp.Matrix <- cbind(index, matrix(unlist(object@isotopes[index]), ncol=4, byrow=TRUE))
  colnames(tmp.Matrix) <- c("Index","IsoCluster","Type","Charge","Val")

  max.cluster <- max(tmp.Matrix[,"IsoCluster"])
  max.type    <- max(tmp.Matrix[,"Type"])
  
  isotope.Matrix <- matrix(NA, nrow=max.cluster, ncol=(max.type+2));
  invisible(apply(tmp.Matrix,1, function(x) {
    isotope.Matrix[x["IsoCluster"],x["Type"]+2] <<- x["Index"];
    isotope.Matrix[x["IsoCluster"],1] <<- x["Charge"];
   }))

  invisible(apply(isotope.Matrix,1, function(x) {
    list(peaks=peakmat[na.omit(x[-1]),],charge=x[1])
  }))
}
