##Functions for findIsotopes
calcIsotopeMatrix <- function(maxiso,maxcharge){
        M_C12_C13 <- 1.0033

        ## Calculate ISO/CHARGE - Matrix
        IM <- matrix(NA,maxiso,maxcharge)
        for (i in 1:maxiso)
        for (j in 1:maxcharge) IM[i,j] <- (i* M_C12_C13) / j

        return(IM)

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
      
  #for every peak in pseudospectrum
  for ( j in 1:(length(mz) - 1)){
    #create distance matrix
    MI <- spectra[j:cnt, 1] - spectra[j, 1];
    max.index <- max(which(MI < params$IM[params$maxiso, 1] + max(2*params$devppm*mz) + params$mzabs))
    #for every charge
    for(charge in params$maxcharge:1){
      m <- fastMatch(params$IM[, charge], MI[1:max.index], tol = max(2*params$devppm*mz) + params$mzabs)
      #for every match, test of isotopes existing
      if(any(!sapply(m, is.null))){
        #checking first isotopic peak
        if(is.null(m[[1]])){
          next;
        }    
        isolength <- sapply(m, length)
        maxIso    <- which(isolength == 0)[1];
        if(is.na(maxIso)){
          #all Isotope peaks have been found
          maxIso <- params$maxiso + 1;  
        }
        #candidate.matrix
        #first column - how often succeded the isotope intensity test
        #second column - how often could a isotope int test be performed
        candidate.matrix <- matrix(0, nrow=maxIso, ncol=max(isolength)*2);
        #for every sample
        for(sample.index in c(1:ncol(int))){
          if(!is.na(int[j, sample.index])){              
            candidate.matrix[maxIso,1 ] <- candidate.matrix[maxIso,1] + 1
          }
          #for every isotope mz match
          for( iso in 1:(maxIso-1)){
            for( candidate in 1:(ncol(candidate.matrix)/2)){
              if (iso == 1){
                # first isotopic peak
                # check C13 rule
                int.c13 <- int[m[[iso]][candidate]+j-1, sample.index];
                int.c12 <- int[j, sample.index]
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
                int.cx <- int[m[[iso]][candidate]+j-1, sample.index];
                int.c12  <- int[j, sample.index];
                int.available <- all(!is.na(c(int.c12, int.cx)))
                if (int.available) {
                  #filter Cx isotopic peaks muss be smaller than c12
                  if(int.cx < int.c12){
                    candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                    candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1                        
                  }else{
                    candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                  }
                } else {
                  candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
                }#end int.available
              }#end if iso
            }#end for candidate
          }#end for iso
        }#end for sample.index
      } else {#end if is.null
        next;
      }

      ##TODO: Implement intrinsic charges
      #evaluate if candiates are real isotopic peaks
      #First C13 - find best candidate - best candidate = best hit / occurence rate
      index <- which.max(res <- candidate.matrix[1, 
                          seq(1, ncol(candidate.matrix), by=2)] / 
                            candidate.matrix[1, seq(1, ncol(candidate.matrix), by=2) + 1])
      if (length(index) > 0 && res[index] >= params$minfrac){
        #C13 Peak over minfrac threshold
        isomatrix <- rbind(isomatrix, c(spectra[j, 2], spectra[m[[1]][index]+j-1, 2], 1, charge, 0))
        maxhits <- candidate.matrix[1,index*2]
        #Check if we have more than one isotopic candidate
        if(nrow(candidate.matrix) > 2){
          for(z in 2:(nrow(candidate.matrix)-1)){
            #Check z'th isotopic peak
            #Select best candidate - here most hits
            for( zi in seq(1,ncol(candidate.matrix),by=2)){
              if(spectra[m[[z]][zi]+j-1,2] %in% isomatrix[,2]){
                candidate.matrix[z,zi] <- -1;
              }
            }
            index <- which.max(res <- candidate.matrix[z,seq(1,ncol(candidate.matrix),by=2)])
            if(res == -1){
              #Best isotope peak was already used. No more isotope peaks can be applied.
              break;
            }
            if(candidate.matrix[z,index*2] <= candidate.matrix[maxIso,1] & all(candidate.matrix[z,index*2-1] <= maxhits)){
              maxhits <- c(maxhits,candidate.matrix[z,index*2-1]);
              isomatrix<-rbind(isomatrix,c(spectra[j,2],spectra[m[[z]][index]+j-1,2],z,charge,0))
            }else{
              break;
            }
          }
        }
      }   
    }#end for charge
  }#end for j
  return(isomatrix)
}

getderivativeIons <- function(annoID,annoGrp,rules,npeaks){
    derivativeIons <- vector("list", npeaks);
    charge <- 0;
    if(nrow(annoID) < 1){
      return(derivativeIons);
    }
    for(i in 1:nrow(annoID)){
        peakid  <-  annoID[i,1];
        grpid   <-  annoID[i,2];
        ruleid  <-  annoID[i,3];
        if(is.null(derivativeIons[[peakid]])){
            if(charge==0 | rules[ruleid,"charge"]==charge){
                derivativeIons[[peakid]][[1]]<- list( rule_id = ruleid, charge=rules[ruleid,"charge"], nmol= rules[ruleid,"nmol"], name=paste(rules[ruleid,"name"]),mass=annoGrp[grpid,2])
            }
        }else{
            if(charge==0 | rules[ruleid,"charge"]==charge){
                derivativeIons[[peakid]][[(length(derivativeIons[[peakid]])+1)]] <- list( rule_id = ruleid, charge=rules[ruleid,"charge"], nmol= rules[ruleid,"nmol"], name=paste(rules[ruleid,"name"]),mass=annoGrp[grpid,2])
            }
        }
        charge=0;
    }
return(derivativeIons);
}

getIsotopeCluster <- function(object, number=NULL, value="maxo"){
 
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

  if(length(sampnames(object@xcmsSet)) > 1){  ## more than one sample
      gvals <- groupval(object@xcmsSet, value=value);
      groupmat <- object@groupInfo;
      iso.vector <- vector(mode="numeric",length=length(object@isotopes));
      for(i in 1:length(object@pspectra)){
         iso.vector[object@pspectra[[i]]] <- gvals[object@pspectra[[i]],object@psSamples[i]]; 
      }
      peakmat <- cbind(groupmat[,"mz"],iso.vector);
      rownames(peakmat) <- NULL;
      colnames(peakmat) <- c("mz",value);
      if(any(is.na(peakmat[,value]))){
        cat("Warning: peak table contains NA values. To remove apply fillpeaks on xcmsSet.\n");
      }
   } else if(length(sampnames(object@xcmsSet)) == 1){  ## only one sample was 
      peakmat <- object@groupInfo[,c("mz",value)];
   }else { stop("sampnames could not extracted from the xcmsSet.\n"); }

  #collect isotopes

  index <- which(!sapply(object@isotopes,is.null));

  tmp.Matrix <- cbind(index,matrix(unlist(object@isotopes[index]),ncol=4,byrow=TRUE))
  colnames(tmp.Matrix) <- c("Index","IsoCluster","Type","Charge","Val")

  max.cluster <- max(tmp.Matrix[,"IsoCluster"])
  max.type    <- max(tmp.Matrix[,"Type"])
  
  isotope.Matrix <- matrix(NA,nrow=max.cluster,ncol=(max.type+2));
  invisible(apply(tmp.Matrix,1, function(x) {
    isotope.Matrix[x["IsoCluster"],x["Type"]+2] <<- x["Index"];
    isotope.Matrix[x["IsoCluster"],1] <<- x["Charge"];
   }))

  invisible(apply(isotope.Matrix,1, function(x) {
    list(peaks=peakmat[na.omit(x[-1]),],charge=x[1])
  }))
}