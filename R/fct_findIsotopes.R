##Functions for findIsotopes
calcIsotopeMatrix <- function(maxiso,maxcharge){
        M_C12_C13 = 1.0033

        ## Calculate ISO/CHARGE - Matrix
        IM <- matrix(NA,maxiso,maxcharge)
        for (i in 1:maxiso)
        for (j in 1:maxcharge) IM[i,j] <- (i* M_C12_C13) / j

        return(IM)

}

getderivativeIons <- function(annoID,annoGrp,rules,npeaks){
    derivativeIons<-vector("list",npeaks);
    charge=0;
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
    stop("No xsa argument was given"); 
  }else if(!class(object)=="xsAnnotate"){
    stop("xsa is no xsAnnotate object");
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