extractfragments <- function(object, pspecidx=NULL) {

  if(class(object) != "xsAnnotate"){
    stop("Object parameter is no xsAnnotate object.\n");
  }
  
  if (is.null(pspecidx)) {
    pspecidx <- seq(along=object@pspectra)
  }
  
  result <- list();
  ruleLoss <- which(object@ruleset$massdiff > 0);
  for (pspec in pspecidx) {

    sp <- getpspectra(object, grp=pspec)
    rt <- median(sp[, "rt"])
    peakID <- object@pspectra[[pspec]]
    if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
        #single sample
        #intIdx <- intVal;
        intIdx <- "maxo";
    }else{
        ##multiple sample
        intIdx <- colnames(object@groupInfo)[object@psSamples[pspec]+7+length(levels(sampclass(object@xcmsSet)))]
    }
    
    ## All peaks without annotation
    nonannoidx <- grepl("^$", sp[,"adduct"])
    
    ## Annotation alternatives 
    annos <- object@annoGrp[object@annoGrp[, "psgrp"]==pspec, , drop=FALSE]
    ipsidx <- order(annos[, "ips"], decreasing=TRUE)
    annos <- annos[ipsidx, , drop=FALSE]
    
    ## only export pspec for which at least some annos exist
    ## maybe be more specific later ?
    ##
    for (i in seq(length.out=nrow(annos))) {
      annoid <- annos[i, "id"]
      neutralmass <- object@annoGrp[annoid, "mass"]
      

      ## All peaks in sp which *could* be a fragment of neutralmass
      fragpeakidx <- sp[, "mz"] < neutralmass

      ## All peaks in sp annotated with *this* annoid
      annopeakidx <- object@pspectra[[pspec]] %in% object@annoID[
        which(object@annoID[, "grpID"]==annoid & (object@annoID[,"ruleID"] %in%(ruleLoss))), "id"]
      fragmentidx <- object@pspectra[[pspec]] %in% object@annoID[
	which(object@annoID[, "grpID"]==annoid & !(object@annoID[,"ruleID"] %in%(ruleLoss))), "id"]
	annopeakidx <- annopeakidx & !fragmentidx
      ## unannotated OR monoisotopic peaks
      nonisoidx <- sapply(peakID, function(x) {
                            if(x %in% object@isoID[, "isopeak"]) {
                              return(FALSE) 
                            } else {
                              return(TRUE)
                            }})
    
      #nonisoidx <- grepl("(^$)|( 0 +)", sp[, "isotopes"])

      ## all known-to-be-multiply-charges
      multichargeidx <- object@pspectra[[pspec]] %in% sapply(object@annoID[object@annoID[, "grpID"]==annoid, "id"], 
                                                         function(x) {
                                                           ion <- object@derivativeIons[[x]];
                                                           if (length(ion) > 0) { 
                                                             if (ion[[1]]$charge > 1) {
                                                                return(unlist(x)) 
                                                             } else {
                                                                return (NA) 
                                                             } 
                                                           } else {
                                                               return(NA)
                                                           }}) & annopeakidx

      filteredsp <- sp[!multichargeidx & fragpeakidx & nonisoidx & ( nonannoidx |  !annopeakidx ), , drop=FALSE]
      
      o <- order(filteredsp[, "mz"])
      if (length(o) > 0) {
        result[[length(result)+1]] <- list(Pseudospectrum=pspec, AnnotationID=annoid, RentionTime=rt, 
                                           ParentMass=neutralmass, Peaks=filteredsp[o, c("mz", intIdx), drop=FALSE])
      }
    }
  }
  invisible(result)
}  


fragments2metfrag <- function(result, filedir=NULL) {
  if(is.null(filedir)){
    return(invisible(result))
  }else{
    lapply(result, function(x){
      zz <- file(paste(filedir, x$Pseudospectrum, "_", x$AnnotationID, ".mf", sep=""), "w");
      cat("# Mode: 1\n", file=zz)
      cat("# Pseudospectrum:", x$Pseudospectrum, "\n", file=zz)
      cat("# Annotation alternative:", x$AnnotationID, "\n", file=zz)
      cat("# Retentiontime:", x$RentionTime, "\n", file=zz)
      cat("# Parent Mass: ", x$ParentMass, "\n", file=zz)
      write.table(x$Peaks,row.names=FALSE, col.names=FALSE, file=zz)
      close(zz)
      })
    return(NULL);
  }
}

pspec2metfrag <- function(object, pspecidx=NULL, filedir=NULL) {
  result <- extractfragments(object, pspecidx) 
  fragments2metfrag(result, filedir)
}  

## CAMERA pseudo-spectra to MetFusion

pspec2metfusion <- function(object, pspecidx=NULL, filedir=NULL) {
  result <- extractfragments(object, pspecidx) 
  fragments2metfusion(result, filedir)
}

fragments2metfusion <- function(result, filedir=NULL) {
  if(is.null(filedir)){
    return(invisible(result))
  }else{
    lapply(result, function(x){
      zz <- file(paste(filedir, x$Pseudospectrum, "_", x$AnnotationID, ".mf", sep=""), "w");
      cat("# mfDatabaseIDs: \n", file=zz)
      cat("# mfLimit: 10000", "\n", file=zz)
      cat("# mfFormula: ", "\n", file=zz)
      cat("# mbLimit: 100", "\n", file=zz)
      cat("# mfDatabase: ", "chemspider", "\n", file=zz)    # kegg, pubchem, chemspider, beilstein, knapsack, massbank, sdf
      cat("# mfMZabs: 0.01", "\n", file=zz)
      cat("# mfMZppm: 30.0", "\n", file=zz)
      cat("# clustering: true", "\n", file=zz)
      cat("# mbInstruments: CE-ESI-TOF,ESI-IT-MS/MS,ESI-QqIT-MS/MS,ESI-QqQ-MS/MS,ESI-QqTOF-MS/MS,LC-ESI-IT,LC-ESI-ITFTLC-ESI-ITTOF,LC-ESI-Q,LC-ESI-QIT,LC-ESI-QQ,LC-ESI-QTOF", "\n", file=zz)
      cat("# mbCutoff: 5", "\n", file=zz)
      cat("# mfAdduct: ", "Neutral", "\n", file=zz)
      cat("# mfSearchPPM: 30.0", "\n", file=zz)
      cat("# mfParentIon: ", x$ParentMass, "\n", file=zz)
      cat("# mbIonization: ", substr(object@polarity,1,3), "\n", file=zz)  # pos, neg, both
      cat("# mfExactMass: ", x$ParentMass, "\n", file=zz)
      #cat("# Mode: 1\n", file=zz)
      #cat("# Pseudospectrum:", x$Pseudospectrum, "\n", file=zz)
      #cat("# Annotation alternative:", x$AnnotationID, "\n", file=zz)
      #cat("# Retentiontime:", x$RentionTime, "\n", file=zz)
      #cat("# Parent Mass: ", x$ParentMass, "\n", file=zz)
      cat("# peaks: ", "\n", file=zz)
      write.table(x$Peaks,row.names=FALSE, col.names=FALSE, file=zz)
      close(zz)
    })
    return(NULL);
  }
}  
