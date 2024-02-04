setGeneric("plotEICs", function(object,
                                pspec=1:length(object@pspectra),
                                maxlabel=0, sleep=0,
                                ...) standardGeneric("plotEICs"))

setMethod("plotEICs", "xsAnnotate", function(object,
                                             pspec=1:length(object@pspectra),
                                             maxlabel=0, sleep=0, method="bin"){
  #Which different samples need to be addressed for extracting raw data
  smpls <- unique(object@psSamples[pspec])
           
  xeic  <- new("xcmsEIC");
  xeic@eic <- vector("list", length(pspec))
  xeic@rtrange <- matrix(nrow=length(pspec), ncol=2)
  xeic@mzrange <- matrix(nrow=length(pspec), ncol=2)
  #iterator for ps-grp
  cnt <- 0;
  #one second overlap
  rtmargin <- 1; 
  
  for (a in seq(along=smpls)) { ## sample-wise EIC collection
    #read rawData into one xcmsRaw
    xraw <- xcmsRaw(object@xcmsSet@filepaths[smpls[a]], profmethod=method)
    xraw@scantime <- object@xcmsSet@rt$corrected[[smpls[a]]]

    pspecS <- pspec[which(object@psSamples[pspec] == smpls[a])]
    ## getting ALL peaks from the current sample (not that bad)
    peaks <- CAMERA:::getPeaks(object@xcmsSet, smpls[a])  
    invisible(lapply (pspecS, function(pc) {
      cnt <<- cnt + 1
      pidx <- object@pspectra[[pc]]
      pks <- peaks[pidx, , drop=FALSE]
      gks <- object@groupInfo[pidx, , drop=FALSE]
      
      nap <- which(is.na(pks[, 1]))
      pks[nap, ] <- cbind(gks[nap, c(1:6), drop=FALSE], matrix(nrow=length(nap), ncol=5, 0))
      bbox <- c(rtmin = min(pks[, "rtmin"])-rtmargin,
                rtmax = max(pks[, "rtmax"])+rtmargin,
                mzmin = min(pks[, "mzmin"]),
                mzmax = max(pks[, "mzmax"]))
      eic <- xcms:::getEIC(xraw, rtrange=pks[, c("rtmin", "rtmax"), drop=FALSE],
                                 mzrange=pks[, c("mzmin", "mzmax"), drop=FALSE])
      #write resulting bounding box into xcmsEIC
      xeic@rtrange[pc, ] <<- bbox[c("rtmin","rtmax")]
      xeic@mzrange[pc, ] <<- bbox[c("mzmin","mzmax")]
      cat("--> ", pc, " (", cnt, " of ", length(pspec), ")", "\n", sep = "")
      xeic@eic[[pc]] <<- eic@eic[[1]]
    }))
  }
  names(xeic@eic) <- paste("Pseudospectrum ", pspec, sep="")

 if(maxlabel > 0){
  col <- rainbow(maxlabel);
 } else {
  col <- c();
 }
  
 ##
 ## Loop through all pspectra
 ##

  for (ps in seq(along=pspec)) {
    EIC <- xeic@eic[[ps]];
    pidx <- object@pspectra[[pspec[ps]]];
    peaks <- CAMERA:::getPeaks(object@xcmsSet,object@psSamples[pspec[ps]])[pidx,,drop=FALSE]
    grps <- object@groupInfo[pidx, ]
    nap <- which(is.na(peaks[, 1]))
    naps <- rep(FALSE, nrow(peaks))
    
    if(length(nap) > 0){  
      naps[nap] <- TRUE;
      peaks[nap,] <- cbind(grps[nap,c(1:6), drop=FALSE], matrix(nrow=length(nap), ncol=5,0))
    }
    
    main <- paste("Pseudospectrum ", pspec[ps], sep="")
    
    ## Calculate EICs and plot ranges
    neics     <- length(pidx)
    lmaxlabel <- min(maxlabel, neics)
    eicidx <- 1:neics
    maxint <- numeric(length(eicidx))

    for (j in eicidx) {
        maxint[j] <- max(EIC[[j]][, "intensity"])
    }
    o  <- order(maxint, decreasing = TRUE)
    rt <- xeic@rtrange[ps, ];    
    rt.min <- round(mean(peaks[, "rtmin"]), digits=3);
    rt.med <- round(peaks[o[1], "rt"], digits=3);
    rt.max <- round(mean(peaks[, "rtmax"]), digits=3);
    ## Open Plot

    plot(0, 0, type = "n", xlim = rt, ylim = c(0, max(maxint)),xaxs='i',
                   xlab = "Retention Time (seconds)", ylab = "Intensity",
                   main = paste("Extracted Ion Chromatograms for ", main,"\nTime: From",rt.min,"to",rt.max,", mean",rt.med))
  
    ## Plot Peak and surrounding
    lcol <- rgb(0.6, 0.6, 0.6);
    lcol <- c(col, rep(lcol, max(nrow(peaks) - maxlabel, 0)));
    cnt  <- 1;
    for (j in eicidx[o]) {
      pts <- xeic@eic[[ps]][[j]]
      points(pts, type = "l", col = lcol[cnt]);
      cnt <- cnt + 1;
      peakrange <- peaks[,c("rtmin","rtmax"), drop=FALSE]
      ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
      if (naps[j]){ 
        points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3, lty=3)
      } else {
        points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3)
      }
    }
    ## Plot Annotation Legend
    pspectrum <- getpspectra(object, grp=pspec[ps])
    mz <- pspectrum[o, "mz"];
    #Check adduct annotation
    if (lmaxlabel > 0 & "adduct" %in% colnames(pspectrum)) {
        adduct <- sub("^ ", "", pspectrum[o, "adduct"]) #Remove Fronting Whitespaces
        mass   <- sapply(strsplit(adduct, " "), function(x) {x[2]})
        adduct <- sapply(strsplit(adduct, " "), function(x) {x[1]})
        umass <- unique(na.omit(mass[1:maxlabel]));
        adduct[is.na(adduct)] <- "";
        test <- vector("list",length=length(mz));
        mz <- format(pspectrum[o[1:lmaxlabel], "mz"], digits=5);
        
        if(length(umass) > 0){
          for(i in 1:length(umass)){
            ini <- which(mass==umass[i]);
            for(ii in 1:length(ini)){
              firstpart  <- strsplit(adduct[ini[ii]], "M")[[1]][1]
              secondpart <- strsplit(adduct[ini[ii]], "M")[[1]][2]
              masspart   <- mz[ini[ii]];
              test[[ini[ii]]] <- substitute(paste(masspart, " ", firstpart, M[i], secondpart), 
                                            list(firstpart=firstpart, i=i, secondpart=secondpart,
                                                 masspart = masspart))
            }
          }
        }
        
         for(i in 1:lmaxlabel){      
          if(is.null(test[[i]])){
            test[[i]] <- mz[i];
          }
        }
        
        leg <- as.expression(test[1:lmaxlabel]);  
        legend("topright", legend=leg, col=lcol, lty=1)
    }
            
    if (sleep > 0) {
      Sys.sleep(sleep)
    }
  }
})

setGeneric("plotPeakEICs", function(object,
                                    pidx=1:nrow(object@groupInfo),
                                samples=min(object@psSamples):max(object@psSamples),
                                maxlabel=0, sleep=0,
                                ...) standardGeneric("plotPeakEICs"))

setMethod("plotPeakEICs", "xsAnnotate", function(object, pidx=1:nrow(object@groupInfo),
                                             samples=min(object@psSamples):max(object@psSamples),
                                             maxlabel=0, sleep=0, method="bin"){
  #Which different samples need to be addressed for extracting raw data
  if (!all(samples %in% object@psSamples)){
    stop ("Parameter samples must be within the number of samples.\n");
  }
  if (!all(pidx <= nrow(object@groupInfo && pidx > 0))){
    stop ("Peadindexes must be greater 0 and lower as maximum number of peaks.\n")
  }
  xeic  <- new("xcmsEIC");
  xeic@rtrange <- matrix(nrow=nrow(object@groupInfo), ncol=2)
  xeic@mzrange <- matrix(nrow=nrow(object@groupInfo), ncol=2)
  #iterator for ps-grp
  pcpos <- 1;
  #one second overlap
  rtmargin <- 1; 
  
  for (a in seq(along=samples)) { ## sample-wise EIC collection
    #read rawData into one xcmsRaw
    xraw <- xcmsRaw(object@xcmsSet@filepaths[samples[a]], profmethod=method)
    
    ## getting ALL peaks from the current sample (not that bad)
    peaks <- CAMERA:::getPeaks(object@xcmsSet, a)  
    eic   <- sapply (peaks, function(pks) {
      pidx <- object@pspectra[[pc]]
      pks <- peaks[pidx, , drop=FALSE]
      gks <- object@groupInfo[pidx, , drop=FALSE]
      
      nap <- which(is.na(pks[, 1]))
      pks[nap, ] <- cbind(matrix(nrow=length(nap), ncol=ncol(pks), 0))
      bbox <- c(rtmin = min(pks[, "rtmin"])-rtmargin,
                rtmax = max(pks[, "rtmax"])+rtmargin,
                mzmin = min(pks[, "mzmin"]),
                mzmax = max(pks[, "mzmax"]))
      eic <- xcms:::getEIC(xraw, rtrange=pks[, c("rtmin", "rtmax"), drop=FALSE],
                           mzrange=pks[, c("mzmin", "mzmax"), drop=FALSE])
      #write resulting bounding box into xcmsEIC
      xeic@rtrange[pcpos, ] <<- bbox[c("rtmin","rtmax")]
      xeic@mzrange[pcpos, ] <<- bbox[c("mzmin","mzmax")]  
      eic@eic[[1]]
    })
    xeic@eic <- c(xeic@eic, eic)
  }
  names(xeic@eic) <- paste("Pseudospectrum ", pspec, sep="")
  
  if(maxlabel > 0){
    col <- rainbow(maxlabel);
  } else {
    col <- c();
  }
  
  ##
  ## Loop through all pspectra
  ##
  
  for (ps in seq(along=pspec)) {
    EIC <- xeic@eic[[ps]];
    pidx <- object@pspectra[[pspec[ps]]];
    peaks <- CAMERA:::getPeaks(object@xcmsSet,object@psSamples[pspec[ps]])[pidx,,drop=FALSE]
    grps <- object@groupInfo[pidx, ]
    nap <- which(is.na(peaks[, 1]))
    naps <- rep(FALSE, nrow(peaks))
    
    if(length(nap) > 0){  
      naps[nap] <- TRUE;
      peaks[nap,] <- cbind(grps[nap,c(1:6), drop=FALSE], matrix(nrow=length(nap), ncol=5,0))
    }
    
    main <- paste("Pseudospectrum ", pspec[ps], sep="")
    
    ## Calculate EICs and plot ranges
    neics     <- length(pidx)
    lmaxlabel <- min(maxlabel, neics)
    eicidx <- 1:neics
    maxint <- numeric(length(eicidx))
    
    for (j in eicidx) {
      maxint[j] <- max(EIC[[j]][, "intensity"])
    }
    o  <- order(maxint, decreasing = TRUE)
    rt <- xeic@rtrange[ps, ];    
    rt.min <- round(mean(peaks[, "rtmin"]), digits=3);
    rt.med <- round(peaks[o[1], "rt"], digits=3);
    rt.max <- round(mean(peaks[, "rtmax"]), digits=3);
    ## Open Plot
    
    plot(0, 0, type = "n", xlim = rt, ylim = c(0, max(maxint)),xaxs='i',
         xlab = "Retention Time (seconds)", ylab = "Intensity",
         main = paste("Extracted Ion Chromatograms for ", main,"\nTime: From",rt.min,"to",rt.max,", mean",rt.med))
    
    ## Plot Peak and surrounding
    lcol <- rgb(0.6, 0.6, 0.6);
    lcol <- c(col, rep(lcol, max(nrow(peaks) - maxlabel, 0)));
    cnt  <- 1;
    for (j in eicidx[o]) {
      pts <- xeic@eic[[ps]][[j]]
      points(pts, type = "l", col = lcol[cnt]);
      cnt <- cnt + 1;
      peakrange <- peaks[,c("rtmin","rtmax"), drop=FALSE]
      ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
      if (naps[j]){ 
        points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3, lty=3)
      } else {
        points(pts[ptsidx, ], type = "l", col = col[j], lwd=1.3)
      }
    }
    ## Plot Annotation Legend
    pspectrum <- getpspectra(object, grp=pspec[ps])
    mz <- pspectrum[o, "mz"];
    #Check adduct annotation
    if (lmaxlabel > 0 & "adduct" %in% colnames(pspectrum)) {
      adduct <- sub("^ ", "", pspectrum[o, "adduct"]) #Remove Fronting Whitespaces
      mass   <- sapply(strsplit(adduct, " "), function(x) {x[2]})
      adduct <- sapply(strsplit(adduct, " "), function(x) {x[1]})
      umass <- unique(na.omit(mass[1:maxlabel]));
      adduct[is.na(adduct)] <- "";
      test <- vector("list",length=length(mz));
      mz <- format(pspectrum[o[1:lmaxlabel], "mz"], digits=5);
      
      if(length(umass) > 0){
        for(i in 1:length(umass)){
          ini <- which(mass==umass[i]);
          for(ii in 1:length(ini)){
            firstpart  <- strsplit(adduct[ini[ii]], "M")[[1]][1]
            secondpart <- strsplit(adduct[ini[ii]], "M")[[1]][2]
            masspart   <- mz[ini[ii]];
            test[[ini[ii]]] <- substitute(paste(masspart, " ", firstpart, M[i], secondpart), 
                                          list(firstpart=firstpart, i=i, secondpart=secondpart,
                                               masspart = masspart))
          }
        }
      }
      
      for(i in seq(along=lmaxlabel)){
        if(is.null(test[[i]])){
          test[[i]] <- mz[i];
        }
      }
      
      leg <- as.expression(test[1:lmaxlabel]);  
      legend("topright", legend=leg, col=lcol, lty=1)
    }
    
    if (sleep > 0) {
      Sys.sleep(sleep)
    }
  }
})



setGeneric("plotPsSpectrum", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,mzrange = numeric(),
                                              sleep=0, cexMulti=1, ...) standardGeneric("plotPsSpectrum"))
setMethod("plotPsSpectrum", "xsAnnotate", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,mzrange = numeric(),
                                              sleep=0, cexMulti=1, ...){
##
## Loop through all requested pspectra
##
  if (is.na(object@sample) || length(object@sample)<1) {
    gvals <- groupval(object@xcmsSet);
    peakmat <- object@xcmsSet@peaks;
    groupmat <- groups(object@xcmsSet);
    #calculate highest peak
    max_mat  <- apply(gvals, 1, function(x, peakmat){peakmat[x, value]}, peakmat);
  }else{
    peakmat <- getPeaks(object@xcmsSet, index=object@sample);
    maxo  <- peakmat[,value]
  }
  for (psp in pspec) {
    pspectrum <- getpspectra(object, grp=psp);
    pindex<-object@pspectra[[psp]];
    if(is.na(object@sample)){
      intensity <- max_mat[object@psSamples[psp],pindex]
      intensity[which(is.na(intensity))] <- 1; #fix if NA
    }else{
      intensity <- maxo[pindex]
    }
    o <- order(intensity, decreasing=TRUE);
    mz <- pspectrum[o, "mz"]

    if (log) {
      intensity <- log(intensity[o])
    } else {
      intensity <- intensity[o]
    }

    if (length(mzrange) == 0) {
      mzrange <- round(range(mz));
    } else {
      mzrange <- c(min(mzrange), max(mzrange))
    }
    if (mzrange[1] == mzrange[2]){
      #give 10% to range
      mzrange[1] <- mzrange[1] - mzrange[1]*0.1
      mzrange[2] <- mzrange[2] + mzrange[2]*0.1
      masslab <- paste(mzrange[1], "-", mzrange[2], " m/z", sep="")
    }else{
      masslab <- paste(mzrange[1], "-", mzrange[2], " m/z, ", sep="")
    }
    mzborder <- (mzrange[2]-mzrange[1])*0.05
    rtrange <- paste(round(median(pspectrum[o,"rt"]),digits=2), "s")
    sample  <- paste(", Sample:",object@psSamples[psp]);
    if(is.null(title)){
      title.plot <- paste("Pseudospectrum",psp,rtrange,sample);
    }else{
      title.plot <- title;
    }
    index <- which(mz >= mzrange[1]-mzborder & mz <= mzrange[2]+mzborder);
    if(length(index) == 0){
        intrange <- range(intensity);
        intrange[1] <- intrange[1]*0.9
        intrange[2] <- intrange[2]*1.05
    }else{
        intrange <- range(intensity[index]);
        intrange[1] <- intrange[1]*0.9
        intrange[2] <- intrange[2]*1.05
    }
    if (maxlabel>0) {
      ## Also check for labels:  https://stat.ethz.ch/pipermail/r-help/2008-November/178666.html
      ## There is also the spread.labs function in the TeachingDemos package that
      ## uses a different method from the plotrix function and should not move any
      ## labels that are not overlapping.
#       index <- which(mz > mzrange[1]-mzborder & mz < mzrange[2]+mzborder);
      lmaxlabel <- min(maxlabel,length(index))
      if ("adduct" %in% colnames(pspectrum) && lmaxlabel > 0) {
        adduct <- sub("^ ","",pspectrum[o, "adduct"]) #Remove Fronting Whitespaces
        mass <- sapply(strsplit(adduct, " "), function(x) {x[2]})
        adduct <- sapply(strsplit(adduct, " "), function(x) {x[1]})
        umass <- unique(na.omit(mass[index[1:maxlabel]]));
        adduct[is.na(adduct)] <- "";
        test <- vector("list",length=length(mz));
        if(length(umass) > 0){
          legend.txt <- c();
          for(i in 1:length(umass)){
            ini <- which(mass==umass[i]);
            for(ii in 1:length(ini)){
              firstpart <- strsplit(adduct[ini[ii]], "M")[[1]][1]
              secondpart <- strsplit(adduct[ini[ii]], "M")[[1]][2]
              test[[ini[ii]]] <- substitute(paste(firstpart,M[i],secondpart),list(firstpart=firstpart,i=i,secondpart=secondpart))
            }
            legend.txt <- c(legend.txt,bquote(M[.(i)] == .(round(as.numeric(umass[i]),digits=5))))
          }
          leg <- as.expression(legend.txt);
          plot(mz, intensity, type="h",
            xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
            main=title.plot, col="darkgrey",ylim=intrange,...)        
          text(mz[index[1:lmaxlabel]],
               intensity[index[1:lmaxlabel]],
               labels=format(mz[index[1:lmaxlabel]], digits=5),
               cex=0.66*cexMulti)
          sapply(1:lmaxlabel, function (x) { text(mz[index[x]],
                    intensity[index[x]] + strheight("0123456789"),
                    labels=test[[index[x]]],cex=0.66*cexMulti)});
          legend("topleft",legend=leg,cex=0.8,bty="n");
        }else{
          plot(mz, intensity, type="h",
            xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
            main=title.plot, col="darkgrey",ylim=intrange,...)        
          text(mz[index[1:lmaxlabel]],
               intensity[index[1:lmaxlabel]],
               labels=format(mz[index[1:lmaxlabel]], digits=5),
               cex=0.66*cexMulti)
        }
      }else{
          plot(mz, intensity, type="h",
            xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
            main=title.plot, col="darkgrey",ylim=intrange,...)
          if(lmaxlabel > 0){
            text(mz[index[1:lmaxlabel]],
               intensity[index[1:lmaxlabel]],
               labels=format(mz[index[1:lmaxlabel]], digits=5),cex=0.66*cexMulti)
          }
      }
    }else{ 
      plot(mz, intensity, type="h",
        xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
        main=title.plot, col="darkgrey",ylim=intrange,...)        
    }
     
    if (sleep > 0) {
      Sys.sleep(sleep)
    }
  }
})
