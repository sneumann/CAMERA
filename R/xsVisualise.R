setGeneric("plotEICs", function(object,
                                pspecIdx=1:length(object@pspectra),
                                maxlabel=0, sleep=0,
                                ...) standardGeneric("plotEICs"))
setMethod("plotEICs", "xsAnnotate", function(object,
                                             pspecIdx=1:length(object@pspectra),
                                             maxlabel=0, sleep=0){
           #stop("angst oO")
           smpls <- unique(object@psSamples[pspecIdx])
           #Ranges <- matrix(ncol=4, nrow=length(pspecIdx))
           xeic <- new("xcmsEIC")
           xeic@rtrange <- matrix(nrow=length(pspecIdx), ncol=2)
           xeic@mzrange <- matrix(nrow=length(pspecIdx), ncol=2)
           pcpos <- 1
           for (a in 1:length(smpls)) { ## sample-wise EIC collection
               xraw <- xcmsRaw(object@xcmsSet@filepaths[smpls[a]])
               rtmargin <- 0.018 * max(xraw@scantime) 
               pspecS <- pspecIdx[which(object@psSamples[pspecIdx] == smpls[a])]
               peaks <- CAMERA:::getPeaks(object@xcmsSet,smpls[a])  ## getting ALL peaks from the current sample (not that bad)
               eic <- lapply (pspecS, function(pc) {
                   pidx <- object@pspectra[[pc]]
                   pks <- peaks[pidx,]
                   gks <- object@groupInfo[pidx,]
                   nap <- which(is.na(pks[,1]))
                   pks[nap,] <- cbind(gks[nap,c(1:6), drop=FALSE],matrix(nrow=length(nap),ncol=5,0))
                   bbox <- c(rtmin = min(pks[,"rtmin"])-rtmargin,
                             rtmax = max(pks[,"rtmax"])+rtmargin,
                             mzmin = min(pks[,"mzmin"]),
                             mzmax = max(pks[,"mzmax"]))
                   eic <- xcms:::getEIC(xraw,
                                 rtrange=matrix(bbox[c("rtmin","rtmax")],
                                 nrow=length(pidx), ncol=2, byrow=TRUE),
                                 mzrange=pks[,c("mzmin", "mzmax"), drop=FALSE])
                   xeic@rtrange[pcpos, ] <<- bbox[c("rtmin","rtmax")]
                   xeic@mzrange[pcpos, ] <<- bbox[c("mzmin","mzmax")]
                   cat("-->",pcpos,"\n")
                   pcpos <<- pcpos+1
                    ##stop("WTF?")
                   eic@eic[[1]]
                   })
               xeic@eic <- c(xeic@eic, eic)
            }
           names(xeic@eic) <- paste("Pseudospectrum ", pspecIdx, sep="")

  ## Setup Peak and Neighbourhood Colors,
  ## Now: all black
  ## Later: by Annotation 
  ##

 ## creating a lightcolor-vector:
  rgbcol <- matrix(nrow=25,ncol=3) ## generating 25 rgb-values
  colv <- c(0,0.6,1) ## color-intensities
  id<-c(1,1,1)
  for (a in 1:25){
      rgbcol[a,] <- c(colv[id[1]],colv[id[2]],colv[id[3]])
      if (id[3] < 3) {
          id[3] <- id[3] + 1
          }else{
          id[3] <- 1
          if (id[2] < 3) {
              id[2] <- id[2] + 1
              }else{ 
              id[2] <- 1
              id[1] <- id[1] + 1
              }
          }
    } ## kind of binary-addition loop
  rgbcol <- rgbcol[order(apply(rgbcol,1,sum)),]
  col <- rgb(rgbcol[-14,]) ## remove middle-grey and generating colors
  lcol <- rgb(0.6,0.6,0.6)
  ##
  ## Loop through all pspectra
  ##
   # stop("again")
  for (pspec in 1:length(pspecIdx)) {
    EIC<-xeic@eic[[pspec]]
    pidx <- object@pspectra[[pspecIdx[pspec]]]    
    peaks <- CAMERA:::getPeaks(object@xcmsSet,object@psSamples[pspecIdx[pspec]])[pidx,]
    grps <- object@groupInfo[pidx,]
    nap <- which(is.na(peaks[,1]))
    naps <- rep(FALSE, nrow(peaks))
    naps[nap] <- TRUE
    peaks[nap,] <- cbind(grps[nap,c(1:6), drop=FALSE],matrix(nrow=length(nap),ncol=5,0))
    main <- paste("Pseudospectrum ", pspecIdx[pspec], sep="")
    ## Calculate EICs and plot ranges
    neics <- length(pidx)
    lmaxlabel <- min(maxlabel, neics)
    ##col <- c((lmaxlabel+1):2, rep(1, neics-lmaxlabel))
    eicidx <- 1:neics;
    maxint <- numeric(length(eicidx))
    for (j in eicidx) {
        maxint[j] <- max(EIC[[j]][,"intensity"])
    }
    rt <- xeic@rtrange[pspec, ];
    ## Open Plot
   # stop("testing idx 4,5")
    plot(0, 0, type = "n", xlim = rt, ylim = c(0, max(maxint)),
                   xlab = "Retention Time (seconds)", ylab = "Intensity",
                   main = paste("Extracted Ion Chromatograms for ", main))
    ## Plot Peak and surrounding
    o <- order(maxint, decreasing = TRUE)
    for (j in eicidx[o]) {
        pts <- xeic@eic[[pspec]][[j]]
        points(pts, type = "l", col = lcol)
        peakrange <- peaks[,c("rtmin","rtmax"), drop=FALSE]
        ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
        if (naps[j]) points(pts[ptsidx,], type = "l", col = col[j], lwd=1.3, lty=3) else
                     points(pts[ptsidx,], type = "l", col = col[j], lwd=1.3)
    }
    ## Plot Annotation Legend
    pspectrum <- getpspectra(object, grp=pspecIdx[pspec])
    if (lmaxlabel>0 & "adduct" %in% colnames(pspectrum)) {
      adduct <- sapply(strsplit(pspectrum[o[1:lmaxlabel], "adduct"], " "),
                       function(x) {if (length(x)>0) x[1] else ""})
      mz <- format(pspectrum[o[1:lmaxlabel], "mz"], digits=5)
      legend("topright", legend=paste(mz, adduct),
             col=col[o], lty=1)
    }
    if (sleep > 0) {
      Sys.sleep(sleep)
    }
  }
})



setGeneric("plotPsSpectrum", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,
                                              sleep=0) standardGeneric("plotPsSpectrum"))
setMethod("plotPsSpectrum", "xsAnnotate", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,
                                              sleep=0)
      {

          ##
          ## Loop through all requested pspectra
          ##
          if(is.na(object@sample)){
              gvals <- groupval(object@xcmsSet);
              peakmat <- object@xcmsSet@peaks;
              groupmat <- groups(object@xcmsSet);
              #errechne höchsten Peaks
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
            }else{
              intensity <- maxo[pindex]
            }
            o <- order(intensity, decreasing=TRUE)
            mz <- pspectrum[o, "mz"]

              if (log) {
                  intensity <- log(intensity[o])
              } else {
                  intensity <- intensity[o]
              }

              mzrange <- range(mz)
              mzborder <- (mzrange[2]-mzrange[1])*0.05
              if(is.null(title)){
                  title <- paste("Plot m/z values of pseudospectrum",psp);
                  plot(mz, intensity, type="h",
                    xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
                    main=title, col="darkgrey")
                  title <- NULL;

              }else{
                  plot(mz, intensity, type="h",
                      xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
                      main=title, col="darkgrey")
              }
              if (maxlabel>0) {
                  ## Also check for labels:  https://stat.ethz.ch/pipermail/r-help/2008-November/178666.html
                  ## There is also the spread.labs function in the TeachingDemos package that
                  ## uses a different method from the plotrix function and should not move any
                  ## labels that are not overlapping.
                lmaxlabel <- min(maxlabel,length(mz))
                  text(mz[1:lmaxlabel],
                       intensity[1:lmaxlabel],
                       labels=format(mz[1:lmaxlabel], digits=5),
                       cex=0.66)

                  if ("adduct" %in% colnames(pspectrum)) {
                      adduct <- sapply(strsplit(pspectrum[o, "adduct"], " "), function(x) {x[1]})
                      text(mz[1:lmaxlabel],
                           intensity[1:lmaxlabel] + strheight("0123456789"),
                           labels=adduct[1:lmaxlabel],
                           cex=0.66)
                  }
              }

              if (sleep > 0)
                  Sys.sleep(sleep)
          }
      })