setGeneric("plotEICs", function(object,
                                xraw,
                                pspecIdx=1:length(object@pspectra),
                                maxlabel=0, sleep=0,
                                ...) standardGeneric("plotEICs"))
setMethod("plotEICs", "xsAnnotate", function(object,
                                             xraw,
                                             pspecIdx=1:length(object@pspectra),
                                             maxlabel=0, sleep=0)
      {

          ## Expand RT plot width by
          ## 1% of runtime:
          rtmargin <- 0.01 * max(xraw@scantime)

          ##
          ## Calculate EICs of requested
          ## PC groups
          ##
          xeic <- new("xcmsEIC")
          xeic@rtrange <- matrix(nrow=length(pspecIdx), ncol=2)
          xeic@mzrange <- matrix(nrow=length(pspecIdx), ncol=2)

          pcpos <- 1
          xeic@eic <- lapply (pspecIdx, function(pc) {
              pidx <- object@pspectra[[pc]]
              bbox <- c(rtmin = min(object@peaks[pidx,"rtmin"])-rtmargin,
                        rtmax = max(object@peaks[pidx,"rtmax"])+rtmargin,
                        mzmin = min(object@peaks[pidx,"mzmin"]),
                        mzmax = max(object@peaks[pidx,"mzmax"]))
              eic <- xcms:::getEIC(xraw,
                            rtrange=matrix(bbox[c("rtmin","rtmax")],
                            nrow=length(pidx), ncol=2, byrow=TRUE),
                            mzrange=object@peaks[pidx,c("mzmin", "mzmax"), drop=FALSE])

              xeic@rtrange[pcpos, ] <<- bbox[c("rtmin","rtmax")]
              xeic@mzrange[pcpos, ] <<- bbox[c("mzmin","mzmax")]
              pcpos <<- pcpos+1
              eic@eic[[1]]
          })
          names(xeic@eic) <- paste("Pseudospectrum ", pspecIdx, sep="")

          ##
          ## Setup Peak and Neighbourhood Colors,
          ## Now: all black
          ## Later: by Annotation ?
          ##
          lcol <- 1
          for (i in seq(along = lcol)) {
              rgbvec <- pmin(col2rgb(lcol[i])+153,255)
              lcol[i] <- rgb(rgbvec[1], rgbvec[2], rgbvec[3], max = 255)
          }

          ##
          ## Loop through all pspectra
          ##
          peaks <- object@peaks
          for (pspec in 1:length(xeic@eic)) {
              ## Calculate EICs and plot ranges
              neics <- length(xeic@eic[[pspec]])
              lmaxlabel <- min(maxlabel, neics)
              col <- c((lmaxlabel+1):2, rep(1, neics-lmaxlabel))
              eicidx <- 1:neics

              maxint <- numeric(length(eicidx))
              for (j in eicidx) {
                  maxint[j] <- max(xeic@eic[[pspec]][[j]][,"intensity"])
              }
              rtrange <- xeic@rtrange[pspec, ]

              ## Open Plot
              plot(0, 0, type = "n", xlim = rtrange, ylim = c(0, max(maxint)),
                   xlab = "Retention Time (seconds)", ylab = "Intensity",
                   main = paste("Extracted Ion Chromatograms for ", names(xeic@eic)[pspec]))

              ## Plot Peak and surrounding
              o <- order(maxint, decreasing = TRUE)
              col[o] <- col
              for (j in eicidx[o]) {
                  pts <- xeic@eic[[pspec]][[j]]
                  points(pts, type = "l", col = lcol[j])
                  peakrange <- peaks[object@pspectra[[pspecIdx[pspec]]], c("rtmin","rtmax"), drop=FALSE]
                  ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
                  points(pts[ptsidx,], type = "l", col = col[j])
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

              if (sleep > 0)
                  Sys.sleep(sleep)
          }
   })


setGeneric("plotPeaks", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,
                                              sleep=0) standardGeneric("plotPeaks"))
setMethod("plotPeaks", "xsAnnotate", function(object, pspec=1:length(object@pspectra), log=FALSE,
                       value="into", maxlabel=0, title=NULL,
                                              sleep=0)
      {

          ##
          ## Loop through all requested pspectra
          ##

          for (psp in pspec) {
            pspectrum <- getpspectra(object, grp=psp)
              intensity <- pspectrum[, value]
              o <- order(intensity, decreasing=TRUE)
              mz <- pspectrum[o, "mz"]

              if (log) {
                  intensity <- log(intensity[o])
              } else {
                  intensity <- intensity[o]
              }

              mzrange <- range(mz)
              mzborder <- (mzrange[2]-mzrange[1])*0.05
              plot(mz, intensity, type="h",
                   xlim=c(max(0,mzrange[1]-mzborder), mzrange[2]+mzborder),
                   main=title, col="darkgrey")

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
