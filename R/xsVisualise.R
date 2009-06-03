setGeneric("plotEICs", function(object,
                                xraw,
                                pspecIdx=1:length(object@pspectra),
                                sleep=0,
                                ...) standardGeneric("plotEICs"))
setMethod("plotEICs", "xsAnnotate", function(object,
                                             xraw,
                                             pspecIdx=1:length(object@pspectra),
                                             sleep=0)
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
                            mzrange=object@peaks[pidx,c("mzmin", "mzmax")])

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
          col <- 1
          lcol <- col
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
              eicidx <- 1:length(xeic@eic[[pspec]])


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
              for (j in eicidx[order(maxint, decreasing = TRUE)]) {
                  pts <- xeic@eic[[pspec]][[j]]
                  points(pts, type = "l", col = lcol[1])
                  peakrange <- peaks[object@pspectra[[pspecIdx[pspec]]], c("rtmin","rtmax")]
                  ptsidx <- pts[,"rt"] >= peakrange[j,1] & pts[,"rt"] <= peakrange[j,2]
                  points(pts[ptsidx,], type = "l", col = col[1])
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
                  text(mz[1:min(maxlabel,length(mz))],
                       intensity[1:min(maxlabel,length(mz))],
                       labels=format(mz[1:min(maxlabel,length(mz))], digits=5),
                       cex=0.66)

                  if ("adduct" %in% colnames(pspectrum)) {
                      adduct <- sapply(strsplit(pspectrum[o, "adduct"], " "), function(x) {x[1]})
                      text(mz[1:min(maxlabel,length(mz))],
                           intensity[1:min(maxlabel,length(mz))] + strheight("0123456789"),
                           labels=adduct[1:min(maxlabel,length(mz))],
                           cex=0.66)
                  }
              }

              if (sleep > 0)
                  Sys.sleep(sleep)
          }
      })
