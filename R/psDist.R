setGeneric("psDist", function(object1,object2, PSpec1,PSpec2, 
                              method=getOption("BioC")$xcms$specDist.method,
                              mzabs=0.001, mzppm=10,symmetric=FALSE,...) standardGeneric("psDist"))
 setMethod("psDist", signature(object1="xsAnnotate",object2="xsAnnotate"), 
                     function(object1,object2, PSpec1,PSpec2, method=getOption("BioC")$xcms$specDist.method, 
                                                                mzabs=0.001, mzppm=10,symmetric=FALSE,...) {
     if (length(object1@pspectra[[PSpec1]])>1) { 
         peakTable1 <- cbind(object1@pspectra[[PSpec1]], object1@groupInfo[object1@pspectra[[PSpec1]],])
         }else{ peakTable1 <- cbind(object1@pspectra[[PSpec1]], t(as.matrix(object1@groupInfo[object1@pspectra[[PSpec1]],])))}
     ints=NA
     gcol <- ncol(object1@groupInfo)-length(sampnames(object1@xcmsSet))+1+object1@psSamples[PSpec1] ## this + 1 is the into-col of the first sample
     peakTable1 <- matrix(ncol=4, data=c(peakTable1[,1],peakTable1[,c("mz","rt")],peakTable1[,gcol]))
     colnames(peakTable1) <- c("oPeak","mz","rt","into")
     if (length(object2@pspectra[[PSpec2]])>1) { 
         peakTable2 <- cbind(object2@pspectra[[PSpec2]], object2@groupInfo[object2@pspectra[[PSpec2]],])
         }else{ peakTable2 <- cbind(object2@pspectra[[PSpec2]], t(as.matrix(object2@groupInfo[object2@pspectra[[PSpec2]],])))}
     ints=NA
     gcol <- ncol(object2@groupInfo)-length(sampnames(object2@xcmsSet))+1+object2@psSamples[PSpec2] ## this + 1 is the into-col of the first sample in xsAnnotate 2
     peakTable2 <- matrix(ncol=4, data=c(peakTable2[,1],peakTable2[,c("mz","rt")],peakTable2[,gcol]))
     colnames(peakTable2) <- c("oPeak","mz","rt","into")
     method <- match.arg(method, getOption("BioC")$xcms$specDist.methods)
     if (is.na(method)) stop("unknown method : ", method)
     method <- paste("specDist", method, sep=".")
     distance <- do.call(method, alist<-list(peakTable1, peakTable2, ...))
     return(distance)
     })

## for the path of single samples - not yet complete implemented

# setGeneric("psDists", function(xanno,method="meanMZmatch", params=NA,rtabs=5,nSlaves=1,samples=NA...) standardGeneric("psDists"))
# setMethod("psDists","xsAnnotate", function(xanno,method="meanMZmatch", params=NA, rtabs=5,nSlaves=1,samples=NA...) {
#    # stop("halooo?")
#     lxs <- length(samples) ## number of samples
#     if (nSlaves>1) mpi.spawn.Rslaves(nslaves=nslaves)
#     lx <- NA ## length of Feature-group list of each xanno    
#     for (a in 1:lxs){
#         lx[a] <- length(xanno@xDataList[[samples[a]]]$pspectra)
#         }
#     ranges <- list() ## all ranges of mz and rt for each pseudospectrum in each sample
#     mzrange <- c(200000,0) ## min/max mz for the whole xanno
#     rtrange <- c(200000,0) ## min/max rt for the whole x
#     for (a in 1:lxs){
#         range<-matrix(ncol=4,nrow=lx[a])
#         for (gfs in 1:lx[a]) {
#             peaks <- xanno@peaks[which(xanno@peaks[,"sample"] == samples[a]),]    
#             range[gfs,] <- c( range(peaks[xanno@xDataList[[samples[a]]]$pspectra[[gfs]],"mz"]),
#                               range(peaks[xanno@xDataList[[samples[a]]]$pspectra[[gfs]],"rt"]) 
#                             )
#             if (range[gfs,1] < mzrange[1]) mzrange[1] <- range[gfs,1]
#             if (range[gfs,2] > mzrange[2]) mzrange[2] <- range[gfs,2]
#             if (range[gfs,3] < rtrange[1]) rtrange[1] <- range[gfs,3]
#             if (range[gfs,4] > rtrange[2]) rtrange[2] <- range[gfs,4]
#             }
#         ranges[[a]] <- range
#         }
#     if (rtabs==0) rtabs=rtrange[2]
# 
#  ### ich weiss noch nicht genau ob ich das mit den ranges brauche
#     cat("mzrange:",mzrange[1]," - ", mzrange[2],"\n")
#     cat("rtrange:",rtrange[1]," - ", rtrange[2],"\n")
#     cat("rtabs:", rtabs)
#     cat("\n")
#     distances <- list()
#     dpos <- 1
#     commands <- list()
#     A <- NA
#     B <- NA
#     for (s in 1:lxs) {## for each sample s
#         if (s+1<=lxs)
#         for (rs in (s+1) : lxs) { ## for each remeaning sample rs
#             command <- list(peaks1=xanno@peaks[which(xanno@peaks[,"sample"]==samples[s]),],
#                             peaks2=xanno@peaks[which(xanno@peaks[,"sample"]==samples[rs]),],
#                             pspectra1=xanno@xDataList[[xanno@samples[samples[s]]]]$pspectra,
#                             pspectra2=xanno@xDataList[[xanno@samples[samples[rs]]]]$pspectra,
#                             ranges1=ranges[[s]],ranges2=ranges[[rs]],
#                             rtabs=rtabs,mzmethod=method,mzparams=params)
#             commands[[dpos]] <- command
#             A[dpos] <- s
#             B[dpos] <- rs
#             dpos <- dpos + 1            
#             }
#         }
#     if (nSlaves>1) {
#         distances <- xcmsPapply(commands,sdist)
#         }else{
#         for (a in 1:length(commands)) {
#             distances[[a]]<-sdist(commands[[a]]) ## for debug
#             }
#         }
#     if (nSlaves>1) mpi.close.Rslaves()
#     distances$A <- A
#     distances$B <- B
#     distances
#     })
# 
# ## helperfunctions
# 
# getPspec <- function(xanno, index, PS){
#     return(xanno@peaks[which(xanno@peaks[,"sample"]==xanno@samples[index]),][xanno@xDataList[[index]]$pspectra[[PS]],])
#     }
# 
# sdist <- function(command){## compares each pair of PS in two given samples
#     peaks1<-command$peaks1##peaks of first sample
#     peaks2<-command$peaks2##peaks of second sample
#     pspectra1 <- command$pspectra1 ## PSlist of first sample
#     pspectra2 <- command$pspectra2 ## PSlist of second sample
#     rtabs <-command$rtabs## maximum rt distance for two PS
#     mzmethod<-command$mzmethod ## mz distance method
#     mzparams<-command$mzparams ## parameter for mzmethod
#     lx1 <- length(pspectra1)
#     lx2 <- length(pspectra2)
#     
#    # stop("here")
#     ## mz/rt ranges for sample1 (a matrix with 4x #pspectra )
#     minMZ1 <-command$ranges1[,1];maxMZ1 <-command$ranges1[,2]
#     minRT1 <-command$ranges1[,3];maxRT1 <-command$ranges1[,4]
#     minMZ2 <-command$ranges2[,1];maxMZ2 <-command$ranges2[,2]
#     minRT2 <-command$ranges2[,3];maxRT2 <-command$ranges2[,4]
#     
#     distances <- matrix(nrow = lx1, ncol = lx2)
#     for (fgs in 1:lx1) { ## for each FG in sample s
#         for (fgrs in 1:lx2) { ## for each FG in sample rs
#             if ((maxMZ1[fgs] >= minMZ2[fgrs]) &
#                 (minMZ1[fgs] <= maxMZ2[fgrs])) {
#                 if ( (maxRT1[fgs]+rtabs > minRT2[fgrs]) &
#                      (minRT1[fgs]-rtabs < maxRT2[fgrs])) {
#                     if (is.na(mzparams)) distances[fgs,fgrs] <-do.call(mzmethod, c(alist(peaks1[pspectra1[[fgs]],],peaks2[pspectra2[[fgrs]],])))
#                     else                 distances[fgs,fgrs] <-do.call(mzmethod, c(alist(peaks1[pspectra1[[fgs]],],peaks2[pspectra2[[fgrs]],]),mzparams))
#                     }
#                 }
#             }
#         }
#     distances   
#     }