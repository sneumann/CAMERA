xsAnnotate <- function(xs=NULL, sample=NA, nSlaves=1, polarity=NULL){

 ## sample is:
 ### NA for maxInt-way
 ##  1-nSamp for manual way
 ##  -1 for specDist way (to be implemented)
  if(is.null(xs)) { stop("no argument was given"); }
  else if(!class(xs)=="xcmsSet") stop("xs is no xcmsSet object ");
  
  object  <- new("xsAnnotate");
  
  if(length(sampnames(xs)) > 1 && !nrow(xs@groups) > 0) {  
      # more than one sample
      # checking alignment
      stop ('First argument must be a xcmsSet with group information or contain only one sample.') 
  }

  # check multiple sample selection
  if(length(sample) > 1){
    #example sample <- c(1,2,3,4,5)
  
    if(all(sample < length(xs@filepaths) && sample > 0)){
      #all sample values are in allowed range
      object@sample <- sample;
    }else{
        #one or more sample values are lower 0 or higher than number of samples (length(xs@filepaths))
        stop("All values in parameter sample must be lower equal the number of samples and greater than 0.\n")
      }
  }else{
    if(is.null(sample) || is.na(sample)) {
      #automatic sample selection, sample = NA
      if(length(xs@filepaths) == 1){
          #If samplesize == 1 than set sample to 1
          object@sample <- 1;
      }else{
          object@sample   <-  as.numeric(NA);
      }
    }else{
      if(sample == -1){
        #Joes Way
        object@sample <-  sample;
      }else if(length(xs@filepaths) < sample | sample < 1) {
        stop("Parameter sample must be lower equal than number of samples and greater than 0.\n")
      }else{       
        object@sample <-  sample;
      }
    }
  }

  object@groupInfo <- getPeaks_selection(xs);
  runParallel <- list()
  runParallel$enable   <-  0;

  if (nSlaves > 1) {
    ## If MPI is available ...
    rmpi = "Rmpi"
    opt.warn <- options("warn")$warn
    options("warn" = -1) 
    if (require(rmpi,character.only=TRUE) && !is.null(nSlaves)) {
      if (is.loaded('mpi_initialize')) {
        #test if not already slaves are running!
        if(mpi.comm.size() >0){ 
          warning("There are already intialized mpi slaves on your machine.\nCamera will try to uses them!\n");
          runParallel$enable <-1;
          runParallel$mode <- rmpi;
        }else{
          mpi.spawn.Rslaves(nslaves=nSlaves, needlog=FALSE)
          if(mpi.comm.size() > 1){
            #Slaves have successfull spawned
            runParallel$enable <-1;
            runParallel$mode <- rmpi;
          }else{ warning("Spawning of mpi slaves have failed. CAMERA will run without parallelization.\n");}
        }
      }else {
          #And now??
          warning("DLL mpi_initialize is not loaded. Run single core mode!\n");
      }
    } else {
      #try local sockets using snow package
      snow = "snow"
      if (try(require(snow,character.only=TRUE,quietly=TRUE))) {
        cat("Starting snow cluster with",nSlaves,"local sockets.\n")
        snowclust <- makeCluster(nSlaves, type = "SOCK")
        runParallel$enable <- 1
        runParallel$mode <- snow;
        runParallel$cluster <- snowclust
      }
    }
    options("warn" = opt.warn)
    cat("Run cleanParallel after processing to remove the spawned slave processes!\n")
  }

  if(!is.null(polarity)){
    if(is.na(match.arg(polarity, c("positive","negative")))){
      stop("Parameter polarity is unknown: ", graphMethod,"\n")  
    }else{
      object@polarity <- polarity;
    }
  }
  object@runParallel<-runParallel;

  colnames(object@annoID) <-  c("id","grpID","ruleID","parentID");
  colnames(object@annoGrp)<-  c("id","mass","ips","psgrp");
  colnames(object@isoID)  <-  c("mpeak","isopeak","iso","charge")

  #save xcmsSet in the xsAnnotate object
  object@xcmsSet  <-  xs;
  return(object);
}

setMethod("show", "xsAnnotate", function(object){
  #show main information
  cat("An \"xsAnnotate\" object!\n");
  cat("With",length(object@pspectra),"groups (pseudospectra)\n");
  if(!is.null(object@xcmsSet)){
    cat("With",length(sampnames(object@xcmsSet)),"samples and",nrow(object@groupInfo),"peaks\n");
  }else{
    cat("Include no xcmsSet set\n");
  }
    
  #Show polarity if avaiable
  if(length(object@polarity) > 0){
    cat("Polarity mode is set to: ",object@polarity,"\n");
  }
    
  #Show samples selection
  if(is.null(object@sample)){
    cat("Parameter sample not set\n");
  }else{
    if(is.na(object@sample[1])){
      cat(paste("Using automatic sample selection\n"));
    } else if(all(object@sample > -1)){
      cat("Using sample(s): ",paste(object@sample),"\n");
    } else { 
      cat(paste("Using complete measurement\n"));
    }
  }
  
  #Show isotope information
  if(length(object@isotopes) > 0){
    cnt <- nrow(object@isoID)
    cat("Annotated isotopes:", cnt, "\n");
  }

  #Show annotation information
  if(length(object@derivativeIons) > 0){
    cnt <- length(unique(object@annoID[, 1]));
    cat("Annotated adducts & fragments:", cnt, "\n");
  }

  #Show memory information
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n");
  if(!is.null(object@runParallel) && object@runParallel$enable == 1){
    cat("CAMERA runs in parallel mode!\n");
  }
})

###End Constructor###

setGeneric("groupComplete", function(object,...)
  standardGeneric("groupComplete"))

setMethod("groupComplete", "xsAnnotate", function(object, h=NULL,...) {
  
  sample    <- object@sample;
  pspectra  <- list();
  psSamples <- NA;
  
  #generate distance matrix
  nPeaks <- nrow(object@groupInfo)
  
  distMat <- matrix(NA, nrow=nPeaks, ncol=nPeaks)

  # First Part of Score
  # Remove peaks combination which are to far away
  
  rt <- object@xcmsSet@peaks[, c("rtmax","rtmin","rt"),drop=FALSE]
  max.Rt <-  max((rt[, 1]-rt[, 3]), (rt[, 3]-rt[, 2]))
  distRT <- dist(object@groupInfo[,"rt"], method="manhattan")
  if(is.null(h)){
    result <- cutree(hclust(distRT), h=max.Rt/2)    
  }else if(is.numeric(h)){
    result <- cutree(hclust(distRT), h=h)
  } else {
    stop("h must be numeric!\n")
  }
  

  for(i in unique(result)){
    index <- which(result == i)
    pspectra[[length(pspectra)+1]] <- index
  }
  psSamples <- rep(sample, length(pspectra))

  object@pspectra  <- pspectra;
  object@psSamples <- psSamples;
  cat("Created", length(object@pspectra), "pseudospectra.\n")
  return(object)
})

setGeneric("groupDen", function(object,...)
  standardGeneric("groupDen"))

setMethod("groupDen", "xsAnnotate", function(object, bw=5, ...) {
  
  sample    <- object@sample;
  pspectra  <- list();
  psSamples <- NA;

  #generate distance matrix
  nPeaks <- nrow(object@groupInfo)
  rt <- object@groupInfo[ ,"rt"]  
  
  den <- density(rt, bw, from = min(rt)-3*bw, to = max(rt)+3*bw)
  maxden <- max(den$y)
  deny <- den$y
  snum <- 0
  while (deny[maxy <- which.max(deny)] > 0 ) {
    grange <- xcms:::descendMin(deny, maxy)
    deny[grange[1]:grange[2]] <- 0
    gidx <- which(rt >= den$x[grange[1]] & rt <= den$x[grange[2]])
    pspectra[[length(pspectra)+1]] <- gidx
    snum <- snum + length(gidx)
  }
  psSamples <- rep(sample, length(pspectra))
  
  object@pspectra  <- pspectra;
  object@psSamples <- psSamples;
  cat("Created", length(object@pspectra), "pseudospectra.\n")
  return(object)
})




###xsAnnotate generic Methods###
setGeneric("groupFWHM", function(object, sigma=6, perfwhm=0.6, intval="maxo") 
  standardGeneric("groupFWHM"))

setMethod("groupFWHM", "xsAnnotate", function(object, sigma=6, perfwhm=0.6, intval="maxo") {
  # grouping after retentiontime 
  # sigma - number of standard deviation arround the mean (6 = 2 x 3 left and right)
  # perfwhm - 0.3;

  if (!class(object) == "xsAnnotate") {
    stop ("no xsAnnotate object")
  }

  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!")
  }

  sample    <- object@sample;
  pspectra  <- list();
  psSamples <- NA;

  cat("Start grouping after retention time.\n")

  if(object@groupInfo[1, "rt"] == -1) {
     # Like FTICR Data
     warning("Warning: no retention times avaiable. Do nothing\n");
     return(invisible(object));
  }else{
    if(is.na(sample[1]) || length(object@xcmsSet@filepaths) > 1) {
      # grouped peaktable within automatic selection or sub selection
      if(is.na(sample[1])){
        index <- 1:length(object@xcmsSet@filepaths);
      }else{
        index <- sample;
      }
      gvals    <- groupval(object@xcmsSet)[,index,drop=FALSE];
      peakmat  <- object@xcmsSet@peaks;
      groupmat <- groups(object@xcmsSet);

      #calculate highest peaks
      maxo      <- as.numeric(apply(gvals, 1, function(x, peakmat){
        val <- na.omit(peakmat[x, intval]);
        if(length(val) == 0){
          return(NA);
        }else{
          return(max(val))
        }
      }, peakmat));
      
      maxo[which(is.na(maxo))] <- -1;
      maxo      <- cbind(1:length(maxo),maxo);
      
      #highest peak index 
      int.max   <- as.numeric(apply(gvals, 1, function(x, peakmat){which.max(peakmat[x, intval])}, peakmat));

      peakrange <- matrix(apply(gvals, 1, function(x, peakmat) { 
        val <- peakmat[x, intval];
        if(length(na.omit(val)) == 0){
          return(c(0,1));
        } else {
          return(peakmat[x[which.max(val)], c("rtmin", "rtmax")]);
        }
      }, peakmat), ncol=2, byrow=TRUE); 

      colnames(peakrange) <- c("rtmin", "rtmax")

      while(length(maxo) > 0){
          iint   <- which.max(maxo[,2]);
          rtmed  <- groupmat[iint, "rtmed"]; #highest peak in whole spectra
          rt.min <- peakrange[iint, "rtmin"];
          rt.max <- peakrange[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt.max-rt.min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
          #all other peaks whose retensiontimes are in the fwhm of the highest peak
          irt    <- which(groupmat[, 'rtmed'] > (rtmed-hwhm) & groupmat[, 'rtmed'] < (rtmed + hwhm)) 
          if(length(irt) > 0){
              #if peaks are found
              idx <- maxo[irt,1];
              pspectra[[length(pspectra)+1]] <- idx; #create groups
              psSamples[length(pspectra)]  <- index[int.max[maxo[iint,1]]] # saves the sample of the peak which is in charge for this pspectrum
              maxo <- maxo[-irt, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
              groupmat   <- groupmat[-irt, ,drop=FALSE];
              peakrange  <- peakrange[-irt, ,drop=FALSE];
          }else{              
              idx <- maxo[iint,1];
              cat("Warning: Feature ",idx," looks odd for at least one peak. Please check afterwards.\n");
              pspectra[[length(pspectra)+1]] <- idx; #create groups
              psSamples[length(pspectra)]  <- index[int.max[maxo[iint,1]]] # saves the sample of the peak which is in charge for this pspectrum
              maxo <- maxo[-iint, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
              groupmat   <- groupmat[-iint, ,drop=FALSE];
              peakrange  <- peakrange[-iint, ,drop=FALSE];
          }
      }
    }else{
      #One sample experiment
      peakmat <- object@xcmsSet@peaks;
      maxo    <- peakmat[, intval]; #max intensities of all peaks
      maxo    <- cbind(1:length(maxo),maxo);

      while(length(maxo)> 0){
          iint   <- which.max(maxo[,2]);
          rtmed  <- peakmat[iint, "rt"]; #highest peak in whole spectra
          rt.min <- peakmat[iint, "rtmin"];
          rt.max <- peakmat[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt.max - rt.min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
          #all other peaks whose retensiontimes are in the fwhm of the highest peak
          irt    <- which(peakmat[, 'rt'] > (rtmed - hwhm) & peakmat[, 'rt'] < (rtmed + hwhm)) 
          if(length(irt)>0){
              #if peaks are found
              idx <- maxo[irt,1];
              pspectra[[length(pspectra)+1]] <- idx; #create groups
              maxo <- maxo[-irt, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
              peakmat <- peakmat[-irt, ,drop=FALSE];
          }else{
              idx <- maxo[iint,1];
              cat("Warning: Feature ",idx," looks odd for at least one peak. Please check afterwards.\n");
              pspectra[[length(pspectra)+1]] <- idx; #create groups
              maxo       <- maxo[-iint, ,drop=FALSE]; #set itensities of peaks to NA, due to not to be found in the next cycle
              peakmat  <- peakmat[-iint, ,drop=FALSE];
          }
      }
      psSamples <- rep(sample, length(pspectra))
    }

    object@pspectra  <- pspectra;
    object@psSamples <- psSamples;
    cat("Created", length(object@pspectra), "pseudospectra.\n")
  }
  return(invisible(object)); #return object
})


setGeneric("groupCorr",function(object, cor_eic_th=0.75, pval=0.05, graphMethod="hcs", 
                                calcIso = FALSE, calcCiS = TRUE, calcCaS = FALSE, psg_list=NULL, 
                                xraw=NULL, cor_exp_th=0.75, ...) standardGeneric("groupCorr"));

setMethod("groupCorr","xsAnnotate", function(object, cor_eic_th=0.75, pval=0.05, 
                                             graphMethod="hcs", calcIso = FALSE, calcCiS = TRUE, 
                                             calcCaS = FALSE, psg_list=NULL, xraw=NULL,
                                             cor_exp_th=0.75) {
  
  if (!is.numeric(cor_eic_th) || cor_eic_th < 0 || cor_eic_th > 1){
    stop ("Parameter cor_eic_th must be numeric and between 0 and 1.\n");
  }

  if (!is.numeric(cor_exp_th) || cor_exp_th < 0 || cor_exp_th > 1){
    stop ("Parameter cor_exp_th must be numeric and between 0 and 1.\n");
  }
  
  if (!is.numeric(pval) || pval < 0 || pval > 1){
    stop ("Parameter pval must be numeric and between 0 and 1.\n");
  }
  
  
  checkMethod <- match.arg(graphMethod, c("hcs","lpc"))
  if (is.na(checkMethod)){
    stop("Parameter graphMethod is unknown: ", graphMethod,"\n")
  }
  rm(checkMethod)
  
  if (!is.logical(calcIso)) {
    stop ("Parameter calcIso must be logical.\n");
  }
  
  if (!is.logical(calcCiS)) {
    stop ("Parameter calcCiS must be logical.\n");
  }
  
  if (!is.logical(calcCaS)) {
    stop ("Parameter calcCaS must be logical.\n");
  }
  
  if (!is.null(psg_list) && (!is.numeric(psg_list) || any(psg_list > length(object@pspectra)) || any(psg_list < 0)) ) {
    stop ("If parameter psg_list is used, it must be numeric and contains only indices lower or equal maximum number
          of pseudospectra.\n");
  }
  
  if (!is.null(xraw) &&  !class(xraw) == "xcmsRaw") {
    stop ("Parameter xraw must be null or an xcmsRaw object\n");
  }
  
  npspectra <- length(object@pspectra);

  cat("Start grouping after correlation.\n")
  #Data is not preprocessed with groupFWHM 
  if(npspectra < 1){
    cat("Data was not preprocessed with groupFWHM, creating one pseudospectrum with all peaks.\n")
    #Group all peaks into one group
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    if(is.na(object@sample[1])){
      object@psSamples <- rep(1,nrow(object@groupInfo)); ##TODO: Change if sample=NA or sample=number
    }else{
      object@psSamples <- rep(object@sample,nrow(object@groupInfo));
    }
  }

  #save number of pspectra before groupCorr
  cnt <- length(object@pspectra);
  res <- list();

  # Check LC information and calcCorr was selected
  if(calcCiS && object@xcmsSet@peaks[1,"rt"] != -1){
    if(!is.null(xraw)) {
      #Use provided xcmsRaw
      maxscans <- length(xraw@scantime)
      scantimes<-list();
      scantimes[[1]] <- xraw@scantime
      pdata <- object@groupInfo
      EIC <- CAMERA:::getEICs(xraw, pdata, maxscans)
      
      res[[1]] <- calcCiS(object, EIC=EIC, corval=cor_eic_th, 
                          pval=pval, psg_list=psg_list);
    }else if(is.na(object@sample[1])){
      #Autoselect sample path for EIC correlation    
      index <- rep(0, nrow(object@groupInfo));
      
      for(i in 1:npspectra){
        index[object@pspectra[[i]]] <- object@psSamples[[i]];
      }

      #Generate EIC data
      tmp <- getAllPeakEICs(object, index=index);
      EIC <- tmp$EIC
      scantimes <- tmp$scantimes
      rm(tmp);

      res[[1]] <- calcCiS(object, EIC=EIC, corval=cor_eic_th, pval=pval, psg_list=psg_list);

    } else {
      #Calculate EIC-Correlation for selected sample(s)

      for(i in object@sample){
        index <- rep(i, nrow(object@groupInfo));
        tmp <- getAllPeakEICs(object, index=index);
        EIC <- tmp$EIC
        scantimes <- tmp$scantimes
        rm(tmp);
        #new calcCL function with rcorr from Hmisc
        if(length(res) == 0){
          res[[1]] <- calcCiS(object, EIC=EIC, corval=cor_eic_th, pval=pval, psg_list=psg_list);
        }else{
          res[[1]] <- combineCalc(res[[1]],calcCiS(object, EIC=EIC, corval=cor_eic_th, pval=pval, psg_list=psg_list),method="sum")
        }
      }
      res[[1]][,3] <- res[[1]][,3] / length(object@sample)
    }
  } else if(calcCiS && object@xcmsSet@peaks[1,"rt"] == -1){
    cat("Object contains no retention time data!\n");
  }
  
  # Check if sample size > 3 and calcCaS was selected
  if( length(object@xcmsSet@filepaths) > 3 && calcCaS){
    res[[length(res)+1]] <- calcCaS(object, corval=cor_exp_th, pval=pval);
  }else if(length(object@xcmsSet@filepaths) <= 3 && calcCaS){
    cat("Object has to contain more than 3 samples to calculate correlation accros samples!\n");
  }
  
  #If object has isotope information and calcIso was selected
  if( nrow(object@isoID) > 0 && calcIso){
    res[[length(res)+1]] <- calcIsotopes(object);
  }else if(nrow(object@isoID) == 0 && calcIso){
    cat("Object contains no isotope or isotope annotation!\n");
  }

  #Check if we have at least 2 result matrixes
  if(length(res) > 2){
    #combine the first two to create the result Table
    resMat <- combineCalc(res[[1]], res[[2]], method="sum");
    for( i in 3:length(res)){
      resMat <- combineCalc(resMat, res[[i]], method="sum");
    }         
  }else if(length(res) == 2){
    #combine one time
    resMat <- combineCalc(res[[1]], res[[2]], method="sum")
  } else {
    #Only one matrix
    resMat <- res[[1]];
  }

  #if(nrow(resMat) < 1){
    #Matrix contains no edge
    #Do nothing!
   # cat("No group was seperated.\n")
  #  return(invisible(object));
  #}

  #Perform graph seperation to seperate co-eluting pseudospectra
  object <- calcPC(object, method=graphMethod, ajc=resMat, psg_list=psg_list);                                   

  #Create pc groups based on correlation results
  cat("xsAnnotate has now", length(object@pspectra), "groups, instead of", cnt, "\n"); 

  return(invisible(object));
})


setGeneric("findIsotopes", function(object, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, 
                                    intval=c("maxo","into","intb"),minfrac=0.5, 
                                    isotopeMatrix=NULL,filter=TRUE) standardGeneric("findIsotopes"));

setMethod("findIsotopes", "xsAnnotate", 
  function(object, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, 
                   intval=c("maxo","into","intb"), minfrac=0.5,  isotopeMatrix=NULL,
           filter=TRUE){
  
  #searches in every pseudospectrum after mass differences, 
  #which matches isotope distances
  
  ####Test arguments####
  #test maxcharge  
  if(!is.wholenumber(maxcharge) || maxcharge < 1){
    stop("Invalid argument 'maxcharge'. Must be integer and > 0.\n")
  }
  
  #test maxiso
  if(!is.wholenumber(maxiso) || maxiso < 1){
    stop("Invalid argument 'maxiso'. Must be integer and > 0.\n")
  }
  
  #test ppm
  if(!is.numeric(ppm) || ppm < 0){
    stop("Invalid argument 'ppm'. Must be numeric and not negative.\n")
  }
  
  #test mzabs
  if(!is.numeric(mzabs) || mzabs < 0){
    stop("Invalid argument 'mzabs'. Must be numeric and not negative.\n")
  }
  
  #test intval
  intval <- match.arg(intval)
  
  #test minfrac
  if(!is.numeric(minfrac) || minfrac < 0 || minfrac > 1){
    stop("Invalid argument 'minfrac'. Must be numeric and between 0 and 1.\n")
  }
  
  #test isotopeMatrix
  if(!is.null(isotopeMatrix)){
    if(!is.matrix(isotopeMatrix) || ncol(isotopeMatrix) != 4 || nrow(isotopeMatrix) < 1
       || !is.numeric(isotopeMatrix)){
      stop("Invalid argument 'isotopeMatrix'. Must be four column numeric matrix.\n")
    } else {
      colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
    }
  }else if(maxiso > 8){
     stop("Invalid argument 'maxiso'. Must be lower 9 or provide your own isotopeMatrix.\n")
  }else{
    isotopeMatrix <- calcIsotopeMatrix(maxiso=maxiso)
  }
  ####End Test arguments####
  
  npeaks.global <- 0; #Counter for % bar
  npspectra <- length(object@pspectra);

  # scaling
  devppm <- ppm / 1000000;
  filter <- filter;
  #generate parameter list
  params <- list(maxiso=maxiso, maxcharge=maxcharge, devppm=devppm, mzabs=mzabs, IM=isotopeMatrix, minfrac=minfrac, filter=filter)
  
  #Check if object have been preprocessed with groupFWHM
  if(npspectra < 1) {
    cat("xsAnnotate contains no pseudospectra. Regroup all peaks into one!\n")
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    object@psSamples  <- 1;
  }
  
  #number of peaks in pseudospectra
  ncl <- sum(sapply(object@pspectra, length));

  # get mz,rt and intensity values from peaktable
  if(nrow(groups(object@xcmsSet)) > 0){
    ##multiple sample or grouped single sample
    if(is.na(object@sample[1])){
      index <- 1:length(object@xcmsSet@filepaths);
    }else{
      index <- object@sample;
    }
    cat("Generating peak matrix!\n");
    mint     <- groupval(object@xcmsSet,value=intval)[,index,drop=FALSE];
    imz <- object@groupInfo[, "mz", drop=FALSE];
    irt <- object@groupInfo[, "rt", drop=FALSE];
  }else{
    ##one sample case
    cat("Generating peak matrix!\n");
    imz  <- object@groupInfo[, "mz", drop=FALSE];
    irt  <- object@groupInfo[, "rt", drop=FALSE];
    mint <- object@groupInfo[, intval, drop=FALSE];      
  }
  
  isotope   <- vector("list", length(imz));
    
  isomatrix <- matrix(ncol=5, nrow=0);
  colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")


  cat("Run isotope peak annotation\n % finished: ");
  lp <- -1;

  #look for isotopes in every pseudospectra
  for(i in seq(along = object@pspectra)){
    #get peak indizes for i-th pseudospectrum
    ipeak <- object@pspectra[[i]];

    #Ouput counter
    percentOutput(npeaks.global, length(ipeak), ncl, lp)
    
    #Pseudospectrum has more than one peak
    if(length(ipeak) > 1){
      #peak mass and intensity for pseudospectrum
      mz  <- imz[ipeak];
      int <- mint[ipeak, , drop=FALSE];
      isomatrix <-  findIsotopesPspec(isomatrix, mz, ipeak, int, params)              
    }
  }

  #clean isotopes
  if(is.null(nrow(isomatrix))) {
    isomatrix = matrix(isomatrix, byrow=F, ncol=length(isomatrix)) 
  }
  
  #check if every isotope has only one annotation
  if(length(idx.duplicated <- which(duplicated(isomatrix[, 2]))) > 0){
    peak.idx <- unique(isomatrix[idx.duplicated, 2]);
    for( i in 1:length(peak.idx)){
      #peak.idx has two or more annotated charge
      #select the charge with the higher cardinality
      peak <- peak.idx[i];
      peak.mono.idx <- which(isomatrix[,2] == peak)
      if(length(peak.mono.idx) < 2){
        #peak has already been deleted
        next;
      }
      peak.mono <- isomatrix[peak.mono.idx,1]
      #which charges we have
      charges.list   <- isomatrix[peak.mono.idx, 4];
      tmp <- cbind(peak.mono,charges.list);
      charges.length <- apply(tmp,1, function(x,isomatrix) { 
        length(which(isomatrix[, 1] == x[1] & isomatrix[,4] == x[2])) }, 
                              isomatrix);
      idx <- which(charges.length == max(charges.length));
      if(length(idx) == 1){
        #max is unique
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }else{
        #select this one, which lower charge
        idx <- which.min(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }
    }
  }


  #check if every isotope in one isotope grp, have the same charge
  if(length(idx.duplicated <- which(duplicated(paste(isomatrix[, 1], isomatrix[, 3])))) > 0){
    #at least one pair of peakindex and number of isotopic peak is identical
    peak.idx <- unique(isomatrix[idx.duplicated,1]);
    for( i in 1:length(peak.idx)){
      #peak.idx has two or more annotated charge
      #select the charge with the higher cardinality
      peak <- peak.idx[i];
      #which charges we have
      charges.list   <- unique(isomatrix[which(isomatrix[, 1] == peak), 4]);
      #how many isotopes have been found, which this charges
      charges.length <- sapply(charges.list, function(x,isomatrix,peak) { length(which(isomatrix[, 1] == peak & isomatrix[, 4] == x)) },isomatrix,peak);
      #select the charge which the highest cardinality
      idx <- which(charges.length == max(charges.length));
      if(length(idx) == 1){
        #max is unique
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }else{
        #select this one, which lower charge
        idx <- which.min(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }
    }
  }

  #Combine isotope cluster, if they overlap
  index2remove <- c();

  if(length(idx.duplicated <- which(isomatrix[, 1] %in% isomatrix[, 2]))>0){
    for(i in 1:length(idx.duplicated)){
      index <-  which(isomatrix[, 2] == isomatrix[idx.duplicated[i], 1])
      index2 <- sapply(index, function(x, isomatrix) which(isomatrix[, 1] == isomatrix[x, 1] & isomatrix[,3] == 1),isomatrix)
      if(length(index2) == 0){
        index2remove <- c(index2remove,idx.duplicated[i])
      }
      max.index <- which.max(isomatrix[index,4]);
      isomatrix[idx.duplicated[i], 1] <- isomatrix[index[max.index], 1];
      isomatrix[idx.duplicated[i], 3] <- isomatrix[index[max.index], 3]+1;
    }
  }

  if(length(index <- which(isomatrix[,"iso"] > maxiso)) > 0){
    index2remove <- c(index2remove, index)
  }
  
  if(length(index2remove) > 0){
    isomatrix <- isomatrix[-index2remove,, drop=FALSE];
  }

  isomatrix <- isomatrix[order(isomatrix[,1]),,drop=FALSE]
  #Create isotope matrix within object
  object@isoID <- matrix(nrow=0, ncol=4);
  colnames(object@isoID)  <-  c("mpeak", "isopeak", "iso", "charge");
  
  #Add isomatrix to object
  object@isoID <- rbind(object@isoID, isomatrix[, 1:4]);
  
  # counter for isotope groups
  globalcnt <- 0;
  oldnum    <- 0;
  
  if(nrow(isomatrix) > 0){
    for( i in 1:nrow(isomatrix)){
      if(!isomatrix[i, 1] == oldnum){
          globalcnt <- globalcnt+1; 
          isotope[[isomatrix[i, 1]]] <- list(y=globalcnt, iso=0, charge=isomatrix[i, 4], val=isomatrix[i, 5]);
          oldnum <- isomatrix[i, 1];
      };
      isotope[[isomatrix[i,2]]] <- list(y=globalcnt,iso=isomatrix[i,3],charge=isomatrix[i,4],val=isomatrix[i,5]);
    }
  }
  cnt<-nrow(object@isoID);
  cat("\nFound isotopes:",cnt,"\n");
  object@isotopes <- isotope;
  return(object);
})

setGeneric("findAdducts", function(object, ppm=5, mzabs=0.015, multiplier=3, polarity=NULL, 
                                   rules=NULL, max_peaks=100, psg_list=NULL) standardGeneric("findAdducts"));
setMethod("findAdducts", "xsAnnotate", function(object, ppm=5, mzabs=0.015, multiplier=3, polarity=NULL, 
                                                rules=NULL, max_peaks=100, psg_list=NULL){
  newFragments <- FALSE;
  
  # Scaling ppm factor
  devppm <- ppm / 1000000;
  # counter for % bar
  npeaks.global <- 0;

  # get mz values from peaklist
  imz    <- object@groupInfo[, "mz"];
  #TODO: Change intensity choosing?
  intval <- "maxo"; #max. intensity

  #number of pseudo-spectra
  npspectra <- length(object@pspectra);

  #If groupCorr or groupFHWM have not been invoke, select all peaks in one sample
  if(npspectra < 1){ 
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo)); 
  }

  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##one sample case
    mint <- object@groupInfo[, intval];
  }else{
    ##multiple sample
    if(is.na(object@sample[1])){
      index <- 1:length(object@xcmsSet@filepaths);
    }else{
      index <- object@sample;
    }
    
    cat("Generating peak matrix for peak annotation!\n");
    mint     <- groupval(object@xcmsSet,value=intval)[, index, drop=FALSE];
    peakmat  <- object@xcmsSet@peaks;
    groupmat <- groups(object@xcmsSet);
    
    imz <- groupmat[, "mzmed"];
    irt <- groupmat[, "rtmed"];
    int.val <- c();
    nsample <- length(object@sample);
  }


  # isotopes
  isotopes  <- object@isotopes;
  # adductlist
  derivativeIons <- vector("list", length(imz));
  
  # other variables
  oidscore <- c();
  index    <- c();
  annoID   <- matrix(ncol=4, nrow=0)
  annoGrp  <- matrix(ncol=4, nrow=0)
  colnames(annoID)  <-  colnames(object@annoID)
  colnames(annoGrp) <-  colnames(object@annoGrp)

  ##Examine polarity and rule set
  if(!(object@polarity == "")){
    cat(paste("Polarity is set in xsAnnotate:", object@polarity, "\n"));
    if(is.null(rules)){
      if(!is.null(object@ruleset)){
        rules <- object@ruleset;
      }else{ 
        cat("Ruleset could not read from object! Recalculate\n");
        rules <- calcRules(maxcharge=3, mol=3, nion=2, nnloss=1, nnadd=1, nh=2,
                           polarity=object@polarity, 
                           lib.loc= .libPaths(),newFragments=newFragments);
        object@ruleset <- rules;
      }
    }else{ 
      object@ruleset <- rules;
      cat("Found and use user-defined ruleset!");
    }
  }else{
    if(!is.null(polarity)){
      if(polarity %in% c("positive","negative")){
        if(is.null(rules)){
          rules <- calcRules(maxcharge=3, mol=3, nion=2, nnloss=1, nnadd=1, 
                              nh=2, polarity=polarity, lib.loc= .libPaths(),
                             newFragments=newFragments);
        }else{ cat("Found and use user-defined ruleset!");}
          object@polarity <- polarity;
      }else stop("polarity mode unknown, please choose between positive and negative.")
    }else if(length(object@xcmsSet@polarity) > 0){
      index <- which(sampclass(object@xcmsSet) == object@category)[1] + object@sample-1;
      if(object@xcmsSet@polarity[index] %in% c("positive","negative")){
        if(is.null(rules)){
          rules <- calcRules(maxcharge=3, mol=3, nion=2, nnloss=1, nnadd=1, 
                             nh=2, polarity=object@xcmsSet@polarity[index], 
                             lib.loc= .libPaths(), newFragments=newFragments);
        }else{ cat("Found and use user-defined ruleset!");}
        object@polarity <- polarity;
      }else stop("polarity mode in xcmsSet unknown, please define variable polarity.")
  }else stop("polarity mode could not be estimated from the xcmsSet, please define variable polarity!")
    #save ruleset
    object@ruleset <- rules;
  }

  ##Run as single or parallel mode
  runParallel <- 0;

  if(object@runParallel$enable == 1){
    if(!(is.null(object@runParallel$cluster)) || mpi.comm.size() > 0 ){
      runParallel <- 1;
    }else{
      warning("CAMERA runs in parallel mode, but no slaves are spawned!\nRun in single core mode!\n");
      runParallel <- 0;
    }
  }else{
    runParallel <- 0;
  }

  if("quasi" %in% colnames(rules)){
  #backup for old rule sets
   quasimolion <- which(rules[, "quasi"]== 1) 
  }else{
   quasimolion <- which(rules[, "mandatory"]== 1)
  }
  
  #Remove recognized isotopes from annotation m/z vector
  for(x in seq(along = isotopes)){
    if(!is.null(isotopes[[x]])){
      if(isotopes[[x]]$iso != 0){
        imz[x] <- NA;
      }
    }
  }
  
  #counter for % bar
  npeaks    <- 0; 
  massgrp   <- 0;
  ncl <- sum(sapply(object@pspectra, length));  

  if (runParallel == 1) { ## ... we run in parallel mode
    if(is.null(psg_list)){
      cat('\nCalculating possible adducts in',npspectra,'Groups... \n');
      lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating possible adducts in',length(psg_list),'Groups... \n'); 
      lp <- -1;
      pspectra_list <- psg_list;
    }
    
    argList <- list();
    cnt_peak <- 0;
    if(is.null(max_peaks)){
      max_peaks=100;
    }
    params <- list();
    
    for(j in 1:length(pspectra_list)){
      i <- pspectra_list[j];
      params$i[[length(params$i)+1]] <- i;
      cnt_peak <- cnt_peak+length(object@pspectra[[i]]);
      if(cnt_peak>max_peaks || j == length(pspectra_list)){
        params$pspectra <-object@pspectra;
        params$imz <- imz;
#         params$mint <- mint;
        params$rules <- rules;
        params$mzabs <- mzabs;
        params$devppm <- devppm;
        params$isotopes <- isotopes;
        params$quasimolion <- quasimolion;
        argList[[length(argList)+1]] <- params
        cnt_peak <- 0;
        params <- list();
      }
    }
    #Some informationen for the user
    cat("Parallel mode: There are",length(argList), "tasks.\n")
    
    if(is.null(object@runParallel$cluster)){
      #Use MPI
      result <- xcmsPapply(argList, annotateGrpMPI)
    }else{
      #For snow
      result <- xcms:::xcmsClusterApply(cl=object@runParallel$cluster, 
                                        x=argList, fun=annotateGrpMPI, 
                                        msgfun=msgfun.snowParallel)
    }
    if("typ" %in% colnames(rules)){
      rules.idx <- which(rules[, "typ"]== "A")
      parent <- TRUE;
    }else{
      #backup for old rule sets
      rules.idx <- 1:nrow(rules);
      parent <- FALSE;
    }
    for(ii in 1:length(result)){
      if(length(result[[ii]]) == 0){
        next;
      }
      for(iii in 1:length(result[[ii]])){
        hypothese <- result[[ii]][[iii]];
        if(is.null(hypothese)){
          next;
        }
        charge <- 0;
        old_massgrp <- 0;
        index <- argList[[ii]]$i[[iii]];
        ipeak <- object@pspectra[[index]];
        for(hyp in 1:nrow(hypothese)){
          peakid <- as.numeric(ipeak[hypothese[hyp, "massID"]]);
          if(old_massgrp != hypothese[hyp,"massgrp"]) {
            massgrp <- massgrp+1;
            old_massgrp <- hypothese[hyp,"massgrp"];
            annoGrp <- rbind(annoGrp,c(massgrp,hypothese[hyp,"mass"],sum(hypothese[ which(hypothese[,"massgrp"]==old_massgrp),"score"]),i) ) 
          }
          if(parent){
            annoID <- rbind(annoID, cbind(peakid, massgrp, hypothese[hyp,c("ruleID","parent")]))  
          }else{
            annoID <- rbind(annoID, cbind(peakid, massgrp, hypothese[hyp,c("ruleID")],NA))
          }
          
        }
      }
    }

    derivativeIons <- getderivativeIons(annoID,annoGrp,rules,length(imz));
    cat("\n");
    object@derivativeIons <- derivativeIons;
    object@annoID  <- annoID;
    object@annoGrp <- annoGrp;
    return(object)
  } else {
    ##Single Core Mode
    if(is.null(psg_list)){
      cat('\nCalculating possible adducts in',npspectra,'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating possible adducts in',length(psg_list),'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- psg_list;
      sum_peaks <- sum(sapply(object@pspectra[psg_list],length));
    }
    if("typ" %in% colnames(rules)){
      rules.idx <- which(rules[, "typ"]== "A")
      parent <- TRUE;
    }else{
      #backup for old rule sets
      rules.idx <- 1:nrow(rules);
      parent <- FALSE;
    }
    
    for(j in seq(along = pspectra_list)){
      i <- pspectra_list[j];
      
      #peak index for those in pseudospectrum i
      ipeak <- object@pspectra[[i]];

      #percent output
      npeaks.global <- npeaks.global + length(ipeak);
      perc   <- round((npeaks.global) / ncl * 100)
      perc   <- perc %/% 10 * 10;
      
      if (perc != lp && perc != 0) { 
        cat(perc,' '); 
        lp <- perc;
      }
      if (.Platform$OS.type == "windows"){ 
        flush.console();
      }
      #end percent output

      #check if the pspec contains more than one peak 
      if(length(ipeak) > 1){
        
        hypothese <- annotateGrp(ipeak, imz, rules, mzabs, devppm, isotopes, quasimolion, rules.idx=rules.idx);
        #save results
        if(is.null(hypothese)){
          next;
        }
        charge <- 0;
        old_massgrp <- 0;
        
        #combine annotation hypotheses to annotation groups for one compound mass
        for(hyp in 1:nrow(hypothese)){
          peakid <- as.numeric(ipeak[hypothese[hyp, "massID"]]);
          if(old_massgrp != hypothese[hyp, "massgrp"]) {
            massgrp <- massgrp + 1;
            old_massgrp <- hypothese[hyp, "massgrp"];
            annoGrp <- rbind(annoGrp, c(massgrp, hypothese[hyp, "mass"], 
                                         sum(hypothese[ which(hypothese[, "massgrp"] == old_massgrp), "score"]), i) ) 
          }
          if(parent){
            annoID <- rbind(annoID, c(peakid, massgrp, hypothese[hyp, "ruleID"], ipeak[hypothese[hyp, "parent"]]))
          }else{
            annoID <- rbind(annoID, c(peakid, massgrp, hypothese[hyp, "ruleID"], NA))
          }
        }
      }
    }

    derivativeIons <- getderivativeIons(annoID, annoGrp, rules, length(imz));

    cat("\n");

    object@derivativeIons <- derivativeIons;
    object@annoID  <- annoID;
    object@annoGrp <- annoGrp;
    return(object)
  }
})

annotateDiffreport <- function(object, sample=NA, nSlaves=1, sigma=6, perfwhm=0.6,
  cor_eic_th=0.75, graphMethod="hcs", pval=0.05, calcCiS=TRUE,
  calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5,
  ppm=5, mzabs=0.015, quick=FALSE, psg_list=NULL, rules=NULL,
  polarity="positive", multiplier=3, max_peaks=100, intval="into",
  pval_th = NULL, fc_th = NULL, sortpval=TRUE, ...) {

  if (!class(object)=="xcmsSet") stop ("no xcmsSet object");
  
  #use diffreport from xcms
  diffrep <- diffreport(object, sortpval=FALSE, ...);
  
  if(quick){
    #Quick run, no groupCorr and findAdducts
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa, perfwhm=perfwhm, sigma=sigma);
    xa <- findIsotopes(xa, maxcharge=maxcharge, maxiso=maxiso, ppm=ppm, mzabs=mzabs)
    xa.result <- getPeaklist(xa);
  }else{
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa, perfwhm=perfwhm, sigma=sigma);
    xa <- findIsotopes(xa, maxcharge=maxcharge, maxiso=maxiso, ppm=ppm, mzabs=mzabs)
    if(is.null(psg_list) & is.null(pval_th) & is.null(fc_th)){
      #no restriction for calculation ps-spectra
      #all groups will be calculated, psg_list=NULL
    }else{
      #One value was set
      #generate psg_list
      peaklist <- getPeaklist(xa);
      #Do we have restriction for psg_list?
      if(is.null(psg_list)){
        psg_list <- 1:length(xa@pspectra);
      }
      
      if(!is.null(pval_th)){
        #Find groups, which include features with p-val < pval_th
        index <- which(diffrep[, "pvalue"] < pval_th);
        index.grp <- unique(peaklist[index, "pcgroup"]);
        if(length(index.grp) > 0){
          psg_list  <- psg_list[which(psg_list %in% index.grp)]
        } else {
          cat("No groups found, which satisfy your conditions!\n")
          result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
          return(result);
        }
      }
      
      if(!is.null(fc_th)){
        #Find groups, which include features with fc > fc_th
        index <- which(diffrep[, "fold"] > fc_th);
        index.grp <- unique(peaklist[index, "pcgroup"]);
        if(length(index.grp) > 0){
          psg_list  <- psg_list[which(psg_list %in% index.grp)]
        } else {
          cat("No groups found, which satisfy your conditions!\n")
          result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
          return(result);
        }
      }
      
      if(length(psg_list) < 1){
        #no group satisfy conditions
        cat("No groups found, which satisfy your conditions!\n")
        result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
        return(result);
      }
    }
    
    #Include into psg_list all groups that has been created after groupCorr
    cnt <- length(xa@pspectra);
    xa <- groupCorr(xa,cor_eic_th=cor_eic_th,psg_list=psg_list)
    if(!is.null(psg_list)){
      psg_list <- c(psg_list,(cnt+1):length(xa@pspectra));
    }
    xa <- findAdducts(xa, multiplier=multiplier, ppm=ppm, mzabs=mzabs, 
                      polarity=polarity, rules=rules, psg_list=psg_list);
    
    xa.result <- getPeaklist(xa);
  }
  
  #combines results
  result <- cbind(diffrep, xa.result[, c("isotopes", "adduct", "pcgroup")])
  if(sortpval){
    #return diffreport order if p-value sorting is true
    result <- result[order(result[, "pvalue"]), ];
  }
  return(result);
} 

###End xsAnnotate generic Methods###

###xsAnnotate exported Methods###

#getpspectra
#Generate for one (or more) pseudospectrum a peaklist
#object - xsAnnotate
#grp - pseudospectrum ID
getpspectra <- function(object, grp=NULL){
  
  if(is.null(grp)) {
    cat("Error: No grp number!\n");
  } else if(!(is.numeric(grp)) | any(grp > length(object@pspectra))) {
    cat("Error: Grp is not numeric or contains number greater than maximum number of pseudospectra!\n");
  } else {
    index <- unlist(object@pspectra[grp]);
    #Check if index has peaks
    if(length(index) == 0){
      cat("Error: Pseudospectra selection contains no peaks.\n")
      return(NULL);
    }
    #get ps number for selected peaks  
    grpvec <- unlist(sapply(grp, function(x) { rep(x, length(object@pspectra[[x]])) }))
    
    #extract peaktable
    peaktable <- object@groupInfo[index,]

    adduct    <- vector("character",length(index));
    isotopes  <- vector("character",length(index));
    ions      <- object@derivativeIons;
    iso       <- object@isotopes;

    lions <- ions[index];
    liso  <- iso[index];

    #default polarity set to positive
    polarity <- "+";
    
    if(length(object@polarity) > 0){
      if(object@polarity == "negative"){
        polarity <- "-";
      }
    }
    
    for(i in 1:length(lions)){
      if(!is.null(lions[[i]])){
        if(length(lions[[i]]) > 1){
          names <- c();
          for(ii in 1:length(lions[[i]])){
            names <- paste(names,lions[[i]][[ii]]$name,lions[[i]][[ii]]$mass);
          }
          adduct[i] <- names;
        } else {
          adduct[i] <- paste(lions[[i]][[1]]$name,lions[[i]][[1]]$mass);
        }
      }
        
      if(!is.null(liso[[i]])){
        if(liso[[i]]$iso == 0){
          iso.name <- "[M]";
        }else{
          iso.name <- paste("[M+",liso[[i]]$iso,"]",sep="");
        }
        if(liso[[i]]$charge > 1){
          isotopes[i] <- paste("[",liso[[i]]$y,"] ",iso.name," ",liso[[i]]$charge,polarity,sep="");
        }else{
          isotopes[i] <- paste("[",liso[[i]]$y,"] ",iso.name," ",polarity,sep="");
        }
      }
    }

    #Check if peaktable has only one row
    if(is.null(nrow(peaktable))){ 
      peaktable <- matrix(peaktable, byrow=F, ncol=length(peaktable));
    }
    colnames(peaktable) <- colnames(object@groupInfo)
    
    return(invisible(data.frame(peaktable,isotopes,adduct,psg=grpvec,stringsAsFactors=FALSE)));
  }
}


setGeneric("getPeaklist", function(object, intval="into") standardGeneric("getPeaklist"))
setMethod("getPeaklist", "xsAnnotate", function(object, intval="into") {
  
  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!")
  }

  #generate peaktable
  #Check if xcmsSet contains only one sample
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    #intval is here ignored since all intensity values are already contained
    peaktable <- object@groupInfo;
  }else {
    #Case of xcmsSet with multiple samples
    #Use groupInfo information and replace intensity values
    peaktable <- object@groupInfo;
    
    #get intensity values from xcmsSet
    grpval <- groupval(object@xcmsSet, value=intval);
    
    #get column range for replacement
    grpval.ncol <- ncol(grpval)
    start <- ncol(peaktable) - grpval.ncol +1;
    ende  <- start + grpval.ncol - 1; 
    
    peaktable[, start:ende] <- grpval;
  }

  #allocate variables for CAMERA output
  adduct   <- vector("character", nrow(object@groupInfo));
  isotopes <- vector("character", nrow(object@groupInfo));
  pcgroup  <- vector("character", nrow(object@groupInfo));
   
  #default polarity set to positive
  polarity <- "+";
    
  if(length(object@polarity) > 0){
    if(object@polarity == "negative"){
      polarity <- "-";
    }
  }
  
  #First isotope informationen and adduct informationen
  for(i in seq(along = isotopes)){
    #check if adduct annotation is present for peak i
    if(length(object@derivativeIons) > 0 && !(is.null(object@derivativeIons[[i]]))) {
      #Check if we have more than one annotation for peak i
      if(length(object@derivativeIons[[i]]) > 1) {
        #combine ion species name and rounded mass hypophysis
        names <- paste(object@derivativeIons[[i]][[1]]$name, signif(object@derivativeIons[[i]][[1]]$mass, 6));
        for(ii in 2:length(object@derivativeIons[[i]])) {
          names <- paste(names, object@derivativeIons[[i]][[ii]]$name, signif(object@derivativeIons[[i]][[ii]]$mass, 6));
        }
        #save name in vector adduct
        adduct[i] <- names;
      } else {
        #Only one annotation
        adduct[i] <- paste(object@derivativeIons[[i]][[1]]$name, signif(object@derivativeIons[[i]][[1]]$mass, 6));
      }
    } else {
      #no annotation empty name
      adduct[i] <- ""; 
    }
    
    #Check if we have isotope informationen about peak i
    if(length(object@isotopes) > 0&& !is.null(object@isotopes[[i]])) {
      num.iso <- object@isotopes[[i]]$iso;
      #Which isotope peak is peak i?
      if(num.iso == 0){
        str.iso <- "[M]";
      } else { 
        str.iso <- paste("[M+", num.iso, "]", sep="")
      }
      #Multiple charged?
      if(object@isotopes[[i]]$charge > 1){
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, object@isotopes[[i]]$charge, polarity, sep="");
      }else{
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, polarity, sep="");
      }
    } else { 
      #No isotope informationen available
      isotopes[i] <- ""; 
    }
  }
  
  #Have we more than one pseudospectrum?
  if(length(object@pspectra) < 1){
      pcgroup <- 0;
  } else {
    for(i in seq(along = object@pspectra)){
      index <- object@pspectra[[i]];
      pcgroup[index] <- i;
    }
  }
          
  rownames(peaktable)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
  return(invisible(data.frame(peaktable, isotopes, adduct, pcgroup, stringsAsFactors=FALSE, row.names=NULL)));
})

setGeneric("annotate", function(object, sample=NA, nSlaves=1, sigma=6, perfwhm=0.6, cor_eic_th=0.75, graphMethod="hcs",
  pval=0.05, calcCiS=TRUE, calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5, ppm=5, mzabs=0.015, 
  quick=FALSE, psg_list=NULL,  rules=NULL, polarity="positive", multiplier=3, max_peaks=100 ,intval="into") standardGeneric("annotate"))

setMethod("annotate", "xcmsSet", function(object, sample=NA, nSlaves=1, sigma=6, perfwhm=0.6, cor_eic_th=0.75, graphMethod="hcs",
  pval=0.05, calcCiS=TRUE, calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5, ppm=5, mzabs=0.015, 
  quick=FALSE, psg_list=NULL,  rules=NULL, polarity="positive", multiplier=3, max_peaks=100 ,intval="into") {

  #check intval
  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!\n")
  }

  #check graphMethod
  if(!sum(graphMethod == c("hcs","lpc"))){
      stop("Unknown graphMethod value!\n")
  }

  if(quick){
    #Quick run, no groupCorr and findAdducts
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa, perfwhm=perfwhm, sigma=sigma, intval=intval);
    xa <- findIsotopes(xa, maxcharge=maxcharge, maxiso=maxiso, ppm=ppm, mzabs=mzabs, intval=intval, minfrac=minfrac)

  }else{

    xa  <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa  <- groupFWHM(xa, perfwhm=perfwhm, sigma=sigma, intval=intval);
    xa  <- findIsotopes(xa, maxcharge=maxcharge, maxiso=maxiso, ppm=ppm, mzabs=mzabs, intval=intval, minfrac=minfrac)
    cnt <- length(xa@pspectra);
    xa  <- groupCorr(xa, cor_eic_th=cor_eic_th, graphMethod=graphMethod, calcIso=calcIso, calcCiS=calcCiS, calcCaS=calcCaS, psg_list=psg_list)
    if(!is.null(psg_list)){
      psg_list <- c(psg_list,(cnt+1):length(xa@pspectra));
    }
    xa <- findAdducts(xa, multiplier=multiplier, ppm=ppm, rules=rules, max_peaks=max_peaks, mzabs=mzabs, polarity=polarity, psg_list=psg_list);
  }
  #Kombiniere Resultate

  return(xa);
})

findKentrickMasses <- function(object , masses=c(14,14.01565), nHomologue=4, 
                               error=0.002, filter=TRUE) {
  ##Screen for mass differences
########Variables########
filename.Input  <- ""
filename.Output <- ""

#Mass Difference in Da
#nominal Mass (14 for CH2)
nominalMass <- 14 
exactMass   <- 14.01565

#max differenze is maxHomologue x nominalMass
maxHomologue <- 4

#Error in Da for matching Kendrick mass defect
error <- 0.002 

#filter rt, higher mass must have higher retetion time
filter <- TRUE
######End Variables######

#Read data file
data <- read.delim(filename.Input, sep=",")

if(any(sapply(data[, 1], is.na))){
  stop("No NA values allowed!")
}

#sorting after mass(first column)
data.sorted <- data[order(data[, 1]), ]

kendrickMass <- data.sorted[, 1] * nominalMass/exactMass;
kendrickMassDefect <- ceiling(kendrickMass) - kendrickMass
kendrickMass <- ceiling(kendrickMass)

#generate index vector
allCandidates <- 1:length(kendrickMassDefect);

#generate result list;
results <- list();

# while loop
while(length(allCandidates) > 0){
  difference <- kendrickMass[allCandidates] - kendrickMass[allCandidates[1]]
  index <- which(difference %% nominalMass == 0 & difference > 0 & 
    difference <= maxHomologue*nominalMass)
  
  if(length(index) < 1){
    allCandidates <- allCandidates[-1, drop=FALSE]
    next;
  }
  
  if(filter){
    index <- index[which(data.sorted[allCandidates[index], 2] >= data.sorted[allCandidates[1], 2])]
    if(length(index) < 1){
      allCandidates <- allCandidates[-1, drop=FALSE]
      next;
    }
  }
  
  index.Kendrick <- CAMERA:::fastMatch(kendrickMassDefect[allCandidates[index]],
                                       kendrickMassDefect[allCandidates[1]], 
                                       tol=error)
  if(any(hit <- !sapply(index.Kendrick, is.null))){
    hits <- c(1, index[which(hit)])
    results[[length(results)+1]] <- allCandidates[hits];  
  }
    allCandidates <- allCandidates[-1, drop=FALSE]
}


resultMatrix <- matrix(NA, nrow=0, ncol=5)

#generate Results
invisible(lapply(results, function(x){
   resultMatrix <<- rbind(resultMatrix, cbind(data.sorted[x, ], kendrickMass[x],
                                              kendrickMassDefect[x]), rep(NA, 5))
}))
  
write.csv(resultMatrix, filename.Output, 
          row.names=FALSE, na="")

}
#Find NeutralLossSpecs
#Test every pseudospectrum in a xsAnnotate object if the peaks match
#a specific mz difference. Returns a boolean vector, where a hit is marked
#with TRUE.
#object - xsAnnotate object
#mzdiff - mz difference of interest
#mzabs  - allowed absolute error in mz
#mppam  - allowed relative error in ppm (not used)

findNeutralLossSpecs <- function(object, mzdiff=NULL, mzabs=0, mzppm=10) {

  if (!class(object) == "xsAnnotate"){ 
    stop ("Parameter object is no xsAnnotate object\n")
  }
  
  if(is.null(mzdiff)){
    stop ("Parameter mzdiff must be set\n")
  }
  
  if(!is.numeric(mzdiff) || any(mzdiff <= 0)) {
    stop ("Parameter mzdiff must be a positive numeric value\n")
  }
  
  if(!is.numeric(mzabs) || mzabs < 0) {
    stop ("Parameter mzabs must be a positive numeric value\n")
  }
  
  if(!is.numeric(mzppm) || mzppm < 0) {
    stop ("Parameter mzppm must be a positive numeric value\n")
  }
  
  # get mz values from peaklist
  imz    <- object@groupInfo[, "mz"];
  
  #get Isotopes
  isotopes  <- object@isotopes;
  #Remove recognized isotopes from annotation m/z vector
  for(x in seq(along = isotopes)){
    if(!is.null(isotopes[[x]])){
      if(isotopes[[x]]$iso != 0){
        imz[x] <- NA;
      }
    }
  }
  
  #Calculate mzrange
  nlmin <- mzdiff-mzabs
  nlmax <- mzdiff+mzabs

  
  hits <- sapply(1:length(object@pspectra), function(j) {
    spec  <-  object@pspectra[[j]]

    mz  <-  CAMERA:::naOmit(imz[spec]) #mz vector
    nl  <-  as.matrix(dist(mz, method="manhattan")) #distance matrix
    hit <- FALSE;
    for(x in seq_along(mz)){
      mzadd <- mz[x]*mzppm*10^-6
      for(y in seq_along(nlmin)){
        if(any(which(nl[x,] > (nlmin[y]-mzadd) & nl[x,] < (nlmax[y] + mzadd)))){
          hit <- TRUE;
          break;
        }
      }
    }
    hit
  })
  return(hits)
}

#Find NeutralLosses
#Test every pseudospectrum in a xsAnnotate object if two peaks match
#a specific mz difference. Returns a artifical xcmsSet
#object - xsAnnotate object
#mzdiff - mz difference of interest
#mzabs  - allowed absolute error in mz
#mppam  - allowed relative error in ppm (not used)
findNeutralLoss <- function(object, mzdiff=NULL, mzabs=0, mzppm=10) {

  if (!class(object) == "xsAnnotate"){ 
    stop ("Parameter object is no xsAnnotate object\n")
  }
  
  if(is.null(mzdiff)){
    stop ("Parameter mzdiff must be set\n")
  }
  
  if(!is.numeric(mzdiff) || any(mzdiff <= 0)) {
    stop ("Parameter mzdiff must be a positive numeric value\n")
  }
  
  if(!is.numeric(mzabs) || mzabs < 0) {
    stop ("Parameter mzabs must be a positive numeric value\n")
  }
  
  if(!is.numeric(mzppm) || mzppm < 0) {
    stop ("Parameter mzppm must be a positive numeric value\n")
  }
  
  #Calculate mzrange
  nlmin <- mzdiff-mzabs
  nlmax <- mzdiff+mzabs

  #Create tempory xcmsSet object to store the peak matches
  xs <- new("xcmsSet");

  peaks <- matrix(ncol=ncol(object@groupInfo)+1, nrow=0);
 
  sampnames(xs) <- c("NeutralLoss")
  sampclass(xs) <- c("NeutralLoss")

  # get mz values from peaklist
  imz    <- object@groupInfo[, "mz"];
  
  #get Isotopes
  isotopes  <- object@isotopes;
  #Remove recognized isotopes from annotation m/z vector
  for(x in seq(along = isotopes)){
    if(!is.null(isotopes[[x]])){
      if(isotopes[[x]]$iso != 0){
        imz[x] <- NA;
      }
    }
  }
  
  #Loop over all pseudospectra
  for(j in seq(along = object@pspectra)) {
    
    #get specific pseudospectrum j
    spec  <-  object@pspectra[[j]]
    #lastcol <- ncol(spec);
    
    #get mz values
    mz  <-  imz[spec] 
    
    #generate distance matrix
    nl  <-  as.matrix(dist(mz, method="manhattan"))
    
    hits <- c();
    for(x in seq_along(mz)){
      mzadd <- mz[x]*mzppm*10^-6
      if(is.na(mzadd)){
        next;
      }
      for(y in seq_along(nlmin)){
        if(length( index <- which(nl[x,] > (nlmin[y]-mzadd) & 
          nl[x,] < (nlmax[y] + mzadd))) > 0 ){
          #found hit
          hits <- rbind(hits,cbind(x,index))
        }
      }
    }
    if(is.null(hits)){
      next;
    }
    ## Remove hit from lower diagonal matrix
    index <- which(apply(hits,1, function(x) x[1] > x[2]))
    if(length(index) > 0){
      hits <- hits[-index,,drop=FALSE]
    }
    # Remove duplicated
    hits <- unique(hits)
    # loop over all hits
    for (f in seq_len(nrow(hits))) {
      # hit f
      hit <- hits[f,]
        
      

      #For each peak-pair save the greater peak(higher mz value) to peak table
      if (mz[hit[1]] > mz[hit[2]]) {
        newpeak <- object@groupInfo[spec[hit[1]],, drop=FALSE]
      } else {
        newpeak <- object@groupInfo[spec[hit[2]],, drop=FALSE]
      }  
      peaks <- rbind(peaks, cbind(newpeak,j))
    }
  }
  colnames(peaks) <- c(colnames(object@groupInfo),"compound spectrum");
  #Store peaks into the xcmsSet
  peaks(xs) <- as.matrix(peaks)
  invisible(xs)
}

msgfun.snowParallel <- function(x,i){
  cat("Sending task # ",i,"\n");
  flush.console();
}

cleanParallel <- function(object){
  ##testing objects
  if (!class(object) == "xsAnnotate"){
    stop ("xsa.pos is no xsAnnotate object")
  }
  
  if(object@runParallel$enable == 1){
    if(object@runParallel$mode == "Rmpi"){
      mpi.close.Rslaves()
    } else if(object@runParallel$mode == "snow" & !(is.null(object@runParallel$cluster))){
      stopCluster(object@runParallel$cluster);
    }
    cat("Slaves were stopped!\n");
  }
}
  
###End xsAnnotate exported Methods###


###xsAnnotate intern Function###

combinexsAnnos <- function(xsa.pos, xsa.neg, pos=TRUE, tol=2, ruleset=NULL){
  # two xsAnnotate objects (pos,neg)
  # pos: returns pos. (true) or neg. (false) peaklist;  default = TRUE
  # tol: max allowed time difference between pos and neg pseudospectra
  
  ##testing objects
  if (!class(xsa.pos) == "xsAnnotate"){
    stop ("xsa.pos is no xsAnnotate object")
  }
  if (!class(xsa.neg) == "xsAnnotate"){ 
    stop ("xsa.neg is no xsAnnotate object")
  }
  if (xsa.pos@polarity != "positive" & xsa.neg@polarity != "negative"){
    stop ("xsAnnotate object have unknown polarities.\nOnly positive or negative are allowed!")
  }
  
  ##1. Step
  #get all rts for every pseudospectra
  rt1 <- sapply(xsa.pos@pspectra, function(x) { 
          mean(xsa.pos@groupInfo[x, "rt"]) 
    })
  
  rt2 <- sapply(xsa.neg@pspectra, function(x) { 
          mean(xsa.neg@groupInfo[x, "rt"]) 
    })
  
  #find matching pseudospectra
  ps.match <- CAMERA:::fastMatch(rt1, rt2, tol=tol);
  rules.pos <- xsa.pos@ruleset;
  rules.neg <- xsa.neg@ruleset;
  
  #ruleset
  if(is.null(ruleset)){
    #generate M+H,M-H rule
    if(pos){
      name <- "M+H/M-H";
    }else{
      name <- "M-H/M+H";
    }
    ruleset <- data.frame(name, 1, 1, 1, 1, 2.014552, 1, 1); 
    colnames(ruleset) <- c("name", "nmol.pos", "nmol.neg", "charge.pos", "charge.neg", 
                           "massdiff", "ID.pos", "ID.neg");
  }else{
    #generate rules stated in ruleset
    if(!is.matrix(ruleset) || ncol(ruleset) != 2){
      stop("Ruleset is no matrix or number of cols is unequal two");
    }else{
      #check provided ruleset
      name <- c(); 
      nmol.pos <- c(); 
      nmol.neg <- c();
      charge.pos <- c(); 
      charge.neg <- c();
      massdiff   <- c();
         apply(ruleset,1, function(x) {
            name <<- c(name,paste(rules.pos[x[1],"name"],rules.neg[x[2],"name"],sep="/"));
            nmol.pos <<- c(nmol.pos,rules.pos[x[1],"nmol"]);
            nmol.neg <<- c(nmol.neg,rules.neg[x[2],"nmol"]);  
            charge.pos <<- c(charge.pos,rules.pos[x[1],"charge"]);
            charge.neg <<- c(charge.neg,rules.neg[x[2],"charge"]);
            massdiff <<- c(massdiff, rules.pos[x[1],"massdiff"] - rules.neg[x[2],"massdiff"]);
            return();
        })
      #check charges
      if(! all(charge.pos == abs(charge.neg))){
        stop("Ruleset can't combine rules with unequal charges")
      }
      ruleset <- data.frame(name,nmol.pos,nmol.neg,charge.pos,charge.neg,massdiff,ruleset);
      colnames(ruleset) <- c("name", "nmol.pos", "nmol.neg","charge.pos","charge.neg", "massdiff", "ID.pos","ID.neg");
    }
  }

  #allocate variable for matching results 
  endresult <- matrix(ncol=6, nrow=0);
  colnames(endresult) <- c("grp.pos", "peak.pos", "grp.neg", "peak.neg", "mass", "rule")
   
  #Annotations from positive Samples
  annoID.pos  <- xsa.pos@annoID; 
  annoGrp.pos <- xsa.pos@annoGrp;
  #Annotations from negative Samples
  annoID.neg  <- xsa.neg@annoID;
  annoGrp.neg <- xsa.neg@annoGrp;
  
  #Groups which can be deleted
  grp2del  <- c();
  grp2save <- c();
  peaklist <- c(); #set global variable

  cat ("Run combining of xsAnnotates\n % finished: ");
  lp <- -1;
  ncl<-sum(sapply(ps.match,length));
  npeaks.global <- 0;
  
  #for every pseudospectra
  for(i in seq(along=ps.match)){
    percentOutput(npeaks.global, length(ps.match[[i]]), ncl, lp)
    #check for match
    if(is.null(ps.match[[i]])){
      #no corresponding match in neg sample
      next;
    }

    #matrix stores matches
    result.matrix <- matrix(ncol=4, nrow=0);
    colnames(result.matrix) <- c("pos.peak.ID", "neg.peak.ID", "j", "rule");
    
    #get all m/z values from pos pseudospectrum
    grp.pos <- getpspectra(xsa.pos, i);
    m.pos <- grp.pos[, "mz"]

    #remove m/z values of annotated isotope peaks
    for(k in 1:nrow(grp.pos)){
      if(!is.null(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]])){
        if(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]]$iso > 0){
          m.pos[k] <- NA;
        }
      }
    }
    
    na.ini.pos <- which(!is.na(m.pos)); #index of non NA values
    
    if(length(na.ini.pos) < 1){
      next;
    }
    #get all m/z hypotheses (if exists)
    if(length(index <- which(annoGrp.pos[, 4] == i)) > 0){
       masslist.pos <- annoGrp.pos[index, 2];
    }else{
       #no annotation exists
       masslist.pos <- NULL;
    }
    
    #for every matching neg pseudospectra
    for(j in seq(along=ps.match[[i]])){
      #get all m/z values from neg pseudospectrum
      grp.neg <- getpspectra(xsa.neg, ps.match[[i]][j]);
      m.neg   <- grp.neg[,"mz"]
      #remove masses of annotated isotope peaks
      for(k in 1:nrow(grp.neg)){
        if(!is.null(xsa.neg@isotopes[[xsa.neg@pspectra[[ps.match[[i]][j]]][k]]])){
          if(xsa.neg@isotopes[[xsa.neg@pspectra[[ps.match[[i]][j]]][k]]]$iso > 0){
            m.neg[k] <- NA;
            next;
          }
        }
      }

      na.ini.neg <- which(!is.na(m.neg));#index of non NA values
      if(length(na.ini.neg) < 1){
        next;
      }
      #get all m/z hypotheses (if exists)
      if(length(index <- which(annoGrp.pos[, 4] == ps.match[[i]][j])) > 0){
        masslist.neg <- annoGrp.pos[index, 2];
      }else{
        #no annotation exists
        masslist.neg <- NULL;
      }

      #match rules against m/z values
      results <- matrix(NA, nrow=length(m.pos), ncol=nrow(ruleset));
      results[na.ini.pos, ] <- CAMERA:::combineHypothese(naOmit(m.pos), naOmit(m.neg), 
                                                ruleset, na.ini.pos, na.ini.neg);
      
      if(!all(is.na(results))){
        #Found Hits between pos and neg datasets with ruleset
        for(l in 1:ncol(results)){
          index <- which(!is.na(results[, l]));
          if(length(index) > 0){
            #peak index of matching peaks
            #first column pos, second neg
            tmp <- matrix(c(xsa.pos@pspectra[[i]][index], 
                            xsa.neg@pspectra[[ps.match[[i]][j]]][results[index, l]]),
                          ncol=2);
            result.matrix <- rbind(result.matrix, cbind(tmp,ps.match[[i]][j],l));
          }
        }
      }
    }
    
    if(! nrow(result.matrix) > 0){
      next; #Found no hit
    }

    for(ii in 1:nrow(result.matrix)){
      if(length(index <- which(annoID.pos[, "id"] == result.matrix[ii, 1] )) > 0){
        #Peak has annotation(s)
        if(all(annoID.pos[index, "ruleID"] == ruleset[result.matrix[ii,4], "ID.pos"])){
          #Peak has only one annotation and could be verified
          mass <- 2;
          endresult <- rbind(endresult, c(i, result.matrix[ii, 1], result.matrix[ii, 3],
                                          result.matrix[ii, 2], mass, result.matrix[ii, 4]));
          grp2save  <- c(grp2save, annoID.pos[index, "grpID"]);
        }else{
          #Peak has more than one annotation or verfication goes wrong
          grp.save <- which(annoID.pos[index, "ruleID"] == ruleset[result.matrix[ii, 4], "ID.pos"]);
          if(length(grp.save) > 0 ){
            #Save verified annotation and remove from index
            grp2save  <- c(grp2save, annoID.pos[index[grp.save], "grpID"]);
            index <- index[-grp2save];
            mass <- 2;
            endresult <- rbind(endresult, c(i, result.matrix[ii, 1], result.matrix[ii, 3],
                                            result.matrix[ii, 2], mass,result.matrix[ii, 4]));  
          }else{
            #Found new annotation
            mass <- 1;
            endresult <- rbind(endresult, c(i, result.matrix[ii, 1], result.matrix[ii, 3],
                                            result.matrix[ii, 2], mass, result.matrix[ii, 4]));
          }  
          #delete all other hypotheses
          grp2del  <- c(grp2del, annoID.pos[index, "grpID"]);
        }
      }else{
        #Peak has no annotation
        #Add annotation according to ruleset
        mass <- 1;
        endresult <- rbind(endresult, c(i, result.matrix[ii, 1], result.matrix[ii, 3],
                                        result.matrix[ii, 2], mass, result.matrix[ii, 4]));
      }
    }
  }#end for i loop

  #Remove grp2del groups, if they are in grp2save
  grp2del <- unique(grp2del);
  grp2save <- unique(grp2save);
  if(length(grp2del) > 0 & length(grp2save) > 0){
    index <- which(grp2del %in% grp2save);
    if(length(index) > 0){
      grp2del <- grp2del[-index];
    }
  }
  cat("\nGenerate Peaklist....\n");
  if(pos){
    #return postive peaklist
    if(length(grp2del) > 0){
      index <- which(annoID.pos[,2] %in% grp2del);
      #annotation to delete
      if(length(index) > 0){
        annoID.pos <- annoID.pos[-index,];
      }
    }

    #new vector to peaklist
    add.adducts <- vector("character", nrow(xsa.pos@groupInfo));
    old.grpid   <- max(annoGrp.pos[, 1]); #counter for grp-ID

    if(nrow(endresult) > 0){
      for(i in 1:nrow(endresult)){
        if(endresult[i, 5] == 1){
          old.grpid <- old.grpid + 1;
          rule.ID <- ruleset[endresult[i,6], "ID.pos"];
          annoID.pos  <- rbind(annoID.pos, c(endresult[i, "peak.pos"], old.grpid, rule.ID, NA))
          mass <- (xsa.pos@groupInfo[endresult[i,"peak.pos"],"mz"]*rules.pos[rule.ID,"charge"] - rules.pos[rule.ID,"massdiff"]) / rules.pos[rule.ID,"nmol"];
          annoGrp.pos <- rbind(annoGrp.pos, c(old.grpid,mass,2,endresult[i,1]))
        }
        add.adducts[endresult[i,"peak.pos"]]<-paste("Found",ruleset[endresult[i,6],1]);
      }
    }
    #save results
    xsa.pos@derivativeIons <- getderivativeIons(annoID.pos, annoGrp.pos, rules.pos, 
                                                length(xsa.pos@isotopes));
    peaklist <- getPeaklist(xsa.pos);
    index    <- ncol(peaklist);
    peaklist<-cbind(peaklist,add.adducts);
    colnames(peaklist)[index+1] <- "neg. Mode"
#     return(peaklist);

  }else{
    #return negative peaklist
    if(length(grp2del) > 0){
      index <- which(annoID.neg[,2] %in% grp2del);
      #annotation to delete
      if(length(index) > 0){
        annoID.neg <- annoID.neg[-index,];
      }
    }

    #new vector to peaklist
    add.adducts <- vector("character", nrow(xsa.neg@groupInfo));
    old.grpid   <- max(annoGrp.neg[, 1]);

    if(nrow(endresult) > 0){
      for(i in 1:nrow(endresult)){
        if(endresult[i, 5] == 1){
          old.grpid <- old.grpid + 1;
          rule.ID <- ruleset[endresult[i,6],"ID.neg"];
          annoID.neg  <- rbind(annoID.neg,  c(endresult[i,"peak.neg"],old.grpid,rule.ID,NA))
          mass <- (xsa.neg@groupInfo[endresult[i,"peak.neg"],"mz"]* abs(rules.neg[rule.ID,"charge"]) - rules.neg[rule.ID,"massdiff"]) / rules.neg[rule.ID,"nmol"];
          annoGrp.neg <- rbind(annoGrp.neg, c(old.grpid,mass,2,endresult[i,1]))
        }
        add.adducts[endresult[i,"peak.neg"]]<-paste("Found",ruleset[endresult[i,6],1]);
      }
    }

    xsa.neg@derivativeIons <- getderivativeIons(annoID.neg,annoGrp.neg,rules.neg,length(xsa.neg@isotopes)); #save results
    peaklist<-getPeaklist(xsa.neg);
    index<-ncol(peaklist)
    peaklist<-cbind(peaklist,add.adducts);
    colnames(peaklist)[index+1] <- "pos. Mode"
  }

  return(peaklist);
}    


combineHypothese <- function(mass.pos, mass.neg, ruleset, ini1, ini2, tol=0.02){
    ##generiere neue Peaklist
    tmp.ruleset <- ruleset[,c("name","nmol.pos","charge.pos","massdiff")];
    colnames(tmp.ruleset) <- c("name","nmol","charge","massdiff")
    ML <- massDiffMatrix(mass.pos,tmp.ruleset)
    ML.match <- apply(ML,2,function(x) {
      m <- fastMatch(x,mass.neg,tol=tol);
      unlist(lapply(m, function(x) { if(is.null(x)){return(NA)}else{return(ini2[x[1]])} })); ##Select first hit, if more than one??
    });
    return(ML.match);
}

##Speed-up function for ignoring NA values
naOmit <- function(x) {
  return (x[!is.na(x)]);
}

# Create complete feature table
# xs - xcmsSet object
# method - groupval parameter method
# value - groupval parameter method
getPeaks_selection <- function(xs, method="medret", value="into"){
  
  if (!class(xs) == "xcmsSet") {
    stop ("Parameter xs is no xcmsSet object\n")
  }
  
  # Testing if xcmsSet is grouped
  if (nrow(xs@groups) > 0 && length(xs@filepaths) > 1) {
    # get grouping information
     groupmat <- groups(xs)
     # generate data.frame for peaktable
     ts <- data.frame(cbind(groupmat, groupval(xs, method=method, value=value)), row.names = NULL)
     #rename column names
     cnames <- colnames(ts)
     if (cnames[1] == 'mzmed') {
       cnames[1] <- 'mz' 
     } else { 
       stop ('Peak information ?!?')
     }
     if (cnames[4] == 'rtmed') {
       cnames[4] <- 'rt' 
     } else {
       stop ('Peak information ?!?')
     }
     colnames(ts) <- cnames
  } else if (length(sampnames(xs)) == 1) { #Contains only one sample?
    ts <- xs@peaks
  } else {
    stop ('First argument must be a xcmsSet with group information or contain only one sample.')
  }
  
  return(as.matrix(ts))
}


# Returns peak table of one! sample from a xcmsSet
# Does not return the feature table (result from group)
# xs - xcmsSet object
# index - sample index
# index = -1 return complete peak table overall sample
getPeaks <- function(xs, index=1){
 
  if (!class(xs) == "xcmsSet") {
    stop ("Parameter xs is no xcmsSet object\n")
  }
  if (!is.numeric(index) && index > length(sampnames(xs)) | index < -1 | index == 0) {
    stop("Parameter index must be between 1 and number of samples or -1\n")
  }
  
  #Testing if xcmsSet is grouped
  if (nrow(xs@groups) > 0) {
    #Should all peaks returned
    if(index == -1) {
      ts <- xs@peaks;
    } else {
      #get peak indices for sample index
      peaki <- getPeaksIdxCol(xs,NULL)[,index]
      #extract peaks from sample index
      ts <- xs@peaks[peaki,]
    }
  } else if (length(sampnames(xs)) == 1) {# xs not grouped, testing if it only contains one sample
        ts <- xs@peaks
  } else {
    stop ('First argument must be a xcmsSet with group information or contain only one sample.')
  }
  
  return(as.matrix(ts))
}


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#function write percent Output
percentOutput <- function(len.old, len.add, len.max, cnt){
  len <- len.old + len.add;
  perc <- round((len) / len.max * 100)
  perc <- perc %/% 10 * 10;
  if (perc != cnt && perc != 0) { 
    cat(perc,' '); 
    cnt2 <- perc;
    eval.parent(substitute(cnt <- cnt2))
  }
  if (.Platform$OS.type == "windows"){ 
    flush.console();
  }
  eval.parent(substitute(len.old <- len))
}
###End xsAnnotate intern Function###
