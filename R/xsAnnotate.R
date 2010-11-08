xsAnnotate <- function(xs=NULL, sample=NA, nSlaves=1){

 ## sample is:
 ### NA for maxInt-way
 ##  1-nSamp for manual way
 ##  -1 for specDist way (to be implemented)
  if(is.null(xs)) { stop("no argument was given"); }
  else if(!class(xs)=="xcmsSet") stop("xs is no xcmsSet object ");
  
  object  <- new("xsAnnotate");
  
  if(length(sampnames(xs)) > 1){  ## more than one sample
    if (!nrow(xs@groups) > 0) {
      stop ('First argument must be a xcmsSet with group information or contain only one sample.') #checking alignment
    }
    if(is.null(sample) || is.na(sample)) {
      object@sample<-as.numeric(NA);
    }else{
      if(sample == -1){
        #Joes Way
        object@sample=sample;
      }else if(length(xs@filepaths) < sample | sample < 1) {
        stop("Paramete sample must be lower equal than number of samples and greater than 0.\n")
      }else{
        object@sample<-sample;
      }
    }
    object@groupInfo <- getPeaks_selection(xs);
  } else if(length(sampnames(xs)) == 1){ ##only one sample was machen wir denn hiermit? wo ist denn das sample=-1 wichtig?
      object@sample= 1;
      object@groupInfo <- getPeaks_selection(xs)
#       object@groupVal <- matrix(ncol=1,nrow=nrow(peaks(xs)),data=seq(1:nrow(peaks(xs))));
#       rownames(object@groupVal) = sapply(1:nrow(peaks(xs)),function(x){paste(round(peaks(xs)[x,"mz"],2),"/",round(peaks(xs)[x,"rt"]),sep="")});
  }else { stop("Unknown error with a grouped xcmsSet"); }

  runParallel<-0;
  if (nSlaves > 1) {
#     cat("Parallel mode is currently not avaible.\nWill be re-enabled with the next CAMERA version!\n");
    ## If MPI is available ...
    rmpi = "Rmpi"
    if (require(rmpi,character.only=TRUE) && !is.null(nSlaves)) {
      if (is.loaded('mpi_initialize')) {
        #test if not already slaves are running!
        if(mpi.comm.size() >0){ 
          warning("There are already intialized mpi slaves on your machine.\nCamera will try to uses them!\n");
          runParallel<-1;
        }else{
          mpi.spawn.Rslaves(nslaves=nSlaves, needlog=FALSE)
          if(mpi.comm.size() > 1){
            #Slaves have successfull spawned
            runParallel<-1;
          }else{ warning("Spawning of mpi slaves have failed. CAMERA will run without parallelization.\n");}
        }
      }else {
          #And now??
      }
    }
  }

  object@runParallel<-runParallel;

  colnames(object@annoID) <-  c("id","grp_id","rule_id");
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
  cat("With",length(sampnames(object@xcmsSet)),"samples and",nrow(object@groupInfo),"peaks\n");
 
  if(is.na(object@sample)){
    cat(paste("Using automatic sample selection\n"));
  } else if(object@sample>-1){
    cat(paste("Using sample:",object@sample,"\n"));
  } else { cat(paste("Using complete measurement\n"));}
  if(length(object@isotopes)>0){
  #Isotopes Vorhanden
    cnt<-nrow(object@isoID)
    cat("Annotated isotopes:",cnt,"\n");
  }
  if(length(object@derivativeIons)>0){
    #Annotations availible
    cnt<-length(unique(object@annoID[,1]));
    cat("Annotated adducts & fragments:",cnt,"\n");
  }
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n");
  if(object@runParallel==1){
    cat("CAMERA runs in parallel mode!\n");
  }
})
###End Constructor###

###xsAnnotate generic Methods###
setGeneric("groupFWHM", function(object, sigma=6, perfwhm=0.6, intval="maxo") standardGeneric("groupFWHM"))
setMethod("groupFWHM","xsAnnotate", function(object, sigma=6, perfwhm=0.6, intval="maxo") {
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

  if(object@groupInfo[1, "rt"] == -1) {
     warning("Warning: no retention times avaiable. Do nothing\n")
  }else{
    if(is.na(sample)) {
      # grouped peaktable within automatic selection
      gvals    <- groupval(object@xcmsSet);
      peakmat  <- object@xcmsSet@peaks;
      groupmat <- groups(object@xcmsSet);

      #calculate highest peaks
      maxo      <- as.numeric(apply(gvals, 1, function(x, peakmat){max(peakmat[x, intval],na.rm=TRUE)}, peakmat));

      #highest peak index 
      int.max   <- as.numeric(apply(gvals, 1, function(x, peakmat){which.max(peakmat[x, intval])}, peakmat));

      peakrange <- matrix(apply(gvals, 1, function(x,peakmat) { peakmat[x[which.max(peakmat[x, intval])], c("rtmin", "rtmax")]}, peakmat), ncol=2, byrow=TRUE); 
      colnames(peakrange) <- c("rtmin", "rtmax")

      while(!all(is.na(maxo) == TRUE)){
          iint   <- which.max(maxo);
          rtmed  <- groupmat[iint, "rtmed"]; #highest peak in whole spectra
          rt.min <- peakrange[iint, "rtmin"];
          rt.max <- peakrange[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt.max-rt.min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
          #all other peaks whose retensiontimes are in the fwhm of the highest peak
          irt    <- which(groupmat[, 'rtmed'] > (rtmed-hwhm) & groupmat[, 'rtmed'] < (rtmed + hwhm) & !is.na(maxo)) 
          if(length(irt) > 0){
              #if peaks are found
              pspectra[[length(pspectra)+1]] <- irt; #create groups
              psSamples[length(pspectra)]  <- int.max[iint] # saves the sample of the peak which is in charge for this pspectrum
              maxo[irt] <- NA; #set itensities of peaks to NA, due to not to be found in the next cycle
          }
      }
    }else{
      #Group with specific sample, using all sample or only a one sample experiment
      peakmat <- getPeaks(object@xcmsSet, index=sample);
      maxo    <- peakmat[, intval]; #max intensities of all peaks

      while(!all(is.na(maxo) == TRUE)){
          iint   <- which.max(maxo);
          rtmed  <- peakmat[iint, "rt"]; #highest peak in whole spectra
          rt.min <- peakmat[iint, "rtmin"];
          rt.max <- peakmat[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt.max - rt.min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
          #all other peaks whose retensiontimes are in the fwhm of the highest peak
          irt    <- which(peakmat[, 'rt'] > (rtmed - hwhm) & peakmat[, 'rt'] < (rtmed + hwhm) & !is.na(maxo)) 
          if(length(irt)>0){
              #if peaks are found
              pspectra[[length(pspectra)+1]] <- irt; #create groups
              maxo[irt] <- NA; #set itensities of peaks to NA, due to not to be found in the next cycle
          }
      }
    psSamples <- rep(sample, length(pspectra))
    }

    object@pspectra  <- pspectra;
    object@psSamples <- psSamples;
    cat("Created", length(object@pspectra), "groups.\n")
  }
  return(invisible(object)); #return object
})

setGeneric("groupCorr",function(object,cor_eic_th=0.75,psg_list=NULL,polarity=NA) standardGeneric("groupCorr"));
setMethod("groupCorr","xsAnnotate", function(object,cor_eic_th=0.75,psg_list=NULL,polarity=NA) {
  if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
  #Is LC data is available?
  if(object@xcmsSet@peaks[1,"rt"] == -1) {
     warning("Warning: no retention times avaiable. Do nothing\n")
  }else {
    npspectra <- length(object@pspectra);
    #Wenn groupFWHM nicht vorher aufgerufen wurde!
    if(npspectra<1){
      npspectra<-1;
      object@pspectra[[1]]<-seq(1:nrow(object@groupInfo));
      if(is.na(object@sample)){
        object@psSamples <- 1;
      }else{
        object@psSamples <- object@sample;
      }
      cat("Calculating peak correlations for all peaks in one big group on sample ",object@psSamples,".\nUse groupFWHM before, to reduce runtime.\n");
    }
    if(is.na(object@sample)){
        index <- rep(0,nrow(object@groupInfo));
        for(i in 1:npspectra){
          index[object@pspectra[[i]]] <- object@psSamples[[i]];
        }
    }else if(object@sample > 0){
        index <- rep(object@sample,nrow(object@groupInfo));
    }else{
      #sample:-1 way, @joe fix pls
      stop("NYI");
    }
#     tmp <- getAllEICs(object@xcmsSet,index=index);
    tmp <- getAllPeakEICs(object@xcmsSet,index=index);
    EIC <- tmp$EIC
    scantimes <- tmp$scantimes
    cnt<-length(object@pspectra);
    rm(tmp);
    if(nrow(object@isoID)>0){
      cat("Isotope annotation found, used as grouping information.\n")
    }
    #new calcCL function with rcorr from Hmisc
    tmp <- calcCL3(object, EIC=t(EIC), scantimes=scantimes, cor_eic_th=cor_eic_th,psg_list=psg_list);
    if(is.null(tmp)){
      #found no subgroups
      cat("No group was seperated.\n")
      return(invisible(object));
    }
    CL<-tmp$CL;
    rm(tmp);
    if(!is.na(polarity)){
        if(polarity %in% c("positive","negative")){
          object@polarity<-polarity
        }else{cat("Unknown polarity parameter.\n"); }
    }
    object <- calc_pc(object, CL, psg_list=psg_list, psSamples=object@psSamples)
    cat("xsAnnotate has now", length(object@pspectra), "groups, instead of", cnt, "\n"); 
  }
return(invisible(object));
})

setGeneric("findIsotopes", function(object, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, intval="into") standardGeneric("findIsotopes"));
setMethod("findIsotopes","xsAnnotate", function(object, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, intval="into"){
  #searches in every pseudospectrum after mass differences, which matches isotope distances

  if (!class(object) == "xsAnnotate"){
    stop ("no xsAnnotate object");
  }

  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!")
  }

  ncl <- sum(sapply(object@pspectra, length));
  npeaks.global <- 0; #Counter for % bar
  npspectra <- length(object@pspectra);

  # calculate Isotope_Matrix
  IM <- calcIsotopeMatrix(maxiso=maxiso, maxcharge=maxcharge);
  # scaling
  devppm <- ppm / 1000000;

  #Check if object have been preprocessed with groupFWHM
  if(npspectra < 1) { 
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    object@psSamples  <- 1;
  }

  # get mz,rt,into from peaktable
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##one sample case
    imz  <- object@groupInfo[, "mz"];
    irt  <- object@groupInfo[, "rt"];
    mint <- object@groupInfo[, intval];
  }else {
    ##multiple sample
    gvals    <- groupval(object@xcmsSet);
    peakmat  <- object@xcmsSet@peaks;
    groupmat <- groups(object@xcmsSet);

    imz <- groupmat[, "mzmed"];
    irt <- groupmat[, "rtmed"];

    #get intensity values, according groupFWHM sample selection
    if(is.na(object@sample)){
      psspec <- 1:npspectra
      mint <- vector(mode="numeric",length=nrow(gvals));
      #get for every psspec its corresponding intensities
      invisible(sapply(psspec, function(x) {
          pi <- object@pspectra[[x]]
          index <- object@psSamples[[x]]
          mint <<- peakmat[gvals[pi,index],intval]
      }));
    }else if(object@sample == -1){
      ##TODO @ Joe: Shot never occur!
    }else{
      mint <- peakmat[gvals[, object@sample], intval]; #errechne höchsten Peaks
    }
  }

  isotope   <- vector("list", length(imz));
  npspectra <- length(object@pspectra);
  isomatrix <- matrix(NA, ncol=5, nrow=1);


  cat("Run isotope peak annotation\n % finished: ");
  #Suche Isotope in jeder Gruppe
  for(i in 1:npspectra){
    #indizes der peaks aus der gruppe in der peaktable
    ipeak <- object@pspectra[[i]];

    #percent output
    npeaks.global <- npeaks.global + length(ipeak);
    perc   <- round((npeaks.global) / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { 
      cat(perc,' '); 
      lp <- perc;
    }
    if (.Platform$OS.type == "windows"){ 
      flush.console();
    }
    #end percent output

    #hat gruppe mehr als einen Peak, sonst mach nichts
    if(length(ipeak) > 1){
      #masse und intensität der Peaks
      mz  <- imz[ipeak];
      int <- mint[ipeak];
      #matrix der peaks mit allen wichtigen Informationen
      spectra <- matrix(c(mz, int, ipeak), ncol=3)
      spectra <- spectra[order(spectra[, 1]), ];
      cnt <- nrow(spectra);
      #für jeden Peak
      for ( j in 1:(length(mz) - 1)){
        #erzeuge Differenzmatrix
        MI <- spectra[j:cnt, 1] - spectra[j, 1];
        #für alle erlaubte Ladungen
        for(charge in maxcharge:1){
          #Suche Übereinstimmungen der Isotopenabständen mit der Differenzmatrix
          m <- fastMatch(MI,IM[,charge],tol= max(2*devppm*mz)+ mzabs)
          #Für jeden Match, teste welches Isotope gefunden wurde
          if(any(!sapply(m,is.null))){
            #für alle erlaubten Isotopenpeaks
            for( iso in 1:maxiso){
              #wurde der iso-Isotopenpeak gefunden?
              pos <- which(sapply(m, function(x){ 
		if(is.null(x)){
		  return(FALSE);
		} else { 
		  x == iso;
		  }
		}));
              if (length(pos) > 0){
                # Isotop Nr. iso scheint zu existieren
                dev <- (devppm * spectra[pos+j-1,1]) + (devppm + spectra[j,1])
                if (isTRUE(all.equal(spectra[pos+j-1,1],spectra[j,1] + IM[iso,charge] ,tolerance=(dev + mzabs),scale=1))){
                  # Isotop Nr. iso existiert
                  int.available <- all(!is.na(c(spectra[pos+j-1,2],spectra[j,2])))
                  if (iso == 1){
                    #wenn der erste Isotopenpeak gefunden wurde
                    if (int.available){
                      ISO_RULE1 <- (spectra[pos+j-1, 2] < spectra[j, 2] ) ## isotopic rule
                      theo.mass <- spectra[j, 1] * charge; #theoretical mass
                      numC      <- round(theo.mass / 12); #max. number of C in molecule
                      inten.max <- spectra[j, 2] * numC * 0.011; #highest possible intensity
                      inten.min <- spectra[j, 2] * 1    * 0.011; #lowest possible intensity
                      ## here C12/C13 rule, isotopic rule now obsolete?
                      ISO_RULE  <- (spectra[pos+j-1, 2] < inten.max && spectra[pos+j-1, 2] > inten.min)
                    } else {
                      ISO_RULE1 <- TRUE
                      ISO_RULE  <- TRUE
                    }
                  }else{
                    # Sind alle anderen isotopen Peaks da?
		    test <- match(apply(isomatrix[, c(1, 3, 4),drop=FALSE], 1, function(x) {paste(x,collapse=" ")}),apply(matrix(cbind(spectra[j,3],1:(iso-1),charge),ncol=3),1, function(x) {paste(x,collapse=" ")}))
                    if(length(naOmit(test))==(iso-1)){
                      ISO_RULE1 <- TRUE
                      if (int.available) ISO_RULE <- (spectra[pos+j-1,2] < spectra[j,2])
                      else ISO_RULE <- TRUE
                    }else{
                      ISO_RULE1 <- FALSE
                    }
                  }
                  if (!ISO_RULE1) { 
                    break;
                  }
                  if (ISO_RULE1 && ISO_RULE){
                    #Neues Isotope gefunden
                    #TODO: Intrinsische Ladungen betrachten
                    if(!length(which(isomatrix[,1]==spectra[j,3] & isomatrix[,2]==spectra[pos+j-1,3]))>0){
                      if(!length(which(isomatrix[,2]==spectra[j,3])>0)){
                        isomatrix<-rbind(isomatrix,c(spectra[j,3],spectra[pos+j-1,3],iso,charge,0))
                      }
                    }
                  }
                } else { 
                  break;
                }
              } else { 
                break;
              }
            }
          }
        }
      }
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
      charges.length <- apply(tmp,1, function(x,isomatrix) { length(which(isomatrix[, 1] == x[1] & isomatrix[,4] == x[2])) },isomatrix);
      idx <- which(charges.length == max(charges.length));
      if(length(idx) == 1){
        #max is unique
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),]
      }else{
        #select this one, which lower charge
        idx <- which.max(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),]
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
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] == charges.list[-idx]),]
      }else{
        #select this one, which lower charge
        idx <- which.min(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] == charges.list[-idx]),]
      }
    }
  }
  object@isoID <- matrix(nrow=0, ncol=4);
  colnames(object@isoID)  <-  c("mpeak","isopeak","iso","charge");
  object@isoID <- rbind(object@isoID, isomatrix[, 1:4]);
  # Zähler für Isotopengruppen
  globalcnt <- 0;
  oldnum    <- 0;
  if(nrow(isomatrix) > 0){
    for( i in 1:nrow(isomatrix)){
      if(!isomatrix[i, 1] == oldnum){
          globalcnt <- globalcnt+1; 
          isotope[[isomatrix[i, 1]]] <- list(y=globalcnt, iso=0, charge=isomatrix[i, 4], val=isomatrix[i, 5]);
          oldnum <- isomatrix[i, 1];
      };
      isotope[[isomatrix[i,2]]]<-list(y=globalcnt,iso=isomatrix[i,3],charge=isomatrix[i,4],val=isomatrix[i,5]);
    }
  }
  cnt<-nrow(object@isoID);
  cat("Found isotopes:",cnt,"\n");
  object@isotopes <- isotope;
  return(object);
})

setGeneric("findAdducts",function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL,rules=NULL,max_peaks=100,psg_list=NULL) standardGeneric("findAdducts"));
setMethod("findAdducts", "xsAnnotate", function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL,rules=NULL,max_peaks=100,psg_list=NULL){
  # norming
  devppm = ppm / 1000000;
  # get mz values from peaklist
 if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##Ein Sample Fall
    imz <- object@xcmsSet@peaks[,"mz"];
  }else {
    ##Mehrsample Fall
    imz <- object@groupInfo[,"mz"];
  }

  # anzahl peaks in gruppen für % Anzeige
  sum_peaks <- sum(sapply(object@pspectra, length));
  # Isotopen
  isotopes  <- object@isotopes;
  # Adduktliste
  derivativeIons <- vector("list", length(imz));
  # Sonstige Variablen
  oidscore <- c();
  index    <- c();
  annoID   <- matrix(ncol=3, nrow=0)
  annoGrp  <- matrix(ncol=4, nrow=0)
  colnames(annoID) <-  c("id", "grp_id", "rule_id");
  colnames(annoGrp)<-  c("id", "mass", "ips", "psgrp");

  if(!(object@polarity == "")){
    cat(paste("Polarity is set in xsAnnotate:", object@polarity, "\n"));
    if(is.null(rules)){
      if(!is.null(object@ruleset)){
        rules <- object@ruleset;
      }else{ 
        cat("Ruleset could not read from object! Recalculate\n");
        rules <- calcRules(maxcharge=3, mol=3, nion=2, nnloss=1, nnadd=1, nh=2, polarity=object@polarity);
        object@ruleset <- rules;
      }
    }else{ cat("Found and use user-defined ruleset!");}
  }else {
    #Erkenne polarität
    if(!is.null(polarity)){
      if(polarity %in% c("positive","negative")){
        if(is.null(rules)){
          rules <- calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=polarity);
        }else{ cat("Found and use user-defined ruleset!");}
          object@polarity <- polarity;
      }else stop("polarity mode unknown, please choose between positive and negative.")
    }else if(length(object@xcmsSet@polarity)>0){
      index <- which(sampclass(object@xcmsSet) == object@category)[1] + object@sample-1;
      if(object@xcmsSet@polarity[index] %in% c("positive","negative")){
        if(is.null(rules)){
          rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=object@xcmsSet@polarity[index]);
        }else{ cat("Found and use user-defined ruleset!");}
        object@polarity <- polarity;
      }else stop("polarity mode in xcmsSet unknown, please define variable polarity.")
  }else stop("polarity mode could not be estimated from the xcmsSet, please define variable polarity!")
    #save ruleset
    object@ruleset<-rules;
  }
  runParallel <- 0;
  if(object@runParallel == 1){
      if(mpi.comm.size() > 0){
        runParallel <- 1;
      }else{
        warning("CAMERA runs in parallel mode, but no slaves are spawned!\nRun in single core mode!\n");
        runParallel <- 0;
        }
    }else{
      runParallel <- 0;
    }

  quasimolion <- which(rules[, "quasi"]== 1)
  #Entferne Isotope aus dem Intensitätsvector, sollen nicht mit annotiert werden
  if(length(isotopes) > 0){
      for(x in 1:length(isotopes)){
          if(!is.null(isotopes[[x]])){
              if(isotopes[[x]]$iso != 0){
                imz[x] <- NA;
              }
          }
      }
  }
  #Anzahl Gruppen
  npspectra <- length(object@pspectra);

  #Wenn vorher nicht gruppiert wurde, alle Peaks in eine Gruppe stecken
  if(npspectra < 1){ 
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo)); 
  }
  
  #zähler für % Anzeige
  npeaks<-0;massgrp<-0;
  #für alle Gruppen

  if (runParallel == 1) { ## ... we use MPI
    if(is.null(psg_list)){
        cat('\nCalculating possible adducts in',npspectra,'Groups... \n % finished: '); lp <- -1;
        pspectra_list<-1:npspectra;
      }else{
        cat('\nCalculating possible adducts in',length(psg_list),'Groups... \n % finished: '); lp <- -1;
        pspectra_list<-psg_list;
      }
        argList <- list();
        cnt_peak<-0;if(is.null(max_peaks)) {max_peaks==100;};
        params <- list();
        for(j in 1:length(pspectra_list)){
          i<-pspectra_list[j];
          params$i[[length(params$i)+1]] <- i;
          cnt_peak<-cnt_peak+length(object@pspectra[[i]]);
          if(cnt_peak>max_peaks || j == length(pspectra_list)){
            params$pspectra <-object@pspectra;
            params$imz <- imz;
            params$rules <- rules;
            params$mzabs <- mzabs;
            params$devppm <- devppm;
            params$isotopes <- isotopes;
            params$quasimolion <- quasimolion;
            argList[[length(argList)+1]]<-params
            cnt_peak<-0;params<-list();
          }
        }
        result <- xcmsPapply(argList, annotateGrpMPI)
    for(ii in 1:length(result)){
        if(length(result[[ii]])==0)next;
        for(iii in 1:length(result[[ii]])){
          hypothese<-result[[ii]][[iii]];
          if(is.null(hypothese)){next;}
          charge=0;old_massgrp=0;
          index <- argList[[ii]]$i[[iii]];
          ipeak <- object@pspectra[[index]];
          for(hyp in 1:nrow(hypothese)){
            peakid<-ipeak[hypothese[hyp,"massID"]];
            if(old_massgrp != hypothese[hyp,"massgrp"]) {
              massgrp<-massgrp+1;old_massgrp<-hypothese[hyp,"massgrp"];
              annoGrp<-rbind(annoGrp,c(massgrp,hypothese[hyp,"mass"],sum(hypothese[ which(hypothese[,"massgrp"]==old_massgrp),"ips"]),i) ) 
            }
            annoID<-rbind(annoID, c(peakid,massgrp,hypothese[hyp,"ruleID"]))
          }
        }
    }
    derivativeIons<-getderivativeIons(annoID,annoGrp,rules,length(imz));
    cat("\n");
    object@derivativeIons <- derivativeIons;
    object@annoID<-annoID;
    object@annoGrp<-annoGrp;
    return(object)
  } else {
    if(is.null(psg_list)){
      cat('\nCalculating possible adducts in',npspectra,'Groups... \n % finished: '); lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating possible adducts in',length(psg_list),'Groups... \n % finished: '); lp <- -1;
      pspectra_list <- psg_list;
      sum_peaks <- sum(sapply(object@pspectra[psg_list],length));
      }

    for(j in 1:length(pspectra_list)){
      i <- pspectra_list[j];
      #Indizes der Peaks in einer Gruppe
      ipeak <- object@pspectra[[i]];
      #Zähler hochzählen und % ausgeben
      npeaks <- npeaks + length(ipeak);
      perc   <- round((npeaks) / sum_peaks * 100)
      if ((perc %% 10 == 0) && (perc != lp)) { 
        cat(perc,' '); 
        lp <- perc; 
      }
      if (.Platform$OS.type == "windows"){
        flush.console();
      }
      #check if the pspec contains more than one peak 
      if(length(ipeak) > 1){
        hypothese <- annotateGrp(ipeak,imz,rules,mzabs,devppm,isotopes,quasimolion);
        #save results
        if(is.null(hypothese)){
          next;
        }
        charge <- 0;old_massgrp <- 0;
        for(hyp in 1:nrow(hypothese)){
            peakid<-ipeak[hypothese[hyp,"massID"]];
            if(old_massgrp != hypothese[hyp,"massgrp"]) {
            massgrp<-massgrp+1;old_massgrp<-hypothese[hyp,"massgrp"];
            annoGrp<-rbind( annoGrp,c(massgrp,hypothese[hyp,"mass"],sum(hypothese[ which(hypothese[,"massgrp"]==old_massgrp),"ips"]),i) ) }
            annoID<-rbind(annoID, c(peakid,massgrp,hypothese[hyp,"ruleID"]))
        }
      }
    }
    derivativeIons<-getderivativeIons(annoID,annoGrp,rules,length(imz));
    cat("\n");
    object@derivativeIons <- derivativeIons;
    object@annoID<-annoID;
    object@annoGrp<-annoGrp;
    return(object)
  }
})

annotateDiffreport <- function(object, sample=NA,sigma=6, perfwhm=0.6, cor_eic_th=0.75, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, multiplier=3, polarity="positive", nSlaves=1, psg_list=NULL, pval_th=NULL, fc_th=NULL, quick=FALSE,rules=NULL, class1 = levels(sampclass(object))[1], class2 = levels(sampclass(object))[2], filebase = character(), eicmax = 0, eicwidth = 200, sortpval = TRUE, classeic = c(class1,class2), value=c("into", "maxo", "intb"), metlin = FALSE, h=480,w=640, ...) {

  if (!class(object)=="xcmsSet") stop ("no xcmsSet object");
  diffrep <- diffreport(object, class1 = class1, class2 = class2, filebase = filebase, eicmax = eicmax, eicwidth = eicwidth, sortpval = FALSE, classeic = classeic, value=value, metlin = metlin, h=h,w=w, ...);
  if(quick){
    #Quick run, no groupCorr and findAdducts
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa,perfwhm=perfwhm,sigma=sigma);
    xa <- findIsotopes(xa,maxcharge=maxcharge,maxiso=maxiso,ppm=ppm,mzabs=mzabs)
    xa.result<-getPeaklist(xa);
  }else{
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa,perfwhm=perfwhm,sigma=sigma);
    xa <- findIsotopes(xa,maxcharge=maxcharge,maxiso=maxiso,ppm=ppm,mzabs=mzabs)
    if(is.null(psg_list) & is.null(pval_th) & is.null(fc_th)){
      #keine Liste vorgegeben!
      #Werde alle gruppen berechent, psg_list=NULL
    }else{
      #Ein Wert wurde vorgegeben
      #Generiere psg_list
      peaklist<-getPeaklist(xa);
      if(is.null(psg_list)){
        psg_list<-1:length(xa@pspectra);
      }
      if(!is.null(pval_th)){
        #Finde alle Gruppen, welche Feature enthalten, mit p-val < pval_th
        index <- which(diffrep[, "pvalue"] < pval_th);
        index.grp <- unique(peaklist[index, "pcgroup"]);
        if(length(index.grp)>0){
          psg_list  <- psg_list[which(psg_list %in% index.grp)]
        }else{
          cat("No groups found, which satisfy your conditions!\n")
          result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
          return(result);
        }
      }
      if(!is.null(fc_th)){
        #Finde alle Gruppen, welche Feature enthalten, mit fc > pval_th
        index <- which(diffrep[, "fold"] > fc_th);
        index.grp <- unique(peaklist[index, "pcgroup"]);
        if(length(index.grp)>0){
          psg_list  <- psg_list[which(psg_list %in% index.grp)]
        }else{
          cat("No groups found, which satisfy your conditions!\n")
          result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
          return(result);
        }
      }
        if(length(psg_list) < 1){
          #Keine Grp mehr uebrig
          cat("No groups found, which satisfy your conditions!\n")
          result <- cbind(diffrep, peaklist[, c("isotopes","adduct","pcgroup")])
          return(result);
        }
    }
    #Add to psg_list all groups, with has been created after groupCorr
    cnt <- length(xa@pspectra);
    xa <- groupCorr(xa,cor_eic_th=cor_eic_th,psg_list=psg_list, polarity=polarity)
    if(!is.null(psg_list)){
      psg_list <- c(psg_list,(cnt+1):length(xa@pspectra));
    }
    xa <- findAdducts(xa,multiplier=multiplier,ppm=ppm,mzabs=mzabs,polarity=polarity,rules=rules,psg_list=psg_list);
    xa.result<-getPeaklist(xa);
  }
  #Kombiniere Resultate

  result <- cbind(diffrep, xa.result[, c("isotopes","adduct","pcgroup")])
  if(sortpval){
    #Sortierung wieder herstellen
    result <- result[order(result[,"pvalue"]),];
  }
  return(result);
} 

###End xsAnnotate generic Methods###

###xsAnnotate exported Methods###

findFragment <- function (object,ppm=20){
  if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
  #number pseudospectra
  devppm = ppm / 1000000;
  npspectra <- length(object@pspectra);
  if (object@polarity != "positive" & object@polarity != "negative") stop ("xsAnnotate object have wrong polarities.\nOnly pos/neg is allowed!")
  neutralloss <- system.file('lists/neutralloss.csv', package = "CAMERA")[1]
  if (!file.exists(neutralloss)) stop('neutralloss.csv not found.')
  neutralloss <- read.table(neutralloss, header=TRUE, dec=".", sep=",",
                            as.is=TRUE, stringsAsFactors = FALSE);
  colnames(neutralloss)<-c("name","massdiff")
  fragment<-rep("",nrow(object@groupInfo))
  for(i in 1:npspectra){
   print (i);
    index<-object@pspectra[[i]];
    peaktable<-object@groupInfo[index,];
    if(!is.matrix(peaktable)) peaktable<-matrix(peaktable,ncol=11);
    mz <- peaktable[,1];

    ML <- massDiffMatrixNL(mz,neutralloss)
    m <- apply(ML,2,function (x) {fastMatch(mz,x,tol = max(2*devppm*mean(mz,na.rm=TRUE)))})
    indi<-which(sapply(m,function (x) {length(unlist(x))})>0)
    if(length(indi)==0) next;
    for( ii in 1:length(indi)){
      for(iii in 1:length(mz)){
          if(!is.null(m[[indi[ii]]][[iii]])){
            if(length(fragment[index[iii]])>1) {
              fragment[index[iii]]<-paste(mz[m[[indi[ii]]][[iii]]],"-",neutralloss[indi[ii],"name"])
            }else{
              fragment[index[iii]]<-paste(fragment[index[iii]], mz[m[[indi[ii]]][[iii]]],"-",neutralloss[indi[ii],"name"])
            }
        }
      }
    }
  }
  invisible(cbind(getPeaklist(object),fragment))
}

getpspectra <- function(object,grp){
peaks<-object@groupInfo
index<-object@pspectra[[grp]];
peaktable<-peaks[index,]
adduct<-vector("character",length(index));
isotopes<-vector("character",length(index));
ions<-object@derivativeIons;
iso<-object@isotopes;
lions <- ions[index];
liso  <- iso[index];

for(i in 1:length(lions)){
    if(!is.null(lions[[i]])){
        if(length(lions[[i]])>1){
            names<-c();
            for(ii in 1:length(lions[[i]])){
                names<-paste(names,lions[[i]][[ii]]$name,lions[[i]][[ii]]$mass);
            }
            adduct[i]<-names;
        }else{
            adduct[i]<-paste(lions[[i]][[1]]$name,lions[[i]][[1]]$mass);
        }
    }
    if(!is.null(liso[[i]])){
	if(liso[[i]]$charge>1){
	  isotopes[i]<-paste("[",liso[[i]]$y,"] ",liso[[i]]$iso," ",liso[[i]]$charge,"+",sep="");
	}else{
	  isotopes[i]<-paste("[",liso[[i]]$y,"] ",liso[[i]]$iso," ","+",sep="");
	}
    }
}
if(is.null(nrow(peaktable))) peaktable = matrix(peaktable,byrow=F,ncol=length(peaktable))
colnames(peaktable)<-colnames(peaks)

return(invisible(data.frame(peaktable,isotopes,adduct,grp,stringsAsFactors=FALSE)));
}


setGeneric("getPeaklist", function(object, intval="into") standardGeneric("getPeaklist"))
setMethod("getPeaklist", "xsAnnotate", function(object, intval="into") {
  
  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!")
  }

  #generate peaktable
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##one sample
    peaktable <- object@groupInfo;
  }else {
    ##multiple sample
    #use groupInfo information and replace intensity values
    peaktable <- object@groupInfo;
    grpval <- groupval(object@xcmsSet, value=intval);
    grpval.ncol <- ncol(grpval)
    start <- ncol(peaktable) - grpval.ncol +1;
    ende  <- start + grpval.ncol - 1; 
    peaktable[,start:ende] <- grpval;
  }

  #allocate variables
  adduct   <- vector("character", nrow(object@groupInfo));
  isotopes <- vector("character", nrow(object@groupInfo));
  pcgroup  <- vector("character", nrow(object@groupInfo));

  for(i in 1:length(isotopes)){
    if(length(object@derivativeIons)>0 && !(is.null(object@derivativeIons[[i]]))){
        if(length(object@derivativeIons[[i]])>1){
            names<-paste(object@derivativeIons[[i]][[1]]$name,signif(object@derivativeIons[[i]][[1]]$mass,6));
            for(ii in 2:length(object@derivativeIons[[i]]))
            {
                    names<-paste(names,object@derivativeIons[[i]][[ii]]$name,signif(object@derivativeIons[[i]][[ii]]$mass,6));
            }
            adduct[i]<-names;
        }else{
            adduct[i]<-paste(object@derivativeIons[[i]][[1]]$name,signif(object@derivativeIons[[i]][[1]]$mass,6));
        }
    }else { adduct[i]<-""; }
    if(length(object@isotopes)>0 && !is.null(object@isotopes[[i]])){
        num_iso<-object@isotopes[[i]]$iso;
        if(num_iso==0){
            str_iso <- "[M]";
        }else { str_iso<-paste("[M+",num_iso,"]",sep="")}
    if(object@isotopes[[i]]$charge>1){
      isotopes[i] <- paste("[",object@isotopes[[i]]$y,"]",str_iso,object@isotopes[[i]]$charge,"+",sep="");
    }else{
      isotopes[i] <- paste("[",object@isotopes[[i]]$y,"]",str_iso,"+",sep="");
    }
      }else { isotopes[i]<-""; }
  }
  if(length(object@pspectra) < 1){
      pcgroup <- 0;
  }else{
    for(i in 1:length(object@pspectra)){
        index<-object@pspectra[[i]];
        pcgroup[index]<-i;
    }
  }
  rownames(peaktable)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
  return(invisible(data.frame(peaktable,isotopes,adduct,pcgroup,stringsAsFactors=FALSE,row.names=NULL)));
})

annotate<-function(object, sigma=6, perfwhm=0.6, cor_eic_th=0.75, maxcharge=3, maxiso=4, ppm=5, mzabs=0.015, multiplier=3, sample=NA, quick=FALSE, psg_list=NULL, polarity="positive", nSlaves=1, max_peaks=100){
  if (!class(object)=="xcmsSet"){
    stop ("Object is not an xcmsSet object")
  }
  if(quick){
    #Quick run, no groupCorr and findAdducts
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa,perfwhm=perfwhm, sigma=sigma);
    xa <- findIsotopes(xa,maxcharge=maxcharge, maxiso=maxiso, ppm=ppm,mzabs=mzabs)
#     xa.result<-getPeaklist(xa);
  }else{
    xa <- xsAnnotate(object, sample=sample, nSlaves=nSlaves);
    xa <- groupFWHM(xa,perfwhm=perfwhm,sigma=sigma);
    xa <- findIsotopes(xa,maxcharge=maxcharge,maxiso=maxiso,ppm=ppm,mzabs=mzabs)
    cnt <- length(xa@pspectra);
    xa <- groupCorr(xa,cor_eic_th=cor_eic_th,psg_list=psg_list, polarity=polarity)
    if(!is.null(psg_list)){
      psg_list <- c(psg_list,(cnt+1):length(xa@pspectra));
    }
    xa <- findAdducts(xa,multiplier=multiplier,ppm=ppm,mzabs=mzabs,polarity=polarity,psg_list=psg_list);
#     xa.result<-getPeaklist(xa);
  }
  #Kombiniere Resultate

  return(xa);
}

findNeutralLossSpecs <- function(object, mzdiff=NULL, mzabs=0, mzppm=10) {

  nlmin <- mzdiff-mzabs
  nlmax <- mzdiff+mzabs

  hits <- sapply(1:length(object@pspectra), function(j) {
    spec  <-  getpspectra(object, j)

    mz  <-  spec[,1] #mz
    nl  <-  as.matrix(dist(mz, method="manhattan"))
      
    length(which(nl > nlmin & nl < nlmax))>0
  })
}


findNeutralLoss <- function(object, mzdiff=NULL, mzabs=0, mzppm=10) {

  nlmin <- mzdiff-mzabs
  nlmax <- mzdiff+mzabs

  xs <- new("xcmsSet");

  peaks <- matrix(ncol=11, nrow=0);
  colnames(peaks) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                       "into", "intb", "maxo", "sn", "sample");
  sampnames(xs) <- c("NeutralLoss")
  sampclass(xs) <- c("NeutralLoss")

  for(j in 1:length(object@pspectra)) {
    
    spec  <-  getpspectra(object, j)
    
    mz  <-  spec[,1] #mz
    nl  <-  as.matrix(dist(mz, method="manhattan"))
      
    hits <- which(nl > nlmin & nl < nlmax,
                  arr.ind = TRUE)

    if(length(hits) > 0) {
      for (f in 1:nrow(hits)) {
        hit <- as.vector( hits[f,])
        
        ## Remove hit from lower diagonal matrix
        if (hit[1] > hit[2]) {
          next;
        }

        lastcol <- ncol(spec);
        if (mz[hit[1]] > mz[hit[2]]) {
          newpeak <- spec[hit[1],c(1:10,lastcol), drop=FALSE]
        } else {
          newpeak <- spec[hit[2],c(1:10,lastcol), drop=FALSE]
        }
        rownames(newpeak) <- as.character(j)
        peaks <- rbind(peaks, newpeak)
      }

    }
  }
  colnames(peaks) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                       "into", "intb", "maxo", "sn", "pseudospectrum");      
  peaks(xs) <- as.matrix(peaks)
  invisible(xs)
}

###End xsAnnotate exported Methods###

###xsAnnotate intern Function###

# resolve_Adducts <- function(peakgrp){
#     neutralloss <- system.file('lists/neutralloss.csv', package = "CAMERA")[1]
#     if (!file.exists(neutralloss)) stop('neutralloss.csv not found.')
#     neutralloss <- read.table(neutralloss, header=TRUE, dec=".",sep=",",
#                               as.is=TRUE, stringsAsFactors = FALSE);
# }



combine_xsanno <- function(xsa.pos, xsa.neg, pos=TRUE, tol=2, ruleset=NULL){
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
    stop ("xsAnnotate object have bad polarities.\nOnly positive or negative are allowed!")
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
    ruleset <- data.frame(name, 1, 1, 1,1,2.014552, 1, 1); 
    colnames(ruleset) <- c("name", "nmol.pos", "nmol.neg","charge.pos","charge.neg", "massdiff", "ID.pos","ID.neg");
  }else{
    #generate rules stated in ruleset
    if(!is.matrix(ruleset) | ncol(ruleset) != 2){
      stop("Ruleset is no matrix or number of cols are unequal two");
    }else{
          name <- c(); 
          nmol.pos <- c(); nmol.neg <- c();
          charge.pos <-c(); charge.neg <- c();
          massdiff<-c();
         apply(ruleset,1, function(x) {
            name <<- c(name,paste(rules.pos[x[1],"name"],rules.neg[x[2],"name"],sep="/"));
            nmol.pos <<- c(nmol.pos,rules.pos[x[1],"nmol"]);
            nmol.neg <<- c(nmol.neg,rules.neg[x[2],"nmol"]);  
            charge.pos <<- c(charge.pos,rules.pos[x[1],"charge"]);
            charge.neg <<- c(charge.neg,rules.neg[x[2],"charge"]);
            massdiff <<- c(massdiff, rules.pos[x[1],"massdiff"] - rules.neg[x[2],"massdiff"]);
            return();
        })
        ruleset <- data.frame(name,nmol.pos,nmol.neg,charge.pos,charge.neg,massdiff,ruleset);
        colnames(ruleset) <- c("name", "nmol.pos", "nmol.neg","charge.pos","charge.neg", "massdiff", "ID.pos","ID.neg");
    }
  }

  #allocate variable for matching results 
  endresult <- matrix(ncol=6, nrow=0);
  colnames(endresult) <- c("grp_pos","peak_pos","grp_neg","peak_neg","mass","rule")
   
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

  #for every pseudospectra
  for(i in 1:length(ps.match)){
    #check for match
#     cat(i,"\n");
    if(is.null(ps.match[[i]])){
      #no corresponding match in neg sample
      next;
    }

    #matrix stores matches
    result.matrix <- matrix(ncol=4, nrow=0);
    colnames(result.matrix) <- c("pos peak ID","neg peak ID","j","Rule");
    #get all m/z values from pos pseudospectrum
    grp.pos <- getpspectra(xsa.pos, i);
    m.pos   <- NA;

    #remove m/z values of annotated isotope peaks
    for(k in 1:nrow(grp.pos)){
      if(!is.null(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]])){
        if(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]]$iso > 0){
          m.pos <- append(m.pos, NA);
          next;
        }
      }
      m.pos <- append(m.pos, grp.pos[k, 1])
    }
    m.pos <- m.pos[-1]; #Remove NA
    na.ini.pos <- which(!is.na(m.pos)); #index of non NA values

    #get all m/z hypotheses (if exists)
    if(length(index <- which(annoGrp.pos[,4] == i)) >0){
       masslist.pos <- annoGrp.pos[index,2];
    }else{
       #no annotation exists
       masslist.pos <- NULL;
    }
    #for every matching neg pseudospectra
    for(j in 1:length(ps.match[[i]])){
      #get all m/z values from neg pseudospectrum
      grp.neg <- getpspectra(xsa.neg, ps.match[[i]][j]);
      m.neg   <- NA;
      #remove masses of annotated isotope peaks
      for(k in 1:nrow(grp.neg)){
        if(!is.null(xsa.neg@isotopes[[xsa.neg@pspectra[[ps.match[[i]][j]]][k]]])){
          if(xsa.neg@isotopes[[xsa.neg@pspectra[[ps.match[[i]][j]]][k]]]$iso > 0){
            m.neg <- append(m.neg, NA);
            next;
          }
        }
        m.neg <- append(m.neg, grp.neg[k, 1]);
      }
      m.neg <- m.neg[-1]; #Remove NA
      na.ini.neg <- which(!is.na(m.neg));#index of non NA values

      #get all m/z hypotheses (if exists)
      if(length(index <- which(annoGrp.pos[,4] == ps.match[[i]][j])) >0){
       masslist.neg <- annoGrp.pos[index,2];
      }else{
       #no annotation exists
       masslist.neg <- NULL;
      }

      #match rules against m/z values
      results <- matrix(NA,nrow=length(m.pos),ncol=nrow(ruleset));
      results[na.ini.pos,] <- combineHypothese(naOmit(m.pos), naOmit(m.neg), ruleset, na.ini.pos, na.ini.neg);
      if(!all(is.na(results))){
#       if(length(results)>0){
        #Found Hits between pos and neg datasets with ruleset
#         apply(results,1, function(x) {
        
        for(l in 1:ncol(results)){
            index <- which(!is.na(results[, l]));
            if(length(index) > 0){
                  tmp <- matrix(c(xsa.pos@pspectra[[i]][index],xsa.neg@pspectra[[ps.match[[i]][j]]][results[index,l]]),ncol=2);
                result.matrix <- rbind(result.matrix,cbind(tmp,ps.match[[i]][j],l));
            }
        }
      }
    }
    
    if(! nrow(result.matrix) > 0){
      next; #Found no hit
    }

    for(ii in 1:nrow(result.matrix)){
#         cat(ii);
        if(length(index <- which(annoID.pos[,"id"] == result.matrix[ii,1] )) > 0){
          #Peak has annotation(s)
          if(all(annoID.pos[index,"rule_id"] == ruleset[result.matrix[ii,4],"ID.pos"])){
            #Peak has only one annotation and could be verified
            mass <- 2;
            endresult <- rbind(endresult,c(i,result.matrix[ii,1],result.matrix[ii,3],result.matrix[ii,2],mass,result.matrix[ii,4]));
            grp2save  <- c(grp2save,annoID.pos[index,"grp_id"]);
          }else{
            #Peak has more than one annotation or verfication goes wrong
            grp.save <- which(annoID.pos[index ,"rule_id"] == ruleset[result.matrix[ii,4],"ID.pos"]);
            if(length(grp.save) > 0 ){
              #Save verified annotation and remove from index
              grp2save  <- c(grp2save,annoID.pos[index[grp.save],"grp_id"]);
              index <- index[-grp2save];
              mass <- 2;
              endresult <- rbind(endresult,c(i,result.matrix[ii,1],result.matrix[ii,3],result.matrix[ii,2],mass,result.matrix[ii,4]));
      
            }else{
              #Found new annotation
              mass <- 1;
              endresult <- rbind(endresult,c(i,result.matrix[ii,1],result.matrix[ii,3],result.matrix[ii,2],mass,result.matrix[ii,4]));
            }  
            #delete all other hypotheses
            grp2del  <- c(grp2del,annoID.pos[index,"grp_id"]);
          }
        }else{
          #Peak has no annotation
          #Add annotation according to ruleset
            mass <- 1;
            endresult <- rbind(endresult,c(i,result.matrix[ii,1],result.matrix[ii,3],result.matrix[ii,2],mass,result.matrix[ii,4]));
        }
    }
  }
  #Remove grp2del groups, if they are in grp2save
  grp2del <- unique(grp2del);
  grp2save <- unique(grp2save);
  if(length(grp2del) > 0 & length(grp2save) > 0){
    index <- which(grp2del %in% grp2save);
    if(length(index) > 0){
      grp2del <- grp2del[-index];
    }
  }

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
          rule.ID <- ruleset[endresult[i,6],"ID.pos"];
          annoID.pos  <- rbind(annoID.pos,  c(endresult[i,"peak_pos"],old.grpid,rule.ID))
          mass <- (xsa.pos@groupInfo[endresult[i,"peak_pos"],"mz"]*rules.pos[rule.ID,"charge"] - rules.pos[rule.ID,"massdiff"]) / rules.pos[rule.ID,"nmol"];
          annoGrp.pos <- rbind(annoGrp.pos, c(old.grpid,mass,2,endresult[i,1]))
        }
        add.adducts[endresult[i,"peak_pos"]]<-paste("Found",ruleset[endresult[i,6],1]);
      }
    }
    xsa.pos@derivativeIons <- getderivativeIons(annoID.pos,annoGrp.pos,rules.pos,length(xsa.pos@isotopes)); #save results
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
          rule.ID <- ruleset[endresult[i,5],"ID.neg"];
          annoID.neg  <- rbind(annoID.neg,  c(endresult[i,"peak_neg"],old.grpid,rule.ID))
          mass <- (xsa.neg@groupInfo[endresult[i,"peak_neg"],"mz"]* abs(rules.neg[rule.ID,"charge"]) - rules.neg[rule.ID,"massdiff"]) / rules.neg[rule.ID,"nmol"];
          annoGrp.neg <- rbind(annoGrp.neg, c(old.grpid,mass,2,endresult[i,1]))
        }
        add.adducts[endresult[i,"peak_neg"]]<-paste("Found",ruleset[endresult[i,5],1]);
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
 
combineHypothese <- function(mass.pos,mass.neg,ruleset,ini1,ini2,tol=0.02){
    ##generiere neue Peaklist
    tmp.ruleset <- ruleset[,c("name","nmol.pos","charge.pos","massdiff")];
    colnames(tmp.ruleset) <- c("name","nmol","charge","massdiff")
    ML <- massDiffMatrix(mass.pos,tmp.ruleset)
    ML.match <- apply(ML,2,function(x) {
      m <- fastMatch(x,mass.neg,tol=tol);
      unlist(lapply(m, function(x) { if(is.null(x)){return(NA)}else{return(ini2[x[1]])} })); ##Select first hit, if more than one??
    });
    return(ML.match);
#     results<-list();
#     if(length(m)==0) {return(results);}
#     for(i in 1:length(m))
#     {
#         if(is.null(m[[i]]))next;
#         results[[length(results)+1]]<-c(ini1[i],ini2[m[[i]]])
#     }
#     return(results);
}

naOmit <- function(x) {
    return (x[!is.na(x)]);
}

##create peak table diff
getPeaks_selection <- function(xs,method="medret",value="into"){
if (nrow(xs@groups) > 0) {
     groupmat <- groups(xs)
     ts <- data.frame(cbind(groupmat,groupval(xs, method=method, value=value)),row.names = NULL)
     cnames <- colnames(ts)
     if (cnames[1] == 'mzmed') cnames[1] <- 'mz' else stop ('Peak information ?!?')
     if (cnames[4] == 'rtmed') cnames[4] <- 'rt' else stop ('Peak information ?!?')
     colnames(ts) <- cnames
} else if (length(sampnames(xs)) == 1)
        ts <- xs@peaks
    else stop ('First argument must be a xcmsSet with group information or contain only one sample.')
return(as.matrix(ts))
}


# create peak table
getPeaks <- function(xs,index=1){
if (nrow(xs@groups) > 0) {
        if(index== -1)
        {
# 		peaki <- getPeaksIdxCol(xs,NULL)[,1]
                ts <- xs@peaks;
        }
        else
        {
                peaki <- getPeaksIdxCol(xs,NULL)[,index]
                ts <- xs@peaks[peaki,]
        }
#      groupmat <- groups(xs)
#      ts <- data.frame(cbind(groupmat,groupval(xs, "medret", "into")),row.names = NULL)
#      index<-which(sampclass(xs)==category)[1]+sample-1;
#      ts<- cbind(groupmat,groupval(xs, "medret", "into"))
#      cnames <- colnames(ts)
#      if (cnames[1] == 'mzmed') cnames[1] <- 'mz' else stop ('Peak information ?!?')
#      if (cnames[4] == 'rtmed') cnames[4] <- 'rt' else stop ('Peak information ?!?')
#      colnames(ts) <- cnames

# 	cnames <- colnames(ts)
# 	ts <- cbind(ts,which(xs@peaks[,"sample"]==index));
# 	cnames<-append(cnames,"PeakID")
# 	colnames(ts)<- cnames;
} else if (length(sampnames(xs)) == 1)
        ts <- xs@peaks
    else stop ('First argument must be a xcmsSet with group information or contain only one sample.')
return(as.matrix(ts))
}


###End xsAnnotate intern Function###
