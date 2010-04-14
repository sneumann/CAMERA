###Constructor###
setClass("xsAnnotate",
    representation(
                    groupInfo = "matrix" ,
                    pspectra = "list",
                    psSamples="numeric",
                    isotopes="list",
                    derivativeIons="list",
                    formula="list",
                    sample="numeric",
                    xcmsSet="xcmsSet",
                    ruleset="data.frame",
                    annoID="matrix",
                    annoGrp="matrix",
                    isoID="matrix",
                    polarity="character",
                    runParallel="numeric"),
    prototype(
                    groupInfo= matrix(ncol=0,nrow=0),
                    pspectra = list(),
                    psSamples=NULL,
                    isotopes=list(),
                    derivativeIons=list(),
                    formula=list(),
                    sample=NULL,
                    xcmsSet=NULL,
                    ruleset=NULL,
                    annoID=matrix(ncol=3,nrow=0),
                    annoGrp=matrix(ncol=3,nrow=0),
                    isoID=matrix(ncol=4,nrow=0),
                    polarity="",
                    runParallel=NULL)
            );

xsAnnotate <- function(xs=NULL,sample=NA,nSlaves=1){

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
    cat("Parallel mode is currently not avaible.\nWill be re-enabled with the next CAMERA version!\n");
    ## If MPI is available ...
#     rmpi = "Rmpi"
#     if (require(rmpi,character.only=TRUE) && !is.null(nSlaves)) {
#       if (is.loaded('mpi_initialize')) {
#         #test if not already slaves are running!
#         if(mpi.comm.size() >0){ 
#           warning("There are already intialized mpi slaves on your machine.\nCamera will try to uses them!\n");
#           runParallel<-1;
#         }else{
#           mpi.spawn.Rslaves(nslaves=nSlaves, needlog=FALSE)
#           if(mpi.comm.size() > 1){
#             #Slaves have successfull spawned
#             runParallel<-1;
#           }else{ warning("Spawning of mpi slaves have failed. CAMERA will run without parallelization.\n");}
#         }
#       }else {
#           #And now??
#       }
#     }
  }

  object@runParallel<-runParallel;

  colnames(object@annoID) <-  c("id","grp_id","rule_id");
  colnames(object@annoGrp)<-  c("id","mass","ips");
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
setGeneric("groupFWHM", function(object,sigma=6,perfwhm=0.6) standardGeneric("groupFWHM"))
setMethod("groupFWHM","xsAnnotate", function(object,sigma=6,perfwhm=0.6) {
  # Gruppierung nach fwhm
  # sigma - number of standard deviation arround the mean (6 = 2 x 3 left and right)
  # perfwhm - 0.3;
  if (!class(object) == "xsAnnotate") stop ("no xsAnnotate object")
  sample <- object@sample;
  pspectra <- list();
  psSamples <- NA;
  if(object@groupInfo[1, "rt"] == -1) {
     warning("Warning: no retention times avaiable. Do nothing\n")
  }else{
    if(is.na(sample)) {
      #Gruppierte Peaktable with automatic selection 
      gvals <- groupval(object@xcmsSet);
      peakmat <- object@xcmsSet@peaks;
      groupmat <- groups(object@xcmsSet);
      #errechne höchsten Peaks
      maxo      <- as.numeric(apply(gvals, 1, function(x, peakmat){max(peakmat[x, "maxo"],na.rm=TRUE)}, peakmat));
      #index des höchsten peaks
      max_int   <- as.numeric(apply(gvals, 1, function(x, peakmat){which.max(peakmat[x, "maxo"])}, peakmat));
      peakrange <- matrix(apply(gvals, 1, function(x,peakmat) { peakmat[x[which.max(peakmat[x, "maxo"])], c("rtmin", "rtmax")]}, peakmat), ncol=2, byrow=TRUE); 
      colnames(peakrange) <- c("rtmin", "rtmax")
      while(!all(is.na(maxo) == TRUE)){
          iint   <- which.max(maxo);
          rtmed  <- groupmat[iint, "rtmed"]; #highest peak in whole spectra
          rt_min <- peakrange[iint, "rtmin"];
          rt_max <- peakrange[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt_max-rt_min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
          #all other peaks whose retensiontimes are in the fwhm of the highest peak
          irt    <- which(groupmat[, 'rtmed'] > (rtmed-hwhm) & groupmat[, 'rtmed'] < (rtmed + hwhm) & !is.na(maxo)) 
          if(length(irt) > 0){
              #if peaks are found
              pspectra[[length(pspectra)+1]] <- irt; #create groups
              psSamples[length(pspectra)]  <- max_int[iint] # saves the sample of the peak which is in charge for this pspectrum
              maxo[irt] <- NA; #set itensities of peaks to NA, due to not to be found in the next cycle
          }
      }
    }else{
      #Group with specific sample, using all sample or only a one sample experiment
      peakmat <- getPeaks(object@xcmsSet, index=sample);
      maxo    <- peakmat[, 'maxo']; #max intensities of all peaks
      while(!all(is.na(maxo) == TRUE)){
          iint   <- which.max(maxo);
          rtmed  <- peakmat[iint, "rt"]; #highest peak in whole spectra
          rt_min <- peakmat[iint, "rtmin"];
          rt_max <- peakmat[iint, "rtmax"]; #begin and end of the highest peak
          hwhm   <- ((rt_max - rt_min) / sigma * 2.35 * perfwhm) / 2; #fwhm of the highest peak
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
  #check if LC data is available
  if(object@xcmsSet@peaks[1,"rt"] == -1) {
     warning("Warning: no retention times avaiable. Do nothing\n")
  }else {
    tmp <- getAllEICs(object@xcmsSet)
    EIC <- tmp$EIC
    scantimes <- tmp$scantimes
    cnt<-length(object@pspectra);
    rm(tmp);
    if(nrow(object@isoID)>0){
      cat("Isotope annotation found, used as grouping information.\n")
    }
    npspectra <- length(object@pspectra);
    #Wenn groupFWHM nicht vorher aufgerufen wurde!
     if(npspectra<1){
      npspectra<-1;
      object@pspectra[[1]]<-seq(1:nrow(object@groupInfo));
      cat('Calculating peak correlations for 1 big group.\nTry groupFWHM before, to reduce runtime. \n');
     }
    if(object@runParallel==1){
        if(mpi.comm.size() >0){
          tmp<- calcCL2(object, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th)
        }else{
          warning("CAMERA runs in parallel mode, but no slaves are spawned!\nRun in single core mode!\n");
          tmp <- calcCL(object, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th);
          }
    }else{
      tmp <- calcCL(object, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th,psg_list=psg_list)
    }
    if(is.null(tmp)){
      #found no subgroups
      cat("No group was seperated.\n")
      return(invisible(object));
    }
    CL<-tmp$CL;
    CI<-as.matrix(tmp$CI)
    psSamples <- tmp$psSamples;
    cor_matrix<-matrix(NA,ncol=nrow(object@groupInfo),nrow=nrow(object@groupInfo))
    for(i in 1:nrow(CI)){
      cor_matrix[CI[i,1],CI[i,2]]<-CI[i,3]
    }
    rm(tmp);
    if(!is.na(polarity)){
        if(polarity %in% c("positive","negative")){
          object@polarity<-polarity
        }else{cat("Unknown polarity parameter.\n"); }
    }
    object <- calc_pc(object,CL,cor_matrix,psg_list=psg_list,psSamples=psSamples)
    cat("xsAnnotate has now",length(object@pspectra),"groups, instead of",cnt,"\n"); 
  }
return(invisible(object));
})

setGeneric("findIsotopes", function(object,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01) standardGeneric("findIsotopes"));
setMethod("findIsotopes","xsAnnotate", function(object,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01){
  if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object");
 # calculate Isotope_Matrix
  IM <- calcIsotopes(maxiso=maxiso,maxcharge=maxcharge);
  # Normierung
  devppm <- ppm / 1000000;
  # get mz,rt,into from peaktable
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##Ein Sample Fall
    imz <- object@xcmsSet@peaks[,"mz"];
    irt <- object@xcmsSet@peaks[,"rt"];
    mint <- object@xcmsSet@peaks[,"into"];
  }else {
    ##Mehrsample Fall
    #Gibt es Unterschiede
    gvals <- groupval(object@xcmsSet);
    peakmat <- object@xcmsSet@peaks;
    groupmat <- groups(object@xcmsSet)
    imz<-groupmat[,"mzmed"];
    irt<-groupmat[,"rtmed"];
    if(is.na(object@sample)){
      mint <- as.numeric(apply(gvals,1,function(x,peakmat) { max(peakmat[x,"into"])},peakmat)); #errechne höchsten Peaks
    }else if(object@sample== -1){
      ##TODO @ Joe: Was machen wir hier?
    }else{
      #Group mit vorgegebenen Sample
      mint <- peakmat[gvals[,object@sample],"into"]; #errechne höchsten Peaks
    }
  }

  isotope <- vector("list",length(imz));
  npspectra <- length(object@pspectra);
  isomatrix<-matrix(NA,ncol=5);
  
  #wenn vorher nicht groupFWHM aufgerufen wurde, gruppiere alle Peaks in eine Gruppe
  if(npspectra < 1) { npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@groupInfo));}
  
  cat("Run isotope peak annotation\n")
  #Suche Isotope in jeder Gruppe
  for(i in 1:npspectra){
    #indizes der peaks aus der gruppe in der peaktable
    ipeak <- object@pspectra[[i]];
    #hat gruppe mehr als einen Peak, sonst mach nichts
    if(length(ipeak)>1){
      #masse und intensität der Peaks
      mz <- imz[ipeak];int <- mint[ipeak];
      #matrix der peaks mit allen wichtigen Informationen
      spectra<-matrix(c(mz,int,ipeak),ncol=3)
      spectra<-spectra[order(spectra[,1]),];
      cnt<-nrow(spectra);
      #für jeden Peak
      for ( j in 1:(length(mz)-1)){
        #erzeuge Differenzmatrix
        MI <- spectra[j:cnt,1] - spectra[j,1];
        #für alle erlaubte Ladungen
        for(charge in maxcharge:1){
          #Suche Übereinstimmungen der Isotopenabständen mit der Differenzmatrix
          m<- fastMatch(MI,IM[,charge],tol= max(2*devppm*mz)+ mzabs)
          #Für jeden Match, teste welches Isotope gefunden wurde
          if(any(sapply(m,is.null))){
            #für alle erlaubten Isotopenpeaks
            for( iso in 1:maxiso){
              #wurde der iso-Isotopenpeak gefunden?
              pos <- which(sapply(m,function(x){ if(is.null(x))return(FALSE)else{x == iso}}))
              if (length(pos) > 0){
                # Isotop Nr. iso scheint zu existieren
                dev <- (devppm * spectra[pos+j-1,1]) + (devppm + spectra[j,1])
                if (isTRUE(all.equal(spectra[pos+j-1,1],spectra[j,1] + IM[iso,charge] ,tolerance=(dev + mzabs),scale=1))){
                  # Isotop Nr. iso existiert
                  int.available <- all(!is.na(c(spectra[pos+j-1,2],spectra[j,2])))
                  if (iso == 1){
                    #wenn der erste Isotopenpeak gefunden wurde
                    if (int.available)
                    ISO_RULE1 <- (spectra[pos+j-1,2] < spectra[j,2] ) ## isotopic rule
                    else ISO_RULE1 <- TRUE
                    ISO_RULE <- TRUE
                  }else{
                    # Sind alle anderen isotopen Peaks da?
                    test<-match(apply(isomatrix[,c(1,3,4)],1,function(x) {paste(x,collapse=" ")}),apply(matrix(cbind(spectra[j,3],1:(iso-1),charge),ncol=3),1, function(x) {paste(x,collapse=" ")}))
                    if(length(naOmit(test))==(iso-1)){
                      ISO_RULE1 <- TRUE
                      if (int.available) ISO_RULE <- (spectra[pos+j-1,2] < spectra[j,2])
                      else ISO_RULE <- TRUE
                    }else{
                      ISO_RULE1 <- FALSE
                    }
                  }
                  if (!ISO_RULE1) { break;}
                  if (ISO_RULE1 && ISO_RULE){
                    #Neues Isotope gefunden
                    #TODO: Intrinsische Ladungen betrachten
                    if(!length(which(isomatrix[,1]==spectra[j,3] & isomatrix[,2]==spectra[pos+j-1,3]))>0){
                      if(!length(which(isomatrix[,2]==spectra[j,3])>0)){
                        isomatrix<-rbind(isomatrix,c(spectra[j,3],spectra[pos+j-1,3],iso,charge,0))
                      }
                    }
                  }
                } else break;
              } else break;
            }
          }
        }
      }
    }
  }
  #clean isotopes
  isomatrix<-isomatrix[-1,];
  if(is.null(nrow(isomatrix))) { isomatrix = matrix(isomatrix,byrow=F,ncol=length(isomatrix)) }
  object@isoID<-rbind(object@isoID,isomatrix[,1:4]);
  # Zähler für Isotopengruppen
  globalcnt <- 0;oldnum<-0;
  if(nrow(isomatrix)>0){
    for( i in 1:nrow(isomatrix)){
      if(!isomatrix[i,1]==oldnum){
          globalcnt<-globalcnt+1; 
          isotope[[isomatrix[i,1]]]<-list(y=globalcnt, iso=0, charge= isomatrix[i,4], val=isomatrix[i,5]);
          oldnum<-isomatrix[i,1];
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
# Normierung
devppm = ppm / 1000000;
# hole die wichtigen Spalten aus der Peaktable
 if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##Ein Sample Fall
    imz <- object@xcmsSet@peaks[,"mz"];
    irt <- object@xcmsSet@peaks[,"rt"];
    mint <- object@xcmsSet@peaks[,"into"];
  }else {
    ##Mehrsample Fall
    #Gibt es Unterschiede
    gvals <- groupval(object@xcmsSet);
    peakmat <- object@xcmsSet@peaks;
    groupmat <- groups(object@xcmsSet)
    imz<-groupmat[,"mzmed"];
    irt<-groupmat[,"rtmed"];
    if(is.na(object@sample)){
      mint <- as.numeric(apply(gvals,1,function(x,peakmat) { max(peakmat[x,"into"])},peakmat)); #errechne höchsten Peaks
    }else if(sample== -1){
      ##TODO @ Joe: Was machen wir hier?
    }else{
      #Group mit vorgegebenen Sample
      mint <- peakmat[gvals[,sample],"into"]; #errechne höchsten Peaks
    }
  }

# anzahl peaks in gruppen für % Anzeige
sum_peaks<-sum(sapply(object@pspectra,length));
# Isotopen
isotopes <- object@isotopes;
# Adduktliste
derivativeIons<-vector("list",length(imz));
# Sonstige Variablen
oidscore<-c();index<-c();
annoID=matrix(ncol=3,nrow=0)
annoGrp=matrix(ncol=3,nrow=0)
colnames(object@annoID) <-  c("id","grp_id","rule_id");
colnames(object@annoGrp)<-  c("id","mass","ips");

if(!(object@polarity=="")){
  cat(paste("Polarity is set in xsAnnotate:",object@polarity,"\n"));
  if(is.null(rules)){
    if(!is.null(object@ruleset)){
      rules<-object@ruleset;
    }else{ cat("Ruleset could not read from object! Recalculate\n");
      rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=object@polarity);
      object@ruleset<-rules;
    }
  }else{ cat("Found and use user-defined ruleset!");}
}else {

  #Erkenne polarität
  if(!is.null(polarity)){
      if(polarity %in% c("positive","negative")){
          if(is.null(rules)){
            rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=polarity);
          }else{ cat("Found and use user-defined ruleset!");}
          object@polarity=polarity;
      }else stop("polarity mode unknown, please choose between positive and negative.")
  }else if(length(object@xcmsSet@polarity)>0){
      index<-which(sampclass(object@xcmsSet)==object@category)[1]+object@sample-1
      if(object@xcmsSet@polarity[index] %in% c("positive","negative")){
        if(is.null(rules)){
          rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=object@xcmsSet@polarity[index]);
        }else{ cat("Found and use user-defined ruleset!");}
        object@polarity=polarity;
      }else stop("polarity mode in xcmsSet unknown, please define variable polarity.")
  }else stop("polarity mode could not be estimated from the xcmsSet, please define variable polarity!")
  #save ruleset
  object@ruleset<-rules;
}
  runParallel <- 0
  if(object@runParallel==1){
      if(mpi.comm.size() >0){
        runParallel <- 1;
      }else{
        warning("CAMERA runs in parallel mode, but no slaves are spawned!\nRun in single core mode!\n");
        runParallel <- 0;
        }
    }else{
      runParallel <- 0;
    }

  quasimolion<-which(rules[,"quasi"]==1)
  #Entferne Isotope aus dem Intensitätsvector, sollen nicht mit annotiert werden
  if(length(isotopes)>0){
      for(x in 1:length(isotopes)){
          if(!is.null(isotopes[[x]])){
              if(isotopes[[x]]$iso!=0)imz[x]=NA;
          }
      }
  }
  #Anzahl Gruppen
  npspectra <- length(object@pspectra);
  #Wenn vorher nicht gruppiert wurde, alle Peaks in eine Gruppe stecken
  if(npspectra < 1){ npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@groupInfo)); }
  
  #zähler für % Anzeige
  npeaks<-0;massgrp<-0;
  #für alle Gruppen
  if (runParallel==1) { ## ... we use MPI
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
#         if(nSlaves>1){
#         mpi.close.Rslaves()
#         }
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
              annoGrp<-rbind(annoGrp,c(massgrp,hypothese[hyp,"mass"],sum(hypothese[ which(hypothese[,"massgrp"]==old_massgrp),"ips"])) ) 
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
        pspectra_list<-1:npspectra;
      }else{
        cat('\nCalculating possible adducts in',length(psg_list),'Groups... \n % finished: '); lp <- -1;
        pspectra_list<-psg_list;
        sum_peaks<-sum(sapply(object@pspectra[psg_list],length));
      }
    for(j in 1:length(pspectra_list)){
      i<-pspectra_list[j];
      #Indizes der Peaks in einer Gruppe
      ipeak <- object@pspectra[[i]];
      #Zähler hochzählen und % ausgeben
      npeaks<-npeaks+length(ipeak);
      perc <- round((npeaks) / sum_peaks * 100)
      if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
      if (.Platform$OS.type == "windows") flush.console()
      #wenn mehr als ein Peaks in einer Gruppe ist
      if(length(ipeak)>1){
        hypothese<-annotateGrp(object@pspectra,i,imz,rules,mzabs,devppm,isotopes,quasimolion)
        #Speichern
        if(is.null(hypothese)){next;}
        charge=0;old_massgrp=0;
        for(hyp in 1:nrow(hypothese)){
            peakid<-ipeak[hypothese[hyp,"massID"]];
            if(old_massgrp != hypothese[hyp,"massgrp"]) {
            massgrp<-massgrp+1;old_massgrp<-hypothese[hyp,"massgrp"];
            annoGrp<-rbind( annoGrp,c(massgrp,hypothese[hyp,"mass"],sum(hypothese[ which(hypothese[,"massgrp"]==old_massgrp),"ips"])) ) }
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

annotateGrpMPI <- function(params){
library(CAMERA);
result<-list();
  for(ii in 1:length(params$i)){
    result[[ii]]<-CAMERA:::annotateGrp(params$pspectra,params$i[[ii]],params$imz,params$rules,params$mzabs,params$devppm,params$isotopes,params$quasimolion);
  }
return(result);
}

annotateGrp <- function(pspectra,i,imz,rules,mzabs,devppm,isotopes,quasimolion) {
  ipeak <- pspectra[[i]];
  mz <- imz[ipeak];
  na_ini<-which(!is.na(mz))
  ML <- massDiffMatrix(mz[na_ini],rules)
  m <- fastMatch(as.vector(ML),as.vector(ML),tol = max(2*devppm*mean(mz,na.rm=TRUE))+ mzabs)
  c<-sapply(m,length)
  index<-which(c>=2)
  if(length(index)==0){ return(NULL);}
  
  #Erstelle Hypothesen
  hypothese<-create_hypothese(m,index,ML,rules,na_ini)
  if(is.null(nrow(hypothese))){return(NULL);}
  
  #Entferne Hypothesen, welche gegen Isotopenladungen verstossen!
  if(length(isotopes)>0){
      hypothese <-check_isotopes(hypothese,isotopes,ipeak)
  }
  if(nrow(hypothese)<2){return(NULL);};
  
  #Test auf Quasi-Molekülionen
  hypothese <-check_quasimolion(hypothese,quasimolion)
  if(nrow(hypothese)<2){return(NULL);};
  
  #Entferne Hypothesen, welche gegen OID-Score&Kausalität verstossen!
  hypothese <- check_oid_causality(hypothese,rules)
  if(nrow(hypothese)<2){return(NULL);};
  
  #Prüfe IPS-Score
  hypothese <- check_ips(hypothese)
  if(nrow(hypothese)<2){return(NULL);};
  
  return(hypothese);
}

annotateDiffreport <- function(object, sample=NA,sigma=6, perfwhm=0.6, cor_eic_th=0.75, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, multiplier=3, polarity="positive", nSlaves=1, psg_list=NULL, pval_th=NULL, fc_th=NULL, quick=FALSE, class1 = levels(sampclass(object))[1], class2 = levels(sampclass(object))[2], filebase = character(), eicmax = 0, eicwidth = 200, sortpval = TRUE, classeic = c(class1,class2), value=c("into", "maxo", "intb"), metlin = FALSE, h=480,w=640, ...) {

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
    xa <- findAdducts(xa,multiplier=multiplier,ppm=ppm,mzabs=mzabs,polarity=polarity,psg_list=psg_list);
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

getPeaklist<-function(object){
    xs<-object@xcmsSet
    peaklist<-object@groupInfo;
    adduct<-vector("character",nrow(object@groupInfo));
    isotopes<-vector("character",nrow(object@groupInfo));
    pcgroup<-vector("character",nrow(object@groupInfo));
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
    for(i in 1:length(object@pspectra)){
        index<-object@pspectra[[i]];
        pcgroup[index]<-i;
    }
#     if(!is.null(object@grp_info)){
#         isotopes_tmp<-vector("character",nrow(peaklist));adduct_tmp<-vector("character",nrow(peaklist));pcgroup_tmp<-vector("character",nrow(peaklist));
#         for(i in 1:length(isotopes)){
#             isotopes_tmp[object@grp_info[i]]<-isotopes[i];
#             adduct_tmp[object@grp_info[i]]<-adduct[i];
#             pcgroup_tmp[object@grp_info[i]]<-pcgroup[i];
#         }
#         isotopes<-isotopes_tmp;
#         adduct<-adduct_tmp;
#         pcgroup<-pcgroup_tmp;
#     }
    rownames(peaklist)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
    return(invisible(data.frame(peaklist,isotopes,adduct,pcgroup,stringsAsFactors=FALSE,row.names=NULL)));
}

annotate<-function(xs, sigma=6, perfwhm=0.6, cor_eic_th=0.75, maxcharge=3, maxiso=4, ppm=5, mzabs=0.01, multiplier=3, sample=NA, polarity="positive", nSlaves=1, max_peaks=100){
  if (!class(xs)=="xcmsSet"){
    stop ("xs is not an xcmsSet object")
  }
  xs_anno  <- xsAnnotate(xs, sample=sample, nSlaves=nSlaves);
  xs_anno2 <- groupFWHM(xs_anno, sigma=sigma, perfwhm=perfwhm);
  xs_anno3 <- findIsotopes(xs_anno2, maxcharge=maxcharge ,maxiso=maxiso, ppm=ppm, mzabs=mzabs);
  xs_anno4 <- groupCorr(xs_anno3, cor_eic_th=cor_eic_th, polarity=polarity);
  xs_anno5 <- findAdducts(xs_anno4, multiplier=multiplier, ppm=ppm, mzabs=mzabs, polarity=polarity, max_peaks=max_peaks);
  return(xs_anno5);
}

###End xsAnnotate exported Methods###

###xsAnnotate intern Function###

resolve_Adducts <- function(peakgrp){
    neutralloss <- system.file('lists/neutralloss.csv', package = "CAMERA")[1]
    if (!file.exists(neutralloss)) stop('neutralloss.csv not found.')
    neutralloss <- read.table(neutralloss, header=TRUE, dec=".",sep=",",
                              as.is=TRUE, stringsAsFactors = FALSE);
}

getderivativeIons <- function(annoID,annoGrp,rules,npeaks){
    derivativeIons<-vector("list",npeaks);
    charge=0;
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

create_hypothese<-function(m,index,ML,rules,na_ini){
a<-m[index];
b<-a[which(duplicated(a)==FALSE)];
b_length <- sapply(b,length);b_ini <- 1:length(b);
b2 <- b[order(b_length)];b_ini<-b_ini[order(b_length)];
add<-1;
ini_new<-c();
while(length(b2)>0)
{
    ini<-which(sapply(b2,function(x) {all((b2[[1]]) %in% x)})==TRUE);
    if(length(ini)>1){
        ini_new<-append(ini_new,add);
        b2<-b2[-1];
        add<-add+1;
    }else{
    b2<-b2[-1];add<-add+1;}
}
if(length(ini_new>0)){
ini_new<-b_ini[ini_new];
b<-b[-ini_new];}

nrow_b<-length(b);
ncol_b<-sapply(b,length)
nrow_ML<-nrow(ML);
ncol_ML<-ncol(ML);
c<-as.vector(ML);
hypomass<-sapply(b,function(x) {mean(c[x])})
hypo<-matrix(NA,ncol=8)
colnames(hypo)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp")
for(row in 1:nrow_b)
{
    for(col in 1:ncol_b[row])
    {
        adduct<-b[[row]][col]%/%nrow_ML+1;
        mass <- b[[row]][col]%%nrow_ML;if(mass==0){mass<-nrow_ML;adduct<-adduct-1;}
        hypo <- rbind(hypo,c(na_ini[mass],adduct,rules[adduct,"nmol"],rules[adduct,"charge"],hypomass[row],rules[adduct,"oidscore"],rules[adduct,"ips"],row));
    }
}
hypo<-hypo[-1,]
check<-1
hypothese<-matrix(NA,ncol=9)
colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp","check")
hypothese<-cbind(hypo,check);
return(hypothese)
}

check_ips<- function(hypothese){
for(hyp in 1:nrow(hypothese))
{
    if(length( which( hypothese[,"massgrp"] == hypothese[hyp,"massgrp"])) < 2){hypothese[hyp,"check"]=0;}
    if(length(id<-which(hypothese[,"massID"]==hypothese[hyp,"massID"] & hypothese[,"check"]!=0))>1)
    {
        masses<-hypothese[id,"mass"]
        nmasses<-sapply(masses,function(x) { sum(hypothese[which(hypothese[,"mass"]==x),"ips"]) })
        masses<-masses[-which(nmasses==max(nmasses))];
        if(length(masses)>0)
        {
            hypothese[unlist(sapply(masses, function(x) {which(hypothese[,"mass"]==x)})),"check"]=0;
        }
    }
}
hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=9)
colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp","check")
return(hypothese)
}

check_oid_causality <- function(hypothese,rules){
for(hyp in 1:nrow(hypothese)){
    #check nmol
    if(hypothese[hyp,"nmol"]>1){
        check_sure<-TRUE;
        if(hypothese[hyp,"charge"]==1){
            if(length(which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"]))>1){
                for(prof in (hypothese[hyp,"nmol"]-1):1){
                    indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
                    if(length(indi)==0){
                        check_sure<-FALSE;
                        hypothese[hyp,"check"]<-0;next;
                    }
                }
            }
        }else if(abs(hypothese[hyp,"charge"])==hypothese[hyp,"nmol"]){
            if(length(which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"]))>1){
                for(prof in (hypothese[hyp,"nmol"]-1):1){
                    indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
                    if(length(indi)==0){
                        check_sure<-FALSE;
                        hypothese[hyp,"check"]<-0;#next;
                    }
                }
            }
            if(length(indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & abs(hypothese[,"charge"])==hypothese[,"nmol"]))>1){
                massdiff<-rules[hypothese[indi,"ruleID"],"massdiff"]/rules[hypothese[indi,"ruleID"],"charge"]
                if(length(indi_new<-which(duplicated(massdiff)))>0){
                    check_sure<-FALSE;
                    hypothese[hyp,"check"]<-0;
                }
            }
        }
        if(check_sure){hypothese[hyp,"check"]<-1;}
    }else{
        if(hypothese[hyp,"charge"]>1){
            ##todo
        }else{
            #nothing to say
        }
    }
}
hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=9)
colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp","check")
return(hypothese)
}

check_quasimolion <- function(hypothese,quasimolion){
hypomass<-unique(hypothese[,"mass"])
for(mass in 1:length(hypomass))
{
    if(!any(quasimolion %in% hypothese[which(hypothese[,"mass"]==hypomass[mass]),"ruleID"])){
        hypothese[which(hypothese[,"mass"]==hypomass[mass]),"check"]=0;
    }
    else if(is.null(nrow(hypothese[which(hypothese[,"mass"]==hypomass[mass]),]))){
        hypothese[which(hypothese[,"mass"]==hypomass[mass]),"check"]=0;
    }
}
hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=9)
colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp","check")
return(hypothese)
}

check_isotopes <- function(hypothese,isotopes,ipeak){
for(hyp in 1:nrow(hypothese)){
    peakid<-ipeak[hypothese[hyp,1]];
    if(!is.null(isotopes[[peakid]])){
        #Isotope da
        explainable<-FALSE;
        if(isotopes[[peakid]]$charge==abs(hypothese[hyp,"charge"])){
            explainable<-TRUE;
        }
        if(!explainable){
            #delete Rule
            hypothese[hyp,"check"]=0;
        }
    }
}
hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=9)
colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp","check")
return(hypothese)
}

calcRules <- function (maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=NULL){
    name<-c();nmol<-c();charge<-c();massdiff<-c();oidscore<-c();quasi<-c();ips<-c();
    ##Read Tabellen
    ionlist <- system.file('lists/ions.csv', package = "CAMERA")[1]
    if (!file.exists(ionlist)) stop('ionlist.csv not found.')
    ionlist<-read.table(ionlist, header=TRUE, dec=".", sep=",",
                        as.is=TRUE, stringsAsFactors = FALSE);

    neutralloss <- system.file('lists/neutralloss.csv', package = "CAMERA")[1]
    if (!file.exists(neutralloss)) stop('neutralloss.csv not found.')
    neutralloss <- read.table(neutralloss, header=TRUE, dec=".", sep=",",
                              as.is=TRUE, stringsAsFactors = FALSE);

    neutraladdition <- system.file('lists/neutraladdition.csv', package = "CAMERA")[1]
    if (!file.exists(neutraladdition)) stop('neutraladdition.csv not found.')
    neutraladdition <- read.table(neutraladdition,
                                  header=TRUE, dec=".", sep=",",
                                  as.is=TRUE, stringsAsFactors = FALSE);
    ##End Read Tabellen

    ##Erzeuge Regeln
    tmpname<-c();tmpnmol<-c();tmpcharge<-0;tmpmass<-0;tmpips<-0;
    #Molekülionen
    if(polarity=="positive"){
        #Wasserstoff, hard codiert
        for(k in 1:mol){
            if(k==1){str<-"";tmpips<-1;}else{str<-k;tmpips<-0.5};
            name<-append(name,paste("[",str,"M+H]+",sep=""));
            charge<-append(charge,1);
            massdiff<-append(massdiff,1.0076);
            nmol<-append(nmol,k);
            if(k==1) {quasi<-append(quasi,1);} else { quasi<-append(quasi,0); };
            oidscore<-append(oidscore,1);
            ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M+2H]2+",sep=""));
            charge<-append(charge,2);
            massdiff<-append(massdiff,2.0152);
            nmol<-append(nmol,k);quasi<-append(quasi,0);
            oidscore<-append(oidscore,2);
            ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M+3H]3+",sep=""));
            charge<-append(charge,3);
            massdiff<-append(massdiff,3.0228);
            nmol<-append(nmol,k);
            quasi<-append(quasi,0);
            oidscore<-append(oidscore,3);
            ips<-append(ips,tmpips)
            oid<-3;
            for(i in 1:nrow(ionlist)){
                if(ionlist[i,2]<=0) {next;}
                if(ionlist[i,2]==1){
                    name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]2+",sep=""));
                }else{
                    name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]",ionlist[i,2]+1,"+",sep=""));
                }
                charge <- append(charge,ionlist[i,2]+1);
                massdiff <- append(massdiff,ionlist[i,3]+1.0076);
                nmol <- append(nmol,k);
                quasi <- append(quasi,0);
                oidscore <- append(oidscore,oid+i);
                if(tmpips>0.75){
                    ips<-append(ips,0.5)
                }else{
                    ips<-append(ips,tmpips);
                }#Austausch
            }
            oid<-oid+nrow(ionlist);
            coeff<-expand.grid(rep(list(0:nion),nrow(ionlist)))
            if(length(list<-which(ionlist[,2]<=0))>0){
                coeff[,list]<-0;
            }
            coeff<-unique(coeff);
            coeff<-cbind(coeff,rep(0,nrow(coeff)));
            coeff<-coeff[-1,]
            tmp<-NULL;
            for(i in 1:nrow(ionlist)){
                if(ionlist[i,2]<=0)next;
                #Austausch erstmal nur einen pro Ion
                tmp<-rbind(tmp,t(apply(coeff,1,function(x) {x[i]<-x[i]+1;x[nrow(ionlist)+1]<-1;x})));
            }
            coeff<-unique(rbind(coeff,tmp));
            for(i in 1:nrow(coeff)){
                tmpname<-paste("[",str,"M",sep="");
                tmpcharge<-0;tmpmass<-0;
                for(ii in 1:(ncol(coeff)-1)){
                    if(coeff[i,ii]>0){
                        if(coeff[i,ii]>1){
                            tmpname<-paste(tmpname,"+",coeff[i,ii],ionlist[ii,1],sep="");
                        }else{
                            tmpname<-paste(tmpname,"+",ionlist[ii,1],sep="");
                        }
                        tmpcharge<-tmpcharge+coeff[i,ii]*ionlist[ii,2];
                        tmpmass<-tmpmass+coeff[i,ii]*ionlist[ii,3]
                    }
                }
                if(coeff[i,ncol(coeff)]>0){
                    #Austausch hat stattgefunden, einfach bsp 1
                    tmpname<-paste(tmpname,"-H",sep="");
                    tmpcharge<-tmpcharge-1;
                    tmpmass<-tmpmass-1.0076;
                    tmpips<-0.25;
                }
                if(tmpcharge>1){
                    tmpname<-paste(tmpname,"]",tmpcharge,"+",sep="")
                }else{
                    tmpname<-paste(tmpname,"]+",sep="")
                }
                name<-append(name,tmpname)
                charge<-append(charge,tmpcharge)
                massdiff<-append(massdiff,tmpmass)
                nmol <- append(nmol,k);
                oidscore<-append(oidscore,oid+i)
                if(sum(coeff[i,])==1&& k==1){
                    quasi <- append(quasi,1);
                    ips <-append(ips,tmpips);
                }else{
                    quasi <- append(quasi,0);
                    if(tmpips>0.75){
                        ips <- append(ips,0.75);
                    }else{
                        ips <- append(ips,tmpips);
                    }
                }
            }
        }
        oid<-max(oidscore);
        ##Erzeuge Neutral Addition
        index<-which(quasi==1)
        for(i in 1:nrow(neutraladdition)){
            if(length(index2<-which(ionlist[,2]>0))>0){
                for(ii in 1:length(index2)){
                    if(ionlist[index2[ii],2] > 1){
                        name    <-  append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]",abs(ionlist[index2[ii],2]),"+",sep=""));
                    }else{
                        name    <-  append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]+",sep=""));
                    }
                    charge <- append(charge,ionlist[index2[ii],2]);
                    massdiff <- append(massdiff,neutraladdition[i,2]+ionlist[index2[ii],3]);
                    nmol <- append(nmol,1);
                    quasi <- append(quasi,0);
                    oidscore    <-  append(oidscore,oid+1);oid<-oid+1;
                    ips<-append(ips,0.5);
                }
            }
            if(neutraladdition[i,1]=="NaCOOH"){next;}
            name<-append(name,paste("[M+H+",neutraladdition[i,1],"]+",sep=""));
            charge<-append(charge,+1);
            massdiff<- append(massdiff,neutraladdition[i,2]+1.0076);
            nmol<-append(nmol,1);
            quasi<-append(quasi,0)
            oidscore<-append(oidscore,oid+1);oid<-oid+1;
            ips<-append(ips,0.5);
        }
        oid<-max(oidscore);
        ##Erzeuge Neutral loss
        index<-which(quasi==1)
        for(i in 1:nrow(neutralloss)){
            for(ii in 1:maxcharge){
              if(ii > 1){
                name<-append(name,paste("[M+",ii,"H-",neutralloss[i,1],"]",ii,"+",sep=""));
              }else {name<-append(name,paste("[M+H-",neutralloss[i,1],"]+",sep=""));}
              charge<-append(charge,ii);
              massdiff<-  append(massdiff,-neutralloss[i,2]+1.0076*ii);
              nmol<-append(nmol,1);
              quasi<-append(quasi,0)
              oidscore<-append(oidscore,oid+1);oid<-oid+1;
              ips<-append(ips,0.5);
            }
        }
        ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
        if(length(index<-which(ruleset[,"charge"]>maxcharge))>0){
            ruleset<- ruleset[-index,];
        }
    }else if(polarity=="negative"){
        #Wasserstoff, hard codiert
        for(k in 1:mol){
            if(k==1){str<-"";tmpips<-1;}else{str<-k;tmpips<-0.5};
            name<-append(name,paste("[",str,"M-H]-",sep=""));
            charge<-append(charge,-1);massdiff<-append(massdiff,-1.0076);nmol<-append(nmol,k);if(k==1){quasi<-append(quasi,1);}else{quasi<-append(quasi,0);};oidscore<-append(oidscore,1);ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M-2H]2-",sep=""));charge<-append(charge,-2);massdiff<-append(massdiff,-2.0152);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,2);ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M-3H]3-",sep=""));charge<-append(charge,-3);massdiff<-append(massdiff,-3.0228);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,3);ips<-append(ips,tmpips)
            oid<-3;
            for(i in 1:nrow(ionlist)){
                if(ionlist[i,2]>=0){
                    if(ionlist[i,2]>1) {next;}
                    name<-append(name,paste("[",str,"M-2H+",ionlist[i,1],"]-",sep=""));
                    charge <- append(charge,ionlist[i,2]-2);
                    massdiff<- append(massdiff,ionlist[i,3]-(2*1.0076));
                    nmol <- append(nmol,k);
                    quasi <- append(quasi,0);
                    oidscore<-append(oidscore,oid+i);
                    ips<-append(ips,0.25);
                    next;
                }
                if(ionlist[i,2]== -1){
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]2-",sep=""));
                }else{
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]",ionlist[i,2]+1,"-",sep=""));
                }
                charge <- append(charge,ionlist[i,2]-1);
                massdiff<- append(massdiff,ionlist[i,3]-1.0076);
                nmol <- append(nmol,k);
                quasi <- append(quasi,0);
                oidscore<-append(oidscore,oid+i);
                ips<-append(ips,tmpips);
                #Austausch
            }
            oid<-oid+nrow(ionlist);
            coeff<-expand.grid(rep(list(0:nion),nrow(ionlist)))
            if(length(list<-which(ionlist[,2]>=0))>0){
                coeff[,list]<-0;
            }
            coeff<-unique(coeff);
            coeff<-cbind(coeff,rep(0,nrow(coeff)));
            coeff<-coeff[-1,]
            for(i in 1:nrow(coeff)){
                tmpname<-paste("[",str,"M",sep="");
                tmpcharge<-0;tmpmass<-0;
                for(ii in 1:(ncol(coeff)-1)){
                    if(coeff[i,ii]>0){
                        if(coeff[i,ii]>1){
                            tmpname<-paste(tmpname,"+",coeff[i,ii],ionlist[ii,1],sep="");
                        }else{
                            tmpname<-paste(tmpname,"+",ionlist[ii,1],sep="");
                        }
                        tmpcharge<-tmpcharge+coeff[i,ii]*ionlist[ii,2];
                        tmpmass<-tmpmass+coeff[i,ii]*ionlist[ii,3]
                    }
                }
                if(coeff[i,ncol(coeff)]>0){
                    #Austausch hat stattgefunden, einfach bsp 1
                    tmpname<-paste(tmpname,"-H",sep="");
                    tmpcharge<-tmpcharge-1;
                    tmpmass<-tmpmass-1.0076;
                    tmpips<-0.5;
                }
                if(tmpcharge< -1){
                    tmpname<-paste(tmpname,"]",abs(tmpcharge),"-",sep="")
                }else{
                    tmpname<-paste(tmpname,"]-",sep="")
                }
                name<-append(name,tmpname)
                charge<-append(charge,tmpcharge)
                massdiff<-append(massdiff,tmpmass)
                nmol <-append(nmol,k);
                oidscore<-append(oidscore,oid+i)
                if(sum(coeff[i,])==1&& k==1){
                    quasi   <-append(quasi,1);
                }else{
                    quasi <-append(quasi,0);
                }
                ips <-append(ips,tmpips);
            }
        }
        oid<-max(oidscore);
        ##Erzeuge Neutral Addition
        index<-which(quasi==1)
        for(i in 1:nrow(neutraladdition)){
            if(length(index2<-which(ionlist[,2]<0))>0){
                for(ii in 1:length(index2)){
                    if(ionlist[index2[ii],2]< -1){
                        name <- append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]",abs(ionlist[index2[ii],2]),"-",sep=""));
                    }else{
                        name <- append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]-",sep=""));
                    }
                    charge <- append(charge,ionlist[index2[ii],2]);
                    massdiff<- append(massdiff,neutraladdition[i,2]+ionlist[index2[ii],3]);
                    nmol <- append(nmol,1);
                    quasi <- append(quasi,0);
                    oidscore<-append(oidscore,oid+1);oid<-oid+1;
                    ips<-append(ips,0.5);
                }
            }
            name<-append(name,paste("[M-H+",neutraladdition[i,1],"]-",sep=""));
            charge<-append(charge,-1);
            massdiff<- append(massdiff,neutraladdition[i,2]-1.0076);
            nmol<-append(nmol,1);
            quasi<-append(quasi,0)
            oidscore<-append(oidscore,oid+1);oid<-oid+1;
            ips<-append(ips,0.5);
        }
         oid<-max(oidscore);
        ##Erzeuge Neutral loss
        index<-which(quasi==1)
        for(i in 1:nrow(neutralloss)){
            name<-append(name,paste("[M-H-",neutralloss[i,1],"]+",sep=""));
            charge<-append(charge,+1);
            massdiff<-  append(massdiff,-neutralloss[i,2]-1.0076);
            nmol<-append(nmol,1);
            quasi<-append(quasi,0)
            oidscore<-append(oidscore,oid+1);oid<-oid+1;
            ips<-append(ips,0.5);
        }
        ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
        if(length(index<-which(ruleset[,"charge"]< -maxcharge))>0){
            ruleset<- ruleset[-index,];
        }
    }else stop("Unknown error")
return(ruleset);
}

add_same_oidscore <-function(hypo,adducts,adducts_no_oid){
        hypo_new<-matrix(NA,ncol=6)
        colnames(hypo_new)<-c("massID","ruleID","nmol","charge","mass","oidscore")
        ids<-hypo[,"ruleID"];
        for(i in 1:nrow(hypo))
        {
                index<-which(adducts[,"oidscore"]==adducts_no_oid[ids[i],"oidscore"])
                hypo_new<-rbind(hypo_new,matrix(cbind(hypo[i,"massID"],index,adducts[index,"nmol"],adducts[index,"charge"],hypo[i,"mass"],adducts[index,"oidscore"]),ncol=6))
        }
        hypo_new<-hypo_new[-1,];
        return(hypo_new);
}

combine_xsanno <- function(xsa.pos, xsa.neg, pos=TRUE, tol=2){
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
    stop ("xsAnnotate object have bad polarities.\nOnly pos/neg are allowed!")
  }
  
  ##1. Step
  #get all rts for every pseudospectra
  rt1 <- sapply(xsa.pos@pspectra, function(x) { mean(xsa.pos@groupInfo[x, "rt"]) })
  rt2 <- sapply(xsa.neg@pspectra, function(x) { mean(xsa.neg@groupInfo[x, "rt"]) })
  
  #find matching pseudospectra
  m <- CAMERA:::fastMatch(rt1, rt2, tol=tol)
  #ruleset
  #atm only one rule
  ##TODO: Add more
  rule.hh <- data.frame("[M-H+H]", 1, 1, 2.0152, 1, 1, 1)
  colnames(rule.hh) <- c("name", "nmol", "charge", "massdiff", "oidscore", "quasi", "ips")

  #allocate variable for matching results 
  endresult <- matrix(ncol=6, nrow=0);
  colnames(endresult) <- c("grp_pos","peak_pos","grp_neg","peak_neg","mass","check")
  
  #find indicies of M+H and M-H rule
  rules.pos <- xsa.pos@ruleset;
  rules.neg <- xsa.neg@ruleset;
  id.h.pos  <- which(rules.pos[,"oidscore"] == 1 & rules.pos[,"quasi"] == 1);
  id.h.neg  <- which(rules.neg[,"oidscore"] == 1 & rules.neg[,"quasi"] == 1);

  annoID.pos  <- xsa.pos@annoID;
  annoID.neg  <- xsa.neg@annoID;
  annoGrp.pos <- xsa.pos@annoGrp;
  annoGrp.neg <- xsa.neg@annoGrp;

  grp2del<-c();

  #for every pos pseudospectra
  for(i in 1:length(m)){
    #check for match
    if(is.null(m[[i]])){
      #no corresponding match in neg sample
      next;
    }

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
    hypothese.pos <- grp.pos[, "adduct"];
    if(length(index <- which(hypothese.pos != "")) > 0){
      liste <- unlist(strsplit(hypothese.pos[index], " "));
      masslist.pos <- as.numeric(unique(liste[-grep('\\+', liste)]));
      masslist.pos <- naOmit(masslist.pos);
    }else{
      #no annotation exists
      masslist.pos <- NULL;
    }
    cat("i: ",i,"\nj: ");
    #for every matching neg pseudospectra
    for(j in 1:length(m[[i]])){
      #get all m/z values from neg pseudospectrum
      grp.neg <- getpspectra(xsa.neg, m[[i]][j]);
      m.neg<-NA;
      cat(" ",j);
      #remove masses of annotated isotope peaks
      for(k in 1:nrow(grp.neg)){
        if(!is.null(xsa.neg@isotopes[[xsa.neg@pspectra[[m[[i]][j]]][k]]])){
          if(xsa.neg@isotopes[[xsa.neg@pspectra[[m[[i]][j]]][k]]]$iso > 0){
            m.neg <- append(m.neg, NA);
            next;
          }
        }
        m.neg <- append(m.neg, grp.neg[k, 1]);
      }
      m.neg <- m.neg[-1]; #Remove NA
      na.ini.neg <- which(!is.na(m.neg));#index of non NA values

    #get all m/z hypotheses (if exists)
      hypothese.neg <- grp.neg[, "adduct"];
      if(length(index <- which(hypothese.neg != "")) > 0){
        liste <- unlist(strsplit(hypothese.neg[index], " "))
        masslist.neg <- as.numeric(unique(liste[-grep('\\-', liste)]));
        masslist.neg <- naOmit(masslist.neg)
      }else{
        #no annotation exists
        masslist.neg <- NULL;
      }
      
      #match rules against m/z values
      results <- combine_hypothese(naOmit(m.pos), naOmit(m.neg), rule.hh, na.ini.pos, na.ini.neg);

      if(length(results)>0){
        #M+H,M-H gefunden, Pfad: 1
        anno.ids.pos <- xsa.pos@pspectra[[i]][sapply(results, function(x) {x[1]})] #peak index of pos match
        anno.ids.neg <- xsa.neg@pspectra[[m[[i]][j]]][sapply(results, function(x) {x[2]})] #peak index of neg match
        
        for(ii in 1:length(results)){
          mass <- mean(c(grp.pos[results[[ii]][1], 1], grp.neg[results[[ii]][2], 1])); #putative mass for hit
          if(pos){
            #returns pos peaklist
            id.pos <- xsa.pos@pspectra[[i]][results[[ii]][1]]; #peak index of pos match
            if(length(xsa.pos@derivativeIons[[id.pos]]) > 0){ #Pfad: 2
              #pos peak match has adduct annotation
              if(length(index.pos <- which(annoID.pos[, 1] == id.pos & annoID.pos[, 3] == id.h.pos)) > 0){ #Pfad: 3
                #3.1
                #adduct annotation is M+H
                if(length(index.pos) > 1){
                  #peak has more than one M+H annotation ???
                  warning(paste("index.pos are greater than allowed. Debug! i,j:", i, j));
                }
                if(!length(del.hypo <- which(annoID.pos[, 1] == id.pos & annoID.pos[, 3] != id.h.pos))>0) { #Pfad 3.1
                  #Pfad: 3.1.1, no other hypotheses found
                  #confirm M+H hypothese
                  #add 2 to endresult, means "confirm M+H annotation"
                  endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 2));
                }else{
                  #Pfad: 3.1.2, other hypotheses found
                  #confirm M+H hypothese
                  #delete all other hypotheses
                  ##TODO: DEL
                  endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 2));
                  for(iii in 1:length(del.hypo)){
                    del.grp <- annoID.pos[del.hypo[iii], 2];
                    index2.pos <- which(annoID.pos[,2] == del.grp & annoID.pos[, 3] == id.h.pos);
                    if(length(index2.pos) > 0){
                      tmp.pos <- annoID.pos[index2.pos, 1];
                      if(!tmp.pos %in% anno.ids.pos){
                        grp2del <- append(grp2del, del.grp);
                      }
                    }else{
                      grp2del <- append(grp2del, del.grp);
                    }
                  }
                }
              }else{
                #adduct annotation is not M+H
                grp.hyp.pos <- annoID.pos[which(annoID.pos[, 1] == id.pos), 2];
                index.pfad <- annoID.pos[which(annoID.pos[, 2] %in% grp.hyp.pos & annoID.pos[, 3] == id.h.pos), 1];
                if(any(index.pfad %in% anno.ids.pos)) { #Pfad: 3.2
                  #Pfad: 3.2.1
                  next;
                }else{
                  #Pfad: 3.2.2
                  ##Test auf Grp.
                  if(is.null(masslist.neg)){
                      next;
                  }
                  if(length(unlist(fastMatch(test.mass <- annoGrp.pos[annoID.pos[which(annoID.pos[, 1] == id.pos), 2], "mass"], masslist.neg, tol=0.05))) > 0){
                    endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 4));
                  }else{
                    del.hypo <- which(annoID.pos[, 1] == id.pos & annoID.pos[, 3] != id.h.pos);
                    del.grp <- annoID.pos[del.hypo, 2];
                    grp2del <- append(grp2del, del.grp);
                    endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 3));
                  }
                }
              }
            }else{
              #2.1
              #add 1 to endresult, means found new M+H
              endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 1));
            }
          }else if(!pos){
            #Join Results to negative List
            id.neg <- xsa.neg@pspectra[[m[[i]][j]]][results[[ii]][2]];
            if(length(xsa.neg@derivativeIons[[id.neg]]) > 0){ #Pfad: 2
              if(length(index.neg <- which(annoID.neg[, 1] == id.neg & annoID.neg[, 3] == id.h.neg)) > 0){ #Pfad: 3
                #3.1
                if(length(index.neg) > 1){
                  cat("index.pos are greater than allowed. please debug!");
                }
                if(!length(del.hypo <- which(annoID.neg[,1] == id.neg & annoID.neg[, 3] != id.h.neg)) > 0) { #Pfad 3.1
                  #Pfad: 3.1.1
                  #bestätige Hypothese
                  endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 2));
                }else{
                  #Pfad: 3.1.2
                  #delete Hypothese
                  ##TODO: DEL
                  endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 2));
                  for(iii in 1:length(del.hypo)){
                    del.grp <- annoID.neg[del.hypo[iii], 2];
                    index2.neg <- which(annoID.neg[, 2] == del.grp & annoID.neg[, 3] == id.h.neg)
                    if(length(index2.neg) > 0){
                      tmp.neg <- annoID.neg[index2.neg, 1];
                      if(!tmp.neg %in% anno.ids.neg){
                        grp2del <- append(grp2del, del.grp);
                      }
                    }else{
                      grp2del <- append(grp2del, del.grp);
                    }
                  }
                }
              }else{
                grp.hyp.neg <- annoID.neg[which(annoID.neg[, 1] == id.neg), 2];
                index.pfad3 <- annoID.neg[which(annoID.neg[, 2] %in% grp.hyp.neg & annoID.neg[, 3] == id.h.neg), 1];
                if(any(index.pfad3 %in% anno.ids.neg)) { #Pfad: 3.2
                  #Pfad: 3.2.1
                  next;
                }else{
                  #Pfad: 3.2.2
                   ##Test auf Grp.
                  if(is.null(masslist.pos)){
                      next;
                  }
                  ##Test auf Grp.
                  if(length(unlist(fastMatch(test.mass <- annoGrp.neg[annoID.neg[which(annoID.neg[, 1] == id.neg), 2], "mass"], masslist.pos, tol=0.05))) > 0){
                    endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 4));
                  }else{
                    del.hypo <- which(annoID.neg[, 1] == id.neg & annoID.neg[, 3] != id.h.neg);
                    del.grp <- annoID.neg[del.hypo, 2];
                    grp2del <- append(grp2del, del.grp);
                    endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 3));
                  }
                }
              }
            }else{
              #2.1
              endresult <- rbind(endresult, c(i, results[[ii]][1], m[[i]][j], results[[ii]][2], mass, 1));
            }
          }else{
            cat("Error: pos have wrong value. please debug!\n");
          }
        }
      }
    }
    cat("\n")
  }
  if(pos){
    index<-which(annoID.pos[,2] %in% grp2del);
    if(length(index)>0){
      annoID.pos<-annoID.pos[-index,];
    }
    add.adducts<-vector("character",length(xsa.pos@isotopes));
    old.grpid<-max(annoGrp.pos[,1]);
    if(nrow(endresult)>0){
      for(i in 1:nrow(endresult)){
        if(!endresult[i,"check"] %in% c(2,4)){
          old.grpid<-old.grpid+1;
          peakid<-xsa.pos@pspectra[[endresult[i,"grp_pos"]]][endresult[i,"peak_pos"]];
          annoID.pos<-rbind(annoID.pos,c(peakid,old.grpid,id.h.pos))
          annoGrp.pos<-rbind(annoGrp.pos,c(old.grpid,endresult[i,"mass"],2))
          if(endresult[i,"check"]==1){
            add.adducts[peakid]<-paste("Found [M-H]");
          }else{
            add.adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");
          }
        }else{
          peakid<-xsa.pos@pspectra[[endresult[i,"grp_pos"]]][endresult[i,"peak_pos"]];
          add.adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");
        }
      }
    }
    xsa.pos@derivativeIons<-getderivativeIons(annoID.pos,annoGrp.pos,rules.pos,length(xsa.pos@isotopes))
    peaklist<-getPeaklist(xsa.pos);
    index<-ncol(peaklist)
#     if(xsa.pos@sample>0){
      #grouped xsa
#       new.adducts<-vector("character",length(nrow(peaklist)));
#       new.adducts[xsa.pos@grp_info]<-add.adducts;
#       peaklist<-cbind(peaklist,new.adducts);
#     }else{
      peaklist<-cbind(peaklist,add.adducts);
#     }
    colnames(peaklist)[index+1]<-"neg. Mode"
  }else if(!pos){
  index<-which(annoID.neg[,2] %in% grp2del);
  if(length(index)>0){
    annoID.neg<-annoID.neg[-index,];
  }
  add.adducts<-vector("character",length(xsa.neg@isotopes));
  old.grpid<-max(annoGrp.neg[,1]);
  if(nrow(endresult)>0){
    for(i in 1:nrow(endresult)){
      if(!endresult[i,"check"] %in% c(2,4)){
        old.grpid<-old.grpid+1;
        peakid<-xsa.neg@pspectra[[endresult[i,"grp_neg"]]][endresult[i,"peak_neg"]];
        annoID.neg<-rbind(annoID.neg,c(peakid,old.grpid,id.h.neg))
        annoGrp.neg<-rbind(annoGrp.neg,c(old.grpid,endresult[i,"mass"],2))
        if(endresult[i,"check"]==1){
          add.adducts[peakid]<-paste("Found [M+H]");
        }else{
          add.adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");
        }
      }else{
        peakid<-xsa.neg@pspectra[[endresult[i,"grp_neg"]]][endresult[i,"peak_neg"]];
        add.adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");
      }
    }
  }
  xsa.neg@derivativeIons<-getderivativeIons(annoID.neg,annoGrp.neg,rules.neg,length(xsa.neg@isotopes))
  peaklist<-getPeaklist(xsa.neg);
  index<-ncol(peaklist)
#   if(xsa.neg@sample>0){
      #grouped xsa
#       new.adducts<-vector("character",length(nrow(peaklist)));
#       new.adducts[xsa.neg@grp_info]<-add.adducts;
#       peaklist<-cbind(peaklist,new.adducts);
#     }else{
      peaklist<-cbind(peaklist,add.adducts);
#     }
  colnames(peaklist)[index+1]<-"pos. Mode"
}else return(NULL);
return(peaklist);
}

combine_hypothese <- function(m.pos,m.neg,rule_hh,ini1,ini2,tol=0.02){
    ##generiere neue Peaklist
    ML <- massDiffMatrix(m.pos,rule_hh)
    m <- fastMatch(ML,m.neg,tol=tol)
    results<-list();
    if(length(m)==0) {return(results);}
    for(i in 1:length(m))
    {
        if(is.null(m[[i]]))next;
        results[[length(results)+1]]<-c(ini1[i],ini2[m[[i]]])
    }
    return(results);
}

naOmit <- function(x) {
    return (x[!is.na(x)]);
}

calc_pc <-function(object,CL,cor_matrix,psg_list=NULL,psSamples=NULL) {
  
  require(RBGL)
  pspectra<-object@pspectra;
#   gm <- matrix(-1,1,2); colnames(gm) <- c('fromID','toID')
  li <- sapply(CL,function(x) length(x) > 0); #l <- which(li)
  if (!any(li)) { return(object@pspectra) }
  npspectra <- length(object@pspectra);
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##Ein Sample Fall
    imz <- object@xcmsSet@peaks[,"mz"];
  }else {
    ##Mehrsample Fall
    imz <- groups(object@xcmsSet)[,"mzmed"]
  }
  if(is.null(psg_list)){
    cat('\nCalculating graph cross linking in',npspectra,'Groups... \n % finished: '); lp <- -1;
    pspectra_list<-1:npspectra;
    ncl<-sum(sapply(object@pspectra,length));
  }else{
    cat('\nCalculating graph cross linking in',length(psg_list),'Groups... \n % finished: '); lp <- -1;
    pspectra_list<-psg_list;
    ncl<-sum(sapply(object@pspectra[psg_list],length));
  }

  npeaks = 0;
  if(nrow(object@isoID)){
    #Isotope wurden vorher erkannt
    
    idx<-unique(object@isoID[,1]); #ID aller monoisotopischen Peaks
  }else {idx<-NULL};
  for(j in 1:length(pspectra_list)){
    i <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    npeaks <- npeaks + length(pi);
    perc <- round(npeaks / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    if (.Platform$OS.type == "windows") flush.console()
 
    ## create list of connected components
    V <- as.character(pi);
    OG <- new("graphNEL",nodes=V)
    edgemode(OG) <- "undirected"
    gm <- matrix(NA,ncol=2);
    ow <- options("warn");
    options(warn = -1);
    for(k in 1:length(pi)){
      to <- CL[[pi[k]]];
      if(length(to)>0) {
        OG <- addEdge(from=as.character(pi[k]),to=as.character(to),graph=OG,weights=1);
      }
    }
    options(ow);
    NG <- matrix(NA,ncol=2);
    ## verify all correlation graphs
    G <- OG;
    if (length(nodes(G)) > 2) {
      ## decomposition might be necessary
      # G <- removeSelfLoops(G)
      hcs <-  highlyConnSG(G)
      lsg <- sapply(hcs$clusters,function(x) length(x))
      lsg.i <- which(lsg > 1)
      if (length(lsg.i)<1) next;
      for(z in 1:length(hcs$clusters)){
        NG<-rbind(NG,cbind(z,as.numeric(hcs$clusters[[z]])))
      }
      NG<-NG[-1,];#Remove NA
      if(object@polarity %in% c("positive","negative")){
        if(object@polarity == "positive"){
          ##Hold M+H,M+Na
          rules<-data.frame(c("[M+H-M+Na]","[M+H-M+K]","[M+Na-M+Na]"),1,1,c(21.9812,37.9552,15.974),1,1,1)
          colnames(rules)<-c("name","nmol","charge","massdiff","oidscore","quasi","ips");         
        }else {
          #negative Way
          ##Hold M-H,M-Na
          rules<-data.frame(c("[M-H-M-2H+Na]","[M-H-M+Cl]"),1,1,c(21.9812,35.9758),1,1,1)
          colnames(rules)<-c("name","nmol","charge","massdiff","oidscore","quasi","ips");         
        }
        mz<-imz[pi];ix<-order(mz);mz<-mz[ix]
        mm<-matrix(NA,ncol=2)
        for(x in 1:length(mz)){
          for(y in x:length(mz)){
            diff<-mz[y]-mz[x];
            if(diff>39){break} ##Diff zu gross
            if(length(unlist(fastMatch(rules[,"massdiff"],diff,0.02)))>0){
              #One rule fit, not of interest which one
              mm<-rbind(mm,c(pi[ix[x]],pi[ix[y]]));
            }
          }
        }
        if(nrow(mm)>1){
          mm<-mm[-1,]#remove NA        
          mm<-matrix(mm,ncol=2);
          for(x in 1:nrow(mm)){
            if(!mm[x,1] %in% NG[,2]){
                  NG<-rbind(NG,c(max(NG[,1])+1,mm[x,1]));
            }
            grp<-NG[which(NG[,2]==mm[x,1]),1];
            if(!mm[x,2] %in% NG[,2]){
                  NG<-rbind(NG,c(grp,mm[x,2]));
            }else{
              NG[which(NG[,2]==mm[x,2]),1]<-grp;
            }
          }
        }
      }
      ##Hold Isotope together
      if(length(idx)>0){
        iidx<-which(pi %in% idx);
        if(length(iidx)>0){
          #Monoiso. in grp gefunden
          for(h in 1:length(iidx)){
            mindex<-pi[iidx[h]];
            if(!mindex %in% NG[,2]){
                NG<-rbind(NG,c(max(NG[,1])+1,mindex));
            }
            isoindex<-object@isoID[which(object@isoID[,1]==mindex),2];
            for(t in 1:length(isoindex)){
              #if hcs sorted isopeak out
              if(!isoindex[t] %in% NG[,2]) {NG<-rbind(NG,c(0,isoindex[t]))}
            }
            grp<-NG[which(mindex==NG[,2]),1]
            NG[sapply(isoindex,function(x,NG) {which(x==NG[,2])},NG),1]<-grp;#Setze isotope auf gleichen index wie monoiso.
          }
        }
      }

      ## calculate all new pspectra
      grps<-unique(NG[,1]);
      cnts<-unlist(lapply(grps,function(x,NG) { length( which( NG[,1] == x) ) },NG))
      grps<-grps[order(cnts,decreasing = TRUE)]
      for (ii in 1:length(grps)){
        if(ii==1){
          #behalten alte Nummer
          pspectra[[i]] <- sort(NG[which(NG[,1]==grps[ii]),2]);

        } else {
          pspectra[[length(pspectra)+1]] <- sort(NG[which(NG[,1]==grps[ii]),2]);
          psSamples[length(psSamples)+1] <- psSamples[i] 
        }
      }
    } else {
      #Only one peak in the pseudospectra
#       pspectra[[i]] <- pi;
    }
  }
  ##Workarround: peaks without groups
    peaks<-vector("logical",nrow(object@groupInfo))
    npspectra<-length(pspectra)
    for(i in 1:npspectra){
        peaks[pspectra[[i]]]<-TRUE;
    }
    index<-which(peaks==FALSE);
    if(length(index)>0){
      for(i in 1:length(index)){
        pspectra[npspectra+i]<-index[i];
      }
    }
  object@pspectra<-pspectra;
  object@psSamples<-psSamples;
  cat("\n");
  return(object)
}

getAllEICs <- function(xs,file=NULL) {
  ##old CAMERA
  ##index = sample
#   peaki <- getPeaksIdxCol(xs,NULL);
#   if(is.matrix(peaki)) peaki<-peaki[,index]
#   scantimes <- list()
#   maxscans <- 0
#   if (file.exists(xs@filepaths[index])) {
#     cat('Reading raw data file:',xs@filepaths[index],'\n')
#     xraw <- xcmsRaw(xs@filepaths[index],profstep=0)
#     cat('Generating EIC\'s .. \n')
#     maxscans <- length(xraw@scantime)
#     scantimes[[1]] <- xraw@scantime
#     pdata <- as.data.frame(xs@peaks[peaki,])
#     EIC <- array(integer(0),c(nrow(pdata),maxscans,1))
#     EIC[,,1] <- getEICs(xraw,pdata,maxscans)
#   }else stop('Raw data file:',xs@filepaths[index],' not found ! \n')
#   invisible(list(scantimes=scantimes,EIC=EIC));

  ##new CAMERA
  peaki <- getPeaksIdxCol(xs,NULL)
  nfiles <- length(filepaths(xs))
  scantimes <- list()
  maxscans <- 0
  cat('Generating EIC\'s .. \n') 
  if (nfiles > 1) { 
      # cat('Searching maxima .. \n')
      for (f in 1:nfiles){
      #  cat('Reading raw data file:',filepaths(xs)[f]) 
        xraw <- xcmsRaw(filepaths(xs)[f],profstep=0)
#         cat(',', length(xraw@scantime),'scans. \n') 
        maxscans <- max(maxscans,length(xraw@scantime))
        scantimes[[f]] <- xraw@scantime
      }
  
      for (f in 1:nfiles){
        if (file.exists(filepaths(xs)[f])) { 
      #    cat('Reading raw data file:',filepaths(xs)[f],'\n') 
          xraw <- xcmsRaw(filepaths(xs)[f],profstep=0)
      #    cat('Generating EIC\'s .. \n') 
          pdata <- as.data.frame(xs@peaks[peaki[,f],]) # data for peaks from file f
          if (f==1) EIC <- array(integer(0),c(nrow(pdata),maxscans,length(filepaths(xs))))   
          EIC[,,f] <- getEICs(xraw,pdata,maxscans)
        }
        else stop('Raw data file:',filepaths(xs)[f],' not found ! \n')
      }
  }  else { ## create EIC's for single file
       if (file.exists(filepaths(xs)[1])) { 
         #cat('Reading raw data file:',filepaths(xs)[1],'\n') 
          xraw <- xcmsRaw(filepaths(xs)[1],profstep=0)
         #cat('Generating EIC\'s .. \n') 
          maxscans <- length(xraw@scantime)
          scantimes[[1]] <- xraw@scantime
          pdata <- as.data.frame(xs@peaks[peaki,]) 
          EIC <- array(integer(0),c(nrow(pdata),maxscans,1))   
          EIC[,,1] <- getEICs(xraw,pdata,maxscans)
        }  else stop('Raw data file:',filepaths(xs)[f],' not found ! \n') 
  } 
   
  if (!is.null(file)) save(EIC,scantimes,file=file,compress=TRUE)
    else invisible(list(scantimes=scantimes,EIC=EIC)) 
}

getPeaksIdxCol <- function(xs, col=NULL) {
if (nrow(xs@groups) > 0 && length(sampnames(xs)) > 1 )
    if (is.null(col)) m <- groupval(xs, "medret", value='index')
        else m <- groupval(xs, "medret", value='index')[,col]
else if (length(sampnames(xs)) == 1)
        m <- 1:length(xs@peaks[,'into'])
    else stop ('First argument must be a xcmsSet with group information or contain only one sample.')
m
}

getEICs <- function(xraw,peaks,maxscans=length(xraw@scantime)) {
npeaks <- dim(peaks)[1]; scans <- length(xraw@scantime)
eics <- matrix(as.numeric(0),npeaks,maxscans)
for (p in 1:npeaks) {
    eics[p,1:scans] <- as.integer(getEIC(xraw,massrange=c(peaks[p,"mzmin"],peaks[p,"mzmax"]))$intensity)
}
eics
}

getEIC <- function(xraw,massrange = numeric(), timerange = numeric(),scanrange= numeric() ){
if (!is.double(xraw@env$mz) || !is.double(xraw@env$intensity) || !is.integer(xraw@scanindex)) stop('mz/int not double.')
if (length(timerange) >= 2)
{
    timerange <- range(timerange)
    tidx <- which((xraw@scantime >= timerange[1]) & (xraw@scantime <= timerange[2]))
    scanrange <- range(tidx)
}else if (length(scanrange) < 2)
{
    scanrange <- c(1, length(xraw@scantime))
}else scanrange <- range(scanrange)

.Call("getEIC",xraw@env$mz,xraw@env$intensity,xraw@scanindex,as.double(massrange),as.integer(scanrange),
as.integer(length(xraw@scantime)), PACKAGE ='xcms' )

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

calcIsotopes <- function(maxiso,maxcharge){
        M_C12_C13 = 1.0033

        ## Calculate ISO/CHARGE - Matrix
        IM <- matrix(NA,maxiso,maxcharge)
        for (i in 1:maxiso)
        for (j in 1:maxcharge) IM[i,j] <- (i* M_C12_C13) / j

        return(IM)

}

massDiffMatrix <- function(m,adducts){
nadd <- length(adducts[,"massdiff"])
DM <- matrix(NA,length(m),nadd)

for (i in 1:length(m))
    for (j in 1:nadd)
    DM[i,j] <- (abs(adducts[j,"charge"] * m[i]) - adducts[j,"massdiff"]) / adducts[j,"nmol"]    # ((z*m) - add) /n

return(DM)
}

massDiffMatrixNL <- function(m,neutralloss){
nadd <- nrow(neutralloss)
DM <- matrix(NA,length(m),nadd)

for (i in 1:length(m))
    for (j in 1:nadd)
    DM[i,j] <- m[i] - neutralloss[j,"massdiff"]

return(DM)
}

calcCL2 <- function(object,xs, EIC, scantimes, cor_eic_th,nSlaves=2){
  CL <- vector("list",nrow(object@groupInfo));
  CIL <- list();
  ncl<-length(CL);npeaks=0;
  npspectra <- length(object@pspectra);
  peaks<-object@groupInfo;
  cat('Calculating peak correlations... \n');
  #Wenn groupFWHM nicht vorher aufgerufen wurde!
  if(npspectra<1){
    npspectra<-1;object@pspectra[[1]]<-seq(1:nrow(object@groupInfo));
  }
  cormat<-matrix(,ncol=2)
  for(i in 1:npspectra){
          pi <- object@pspectra[[i]];
          cormat<-rbind(cormat,as.matrix(expand.grid(pi,pi)));
  }
  index<-which(cormat[,1]>=cormat[,2]);
  cormat<-cormat[-index,]
  cormat<-cormat[-1,]
  options(show.error.messages = FALSE);
 
   index<-seq(1:nrow(cormat));
   liste<-lapply(index,function (x) {as.numeric(c(cormat[x,1],cormat[x,2]))})
   res <- mpi.parSapply(liste,function(x,EIC=EIC,peaks=peaks,scantimes=scantimes) {
    eicx <-  EIC[x[1],,1];eicy <-  EIC[x[2],,1];
    px <- peaks[x[1],];py <- peaks[x[2],];
    crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"]);
    rti <- which(scantimes[[1]] >= crt[1] & scantimes[[1]] <= crt[2])
    cors <- 0;
    if (length(rti)>1){
      dx <- eicx[rti]; dy <- eicy[rti]
      dx[dx==0] <- NA; dy[dy==0] <- NA;
      if (length(which(!is.na(dx) & !is.na(dy))) >= 4){
      ct <- NULL;
      try(ct <- cor.test(dx,dy,method='pearson',use='complete'));
      if(!is.null(ct) && !is.na(ct)){if(ct$p.value <= 0.05){cors <- ct$estimate}else{ cors <- 0;}}else cors <- 0;                                                            
      }
    }
  return(cors);
  },EIC,peaks,scantimes)
  options(show.error.messages = TRUE);
  index<-as.numeric(which(res>cor_eic_th))
  sapply(index,function(x) {
      CL[[cormat[x,1]]]<<-c(CL[[cormat[x,1]]],as.numeric(cormat[x,2]));
      CL[[cormat[x,2]]]<<-c(CL[[cormat[x,2]]],as.numeric(cormat[x,1]));
      invisible(CIL[[length(CIL)+1]]<<-list(p=c(cormat[x,1],cormat[x,2]),cor=as.numeric(res[x])));
  })
  if (length(CIL) >0){ CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) )
  }else{ return(NULL)}
  colnames(CI) <- c('xi','yi','cors');
  return(invisible(list(CL=CL,CI=CI)));
}

calcCL <-function(object, EIC, scantimes, cor_eic_th, psg_list=NULL){
  xs <- object@xcmsSet;
  if(is.na(object@sample)){
    peaki <- getPeaksIdxCol(xs,col=NULL)
    peaks <- groupval(xs,value="maxo")
  }else if(object@sample == -1){
    ##TODO @Joe: Sollte das hier auftreten?
  }else{
    peaki <- getPeaksIdxCol(xs,col=object@sample)
  }
  Nf <- length(filepaths(xs))
  if(is.vector(peaki)){ peaki <- as.matrix(peaki) }
  Nrow <- nrow(peaki)
  CL <- vector("list",Nrow)
  CIL <- list()
  ncl<-length(CL);
  
  npeaks=0;
  npspectra <- length(object@pspectra);
  
  #Wenn groupFWHM nicht vorher aufgerufen wurde!
  if(npspectra<1){
    npspectra<-1;object@pspectra[[1]]<-seq(1:nrow(object@groupInfo));
    cat('Calculating peak correlations for 1 big group.\nTry groupFWHM before, to reduce runtime. \n% finished: '); lp <- -1;
    pspectra_list<-1;
    object@psSamples<-1;
  }else{
    if(is.null(psg_list)){
      cat('\nCalculating peak correlations in',npspectra,'Groups... \n % finished: '); lp <- -1;
      pspectra_list<-1:npspectra;
    }else{
      cat('\nCalculating peak correlations in',length(psg_list),'Groups... \n % finished: '); lp <- -1;
      pspectra_list<-psg_list;
      ncl<-sum(sapply(object@pspectra[psg_list],length));
    }
  }
  psSamples <- object@psSamples;
  for(j in 1:length(pspectra_list)){
    i <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    #percent output
    npeaks<-npeaks+length(pi);
    perc <- round((npeaks) / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    if (.Platform$OS.type == "windows") flush.console()
    #end percent output

    #select sample f
    if(is.na(object@sample)){
      if(length(pi)>1){
        f <- as.numeric(which.max(apply(peaks[pi,],2,function(x){mean(x,na.rm=TRUE)}))) #errechne höchsten Peaks, oder als mean,median
        psSamples[i] <- f;
      }else{
        f <- which.max(peaks[pi,]);
        psSamples[i] <- f;
      }
    }else {
        f<-object@sample;
        psSamples[i] <- f;
    }
    #end selection

    if(length(pi)>1){
      for(x in 2:length(pi)){
        xi <- pi[x];pxi<-peaki[xi,f];
        for (y in 1:(x-1)){
          yi <- pi[y];pyi<-peaki[yi,f];
          if ( ! (yi %in% CL[[xi]] || yi == xi)){
            cors <-0;
            eicx <-  EIC[xi,,f]
            eicy <-  EIC[yi,,f]
            px <- xs@peaks[pxi,]
            py <- xs@peaks[pyi,]
            crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"])
            rti <- which(scantimes[[1]] >= crt[1] & scantimes[[1]] <= crt[2])
            if (length(rti)>1){
              dx <- eicx[rti]; dy <- eicy[rti]
              dx[dx==0] <- NA; dy[dy==0] <- NA;
              if (length(which(!is.na(dx) & !is.na(dy))) >= 4){
                ct <- NULL
                options(show.error.messages = FALSE)
                try(ct <- cor.test(dx,dy,method='pearson',use='complete'))
                options(show.error.messages = TRUE)
                if (!is.null(ct) && !is.na(ct)){
                  if (ct$p.value <= 0.05) {
                    cors <- ct$estimate;
                  }  else { cors <- 0; }
                } else { cors <- 0; }
              } else { cors <- 0; }
            } else { cors <- 0; }
            if(cors>cor_eic_th){
              CL[[xi]] <- c(CL[[xi]],yi)
              CL[[yi]] <- c(CL[[yi]],xi) ## keep the list symmetric
              CIL[[length(CIL)+1]] <- list(p=c(xi,yi),cor=cors)
            }
          }
        }
      }
    }
  }
  cat("\n");
  if (length(CIL) >0) {
    CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) );
  } else { return(NULL) }
  colnames(CI) <- c('xi','yi','cors')
  return(invisible(list(CL=CL,CI=CI,psSamples=psSamples)))
}
###End xsAnnotate intern Function###
