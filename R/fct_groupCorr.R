##Functions for groupCorr

calc_pc <-function(object,CL,psg_list=NULL,psSamples=NULL) {
  
  pspectra <- object@pspectra;
  li <- sapply(CL, function(x) length(x) > 0);
  if (!any(li)){
    #No Connection Found
    return(object@pspectra);
  }

  npspectra <- length(object@pspectra);
  if(object@sample == 1 && length(sampnames(object@xcmsSet)) == 1){
    ##Ein Sample Fall
    imz <- object@xcmsSet@peaks[, "mz"];
  }else {
    ##Mehrsample Fall
    imz <- groups(object@xcmsSet)[,"mzmed"]
  }

  if(is.null(psg_list)){
    cat('\nCalculating graph cross linking in', npspectra, 'Groups... \n % finished: ');
    lperc <- -1;
    pspectra_list <- 1:npspectra;
    ncl <- sum(sapply(object@pspectra, length));
  }else{
    cat('\nCalculating graph cross linking in',length(psg_list),'Groups... \n % finished: ');
    lperc <- -1;
    pspectra_list <- psg_list;
    ncl <- sum(sapply(object@pspectra[psg_list], length));
  }

  npeaks <- 0;
  if(nrow(object@isoID)){
    #Isotope wurden vorher erkannt
    idx <- unique(object@isoID[, 1]); #ID aller monoisotopischen Peaks
  }else{
    idx<-NULL
  };
  for(j in 1:length(pspectra_list)){
    i  <- pspectra_list[j];#index of pseudospectrum
    pi <- object@pspectra[[i]]; #peak_id in pseudospectrum
    
    ## % output
    npeaks <- npeaks + length(pi); 
    perc <- round(npeaks / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lperc)) { 
      cat(perc, ' ');
      lperc <- perc; 
    }
    if (.Platform$OS.type == "windows"){
      flush.console();
    }
  
    ## create list of connected components
    V  <- as.character(pi);
    
#     time <- proc.time()
#     OG <- new("graphNEL", nodes=V)
#     edgemode(OG) <- "undirected";
#     ow <- options("warn");
#     options(warn = -1);
#     for(k in 1:length(pi)){
#       to <- CL[[pi[k]]];
#       if(length(to)>0) {
#         OG <- addEdge(from=as.character(pi[k]), to=as.character(to), graph=OG, weights=1);
#       }
#     }
#     options(ow);
#     proc.time() - time

#     time <- proc.time()
    L <- matrix(NA,ncol=2,nrow=0);
    colnames(L) <- c("From","To");
    invisible( sapply(1:length(pi), function(x) { 
      to <- which(CL[[pi[x]]] > pi[x]);
      if(length(to) > 0){ 
        L <<- rbind(L,cbind(as.character(pi[x]),as.character(CL[[pi[x]]][to])))
      }
      return(NULL)}))

    W <- rep(1,nrow(L));
    OG <- ftM2graphNEL(L, W, edgemode="undirected");
    index <- which(!V %in% OG@nodes); #Has graph all nodes?
    if(length(index) > 0){ #graph misses nodes
      OG <- addNode(as.character(V[index]),OG)
    }
#     proc.time() - time
    rm(L,W);

    
    NG <- matrix(NA, ncol=2);
    ## verify all correlation graphs
    if (length(nodes(OG)) > 2) {
      ## decomposition might be necessary
      hcs <- highlyConnSG(OG);
      rm(OG);gc();
      #calculate size of hcs cluster
      lsg <- sapply(hcs$clusters, function(x) length(x));
      lsg.i <- which(lsg > 1)
      #Have we at least one hcs with length > 1
      if (length(lsg.i) < 1){
        next;
      }
      #NG[index of cluster, index of peak]
      for(z in 1:length(hcs$clusters)){
        NG <- rbind(NG, cbind(z, as.numeric(hcs$clusters[[z]])))
      }
      NG <- NG[-1,];#Remove NA
      rm(hcs);gc();

      ##Hold primary adducts together
      if(object@polarity %in% c("positive","negative")){
        if(object@polarity == "positive"){
          ##Hold M+H,M+Na
          rules <- data.frame(c("[M+H-M+Na]","[M+H-M+K]","[M+Na-M+Na]"),1,1,c(21.9812,37.9552,15.974),1,1,1)
          colnames(rules) <- c("name","nmol","charge","massdiff","oidscore","quasi","ips");         
        }else {
          #negative Way
          ##Hold M-H,M-Na
          rules <- data.frame(c("[M-H-M-2H+Na]","[M-H-M+Cl]"),1,1,c(21.9812,35.9758),1,1,1)
          colnames(rules) <- c("name","nmol","charge","massdiff","oidscore","quasi","ips");         
        }
        mz <- imz[pi];
        ix <- order(mz);
        mz <- mz[ix];
        mm <- matrix(NA, ncol=2);
        for(x in 1:length(mz)){
          for(y in x:length(mz)){
            diff <- mz[y]-mz[x];
            if(diff > 39){
              break; ##Diff zu gross
            }
            if(length(unlist(fastMatch(rules[, "massdiff"], diff, 0.02))) > 0){
              #One rule fit, not of interest which one
              mm <- rbind(mm, c(pi[ix[x]], pi[ix[y]]));
            }
          }
        }
        if(nrow(mm) > 1){
          mm <- mm[-1, ]; #remove NA        
          mm <- matrix(mm, ncol=2);
          for(x in 1:nrow(mm)){
            if(!mm[x,1] %in% NG[,2]){
              NG <- rbind(NG, c(max(NG[, 1])+1, mm[x, 1]));
            }
            grp <- NG[which(NG[, 2]==mm[x, 1]), 1];
            if(!mm[x,2] %in% NG[,2]){
              NG <- rbind(NG, c(grp, mm[x, 2]));
            }else{
              NG[which(NG[, 2] == mm[x, 2]), 1] <- grp;
            }
          }
        }
      }

      ##Hold Isotope together
      if(length(idx) > 0){
        iidx <- which(pi %in% idx);
        if(length(iidx) > 0){
          #Monoiso. in grp gefunden
          for(h in 1:length(iidx)){
            mindex <- pi[iidx[h]];
            if(!mindex %in% NG[, 2]){
                NG <- rbind(NG, c(max(NG[, 1])+1, mindex));
            }
            isoindex <- object@isoID[which(object@isoID[, 1] == mindex), 2];
            for(t in 1:length(isoindex)){
              #if hcs sorted isopeak out
              if(!isoindex[t] %in% NG[, 2]) {
                NG <- rbind(NG, c(0, isoindex[t]))
              }
            }
            grp <- NG[which(mindex == NG[, 2]), 1];
            NG[sapply(isoindex, function(x, NG) {which(x == NG[, 2])}, NG), 1] <- grp;#Setze isotope auf gleichen index wie monoiso.
          }
        }
      }
     ## calculate all new pspectra
      grps <- unique(NG[, 1]);
      cnts <- unlist(lapply(grps, function(x, NG) { length( which( NG[,1] == x) ) }, NG))
      grps <- grps[order(cnts, decreasing = TRUE)]
      pcollect<-NULL ## collecting all indicies used for old/new pspectra
      for (ii in 1:length(grps)){
        if(ii==1){
          #behalten alte Nummer
          pspectra[[i]] <- sort(NG[which(NG[, 1] == grps[ii]), 2]);
          pcollect <- c(pcollect, pspectra[[i]]);
        } else {
          lp <- length(pspectra)+1
          pspectra[[lp]] <- sort(NG[which(NG[,1]==grps[ii]),2]);
          psSamples[lp] <- psSamples[i]
          pcollect <- c(pcollect, pspectra[[lp]])
        }
      }
      singleton <- pi[which(!(pi %in% pcollect))]
      if (length(singleton)>0){
        npspectra <- length(pspectra)
        for(i in 1:length(singleton)){ ## inserting singleton peaks in own pspectra
          pspectra[npspectra+i] <- singleton[i];
          psSamples[npspectra+i] <- psSamples[j]
        }
      }
    } else {
      #Only one peak in the pseudospectra
#       pspectra[[i]] <- pi;
        rm(OG);gc();
    }
    ##rm not used data
    rm(NG);gc();
}
#   ##Workarround: peaks without groups
#     peaks<-vector("logical",nrow(object@groupInfo))
#     npspectra<-length(pspectra)
#     for(i in 1:npspectra){
#         peaks[pspectra[[i]]] <- TRUE;
#     }
#     index <- which(peaks==FALSE);
#     if(length(index) > 0){
#       for(i in 1:length(index)){
#         pspectra[npspectra+i] <- index[i];
#         psSamples[npspectra+i] <- object@groupInfo[index[i], "sample"];
#       }
#     }
  object@pspectra <- pspectra;
  object@psSamples <- psSamples;
  cat("\n");
  return(object)
}

calcCL3 <- function(object, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th, psg_list=psg_list){
  xs <- object@xcmsSet;

  if(is.na(object@sample)){
    peaki <- getPeaksIdxCol(xs,col=NULL)
#     peaks <- groupval(xs,value="maxo")
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
  CI <- NULL;
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
    f <- psSamples[j]; #?????
    npi <- length(pi);
    #percent output
    npeaks<-npeaks+length(pi);
    perc <- round((npeaks) / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    if (.Platform$OS.type == "windows") flush.console()
    #end percent output
    if(npi < 2) {
      next;
    }
    options(warn = -1);
    res <- rcorr(EIC[,pi],type="pearson")
    options(warn = 0);
    res$r[lower.tri(res$r,diag = TRUE)] <- NA
    res$P[lower.tri(res$P,diag = TRUE)] <- NA
    index <- which(( res$r > cor_eic_th) & (res$P <= 0.05))
    if(length(index) > 0){
     for( x in 1:(length(index))){
      col <- index[x] %/% npi + 1;
      row <- index[x] %%  npi;
      if(row == 0){
        row <- npi;
        col <- col - 1;
      };
      xi <- pi[row];yi <- pi[col];
      CL[[xi]] <- c(CL[[xi]],yi)
      CL[[yi]] <- c(CL[[yi]],xi) ## keep the list symmetric
      CIL[[length(CIL)+1]] <- list(p=c(xi,yi),cor=res$r[row,col])
      };
    }else{
      next;
    }
  }

  if (length(CIL) >0){ CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) )
  }else{ return(NULL)}
  colnames(CI) <- c('xi','yi','cors');
  cat("\n");
  return(invisible(list(CL=CL,CI=CI)))

}

calcCL2 <- function(object,EIC, scantimes, cor_eic_th,psg_list=NULL){
  xs <- object@xcmsSet;
  peaks <- xs@peaks;

  if(is.na(object@sample)){
    peaki <- getPeaksIdxCol(xs, col=NULL);
    #peaks <- groupval(xs, value="maxo");
  }else if(object@sample == -1){
    ##TODO @Joe: Sollte das hier auftreten?
  }else{
    peaki <- getPeaksIdxCol(xs, col=object@sample);
  };

  Nf <- length(filepaths(xs)) #Anzahl Samples

  if(is.vector(peaki)){ 
    peaki <- as.matrix(peaki);
  };

  Nrow <- nrow(peaki);
  
  CL  <- vector("list", nrow(object@groupInfo));
  CIL <- list();
  ncl <- length(CL);

  npeaks    <- 0;
  npspectra <- length(object@pspectra);

  #Wenn groupFWHM nicht vorher aufgerufen wurde!
  if(npspectra < 1){
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    cat('Calculating peak correlations for 1 big group.\nTry groupFWHM before, to reduce runtime. \n% finished: '); 
    lp <- -1;
    pspectra_list    <- 1;
    object@psSamples <- 1;
  }else{
    if(is.null(psg_list)){
      cat('\nCalculating peak correlations in',npspectra,'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating peak correlations in',length(psg_list),'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- psg_list;
      ncl <- sum(sapply(object@pspectra[psg_list], length));
    }
  }
  psSamples <- object@psSamples;

  cormat<-matrix(,ncol=3)
  for(j in 1:length(pspectra_list)){
    i <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    cormat<-rbind(cormat,cbind(as.matrix(expand.grid(pi,pi)),psSamples[j]));
  }

  index<-which(cormat[,1]>=cormat[,2]);
  cormat<-cormat[-index,]
  cormat<-cormat[-1,]
  colnames(cormat) <- c("x1","x2","x3");
  dimnames(cormat)[[1]] <- 1:nrow(cormat);

  options(show.error.messages = FALSE);
 
  for(j in 1:length(pspectra_list)){
    i <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    pindex <- peaki[pi,psSamples[j]];
    rttimes <- peaks[pindex,c("rtmin","rtmax")];
    apply(as.vector(rttimes),1, function(x){
       which(scantimes[[psSamples[j]]] >=x[1] & scantimes[[psSamples[j]]] <=x[1]);
       which(scantimes[[psSamples[j]]] >=x[2] & scantimes[[psSamples[j]]] <=x[2]);
    })
    rt.max <- max(rttimes[,2]);
    rt.min <- min(rttimes[,1]);
    
    rti <- which(scantimes[[psSamples[j]]] >=rt.min & scantimes[[psSamples[j]]] <=rt.max)
    
    pmat <- matrix(NA,nrow=length(pi),ncol=length(rti))
    
    apply(rttimes,1, function (x) {
       which(scantimes[[psSamples[j]]] >=x[1] & scantimes[[psSamples[j]]] <=x[2]);
    })
    rti <- which(scantimes[[psSamples[j]]] >=rttimes[,1] & scantimes[[psSamples[j]]] <=rttimes[,2])
    peak.eic <- EIC[pindex]
    scantimes[[psSamples[j]]][rttimes];
    
  }

   proc <- proc.time()
   res <- mpi.parApply(cormat,1,function(x,EIC=EIC,peaks=xs@peaks,scantimes=scantimes,peaki=peaki) {
    eicx <-  EIC[x[1],,1];eicy <-  EIC[x[2],,1];
    px <- peaks[peaki[x[1],x[3]],];py <- peaks[peaki[x[2],x[3]],];
    crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"]);
    rti <- which(scantimes[[x[3]]] >= crt[1] & scantimes[[x[3]]] <= crt[2])
    cors <- 0;
    if (length(rti)>1){
      dx <- eicx[rti]; dy <- eicy[rti]
      dx[dx==0] <- NA; dy[dy==0] <- NA;
      if (length(which(!is.na(dx) & !is.na(dy))) >= 4){
      ct <- NULL;
      try(ct <- cor.test(dx[-index],dy[-index],method='pearson',use='complete'));
      if(!is.null(ct) && !is.na(ct)){if(ct$p.value <= 0.05){cors <- ct$estimate}else{ cors <- 0;}}else cors <- 0;                                                            
      }
    }
  return(cors);
  },EIC,xs@peaks,scantimes,peaki)
  proc.time() - proc;
 
cor_function <- function(x,eicx,eicy,px,py,scantimes){
    crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"]);
    rti <- which(scantimes[[x[3]]] >= crt[1] & scantimes[[x[3]]] <= crt[2])
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
}
   proc <- proc.time()

it <- iter(as.vector(sapply(1:nrow(cormat),function(x) rep(x,6))))

res2 <- mpi.parApply(cormat,1, cor_function,eicx=EIC[cormat[nextElem(it),1],,1],eicy=EIC[cormat[nextElem(it),2],,1],px=peaks[peaki[cormat[nextElem(it),1],cormat[nextElem(it),3]],],py=peaks[peaki[cormat[nextElem(it),2],cormat[nextElem(it),3]],],scantimes=scantimes)
  proc.time() - proc;

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
#     peaks <- groupval(xs,value="maxo")
  }else if(object@sample == -1){
    ##TODO @Joe: Sollte das hier auftreten?
  }else{
    peaki <- getPeaksIdxCol(xs,col=object@sample)
  }
  Nf <- length(filepaths(xs))
  if(is.vector(peaki)){ peaki <- as.matrix(peaki) }
  Nrow <- nrow(peaki)
  CL <- vector("list",Nrow)
#   CIL <- list()
  CI <- NULL;
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
#     if(is.na(object@sample)){
#       if(length(pi)>1){
#         f <- as.numeric(which.max(apply(peaks[pi,],2,function(x){mean(x,na.rm=TRUE)}))) #errechne höchsten Peaks, oder als mean,median
#         psSamples[i] <- f;
#       }else{
#         f <- which.max(peaks[pi,]);
#         psSamples[i] <- f;
#       }
#     }else {
#         f<-object@sample;
#         psSamples[i] <- f;
#     }
    #end selection
    f <- psSamples[j];

    if(length(pi)>1){
      for(x in 2:length(pi)){
        xi <- pi[x];pxi<-peaki[xi,f];
        for (y in 1:(x-1)){
          yi <- pi[y];pyi<-peaki[yi,f];
          if ( ! (yi %in% CL[[xi]] || yi == xi)){
            cors <-  0;
            eicx <-  EIC[xi,,1]
            eicy <-  EIC[yi,,1]
            px <- xs@peaks[pxi,]
            py <- xs@peaks[pyi,]
            crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"])
            rti <- which(scantimes[[f]] >= crt[1] & scantimes[[f]] <= crt[2])
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
#                CIL[[length(CIL)+1]] <- list(p=c(xi,yi),cor=cors)
            }
          }
        }
      }
    }
  }
  cat("\n");
#   if (length(CIL) >0) {
#     CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) );
#   } else { return(NULL) }
#   colnames(CI) <- c('xi','yi','cors')
  return(invisible(list(CL=CL,CI=CI)))
}


getAllPeakEICs <- function(xs,index=NULL){
  peaki <- getPeaksIdxCol(xs,NULL)
  nfiles <- length(filepaths(xs))
  scantimes <- list()
  maxscans <- 0
  if(is.null(index)){
    cat("Missing Index, generate all EICs from sample 1.\n");
    index <- rep(1,nrow(peaki));
  }
  cat('Generating EIC\'s .. \n') 
  if (nfiles > 1) { 
    for (f in 1:nfiles){
      xraw <- xcmsRaw(filepaths(xs)[f], profstep=0);
      maxscans <- max(maxscans, length(xraw@scantime));
      scantimes[[f]] <- xraw@scantime;
    }
    EIC <- array(integer(0), c(nrow(peaki),maxscans))
    na.flag <- 0;
    for (f in 1:nfiles){
      if (file.exists(filepaths(xs)[f])) { 
        xraw <- xcmsRaw(filepaths(xs)[f], profstep=0);
        idx.peaks <- which(index == f);
        if(length(idx.peaks) > 0){
            pdata <- as.data.frame(xs@peaks[peaki[idx.peaks,f],]) # data for peaks from file f
            if(length(which(is.na(pdata[,1]))) >0){
              na.flag <- 1;
            }
            if(length(idx.peaks)==1){
              pdata <- t(pdata);
            }
            EIC[idx.peaks,] <- getEIC4Peaks(xraw,pdata,maxscans)
          }
        }
        else stop('Raw data file:',filepaths(xs)[f],' not found ! \n')
      }
    if(na.flag ==1){
      cat("Found NA peaks in selected samples. Those will be seperated in new pcgroups in each case.\nUse fillpeaks if not desired!\n");
    }
  }  else { ## create EIC's for single file
       if (file.exists(filepaths(xs)[1])) { 
          xraw <- xcmsRaw(filepaths(xs)[1],profstep=0)
          maxscans <- length(xraw@scantime)
          scantimes[[1]] <- xraw@scantime
          pdata <- as.data.frame(xs@peaks[peaki,]) 
          EIC <- array(NA,c(nrow(pdata),maxscans))   
          EIC[,] <- getEIC4Peaks(xraw,pdata,maxscans)
        }  else stop('Raw data file:',filepaths(xs)[f],' not found ! \n') 
  } 
  invisible(list(scantimes=scantimes,EIC=EIC)); 
}

getEIC4Peaks <- function(xraw,peaks,maxscans=length(xraw@scantime)){
  if (!is.double(xraw@env$mz) || !is.double(xraw@env$intensity) || !is.integer(xraw@scanindex)) stop('mz/int not double.')
  npeaks <- dim(peaks)[1]; 
  scans  <- length(xraw@scantime);
  eics <- matrix(NA,npeaks,maxscans);
  for (p in 1:npeaks) {
    timerange       <- c(peaks[p,"rtmin"],peaks[p,"rtmax"]);
    tidx <- which((xraw@scantime >= timerange[1]) & (xraw@scantime <= timerange[2]));
    if(length(tidx)>0){
      scanrange <- range(tidx);
    }else{
      scanrange <- 1:scans;
    }
    massrange <- c(peaks[p,"mzmin"],peaks[p,"mzmax"]);
    eic <- .Call("getEIC",xraw@env$mz,xraw@env$intensity,xraw@scanindex,as.double(massrange),
      as.integer(scanrange),as.integer(length(xraw@scantime)), PACKAGE ='xcms' )$intensity;
    eic[eic==0] <- NA;
    eics[p,scanrange[1]:scanrange[2]] <- eic; 
  }
eics
}

getAllEICs <- function(xs,index=NULL,file=NULL) {
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
  if(is.null(index)){
    cat("Missing Index, generate all EICs from sample 1.\n");
    index <- rep(1,nrow(peaki));
  }
  cat('Generating EIC\'s .. \n') 
  if (nfiles > 1) { 
#       cat('Searching maxima .. \n')
      for (f in 1:nfiles){
#         cat('Reading raw data file:',filepaths(xs)[f]) 
        xraw <- xcmsRaw(filepaths(xs)[f],profstep=0)
#         cat(',', length(xraw@scantime),'scans. \n') 
        maxscans <- max(maxscans,length(xraw@scantime))
        scantimes[[f]] <- xraw@scantime
      }
      EIC <- array(integer(0),c(nrow(peaki),maxscans,1))     
      for (f in 1:nfiles){
        if (file.exists(filepaths(xs)[f])) { 
#           cat('Reading raw data file:',filepaths(xs)[f],'\n') 
          xraw <- xcmsRaw(filepaths(xs)[f],profstep=0)
      #    cat('Generating EIC\'s .. \n') 
          idx.peaks <- which(index == f);
          if(length(idx.peaks)>0){
            pdata <- as.data.frame(xs@peaks[peaki[idx.peaks,f],]) # data for peaks from file f
            if(length(idx.peaks)==1){
              pdata <- t(pdata);
            }
  #           if (f==1) EIC <- array(integer(0),c(nrow(pdata),maxscans,length(filepaths(xs))))   
            EIC[idx.peaks,,1] <- getEICs(xraw,pdata,maxscans)
          }
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
