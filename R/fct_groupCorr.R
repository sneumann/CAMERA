##Functions for groupCorr

##Fast matrix generation
create.matrix <- function(dim1,dim2) {
  x <- matrix()
  length(x) <- dim1*dim2
  dim(x) <- c(dim1,dim2)
  x
}

##Calculate Correlation in Samples
##Needs: xsa object, EIC correlation matrix (from getAllPeakEICs), parameters
setGeneric("calcCiS", function(object, ...) standardGeneric("calcCiS"))
setMethod("calcCiS", "xsAnnotate", function(object, EIC=EIC, corval=0.75, 
                                              pval=0.05, psg_list=NULL ){

  #Columns Peak 1, Peak 2, correlation coefficienct, Pseudospectrum Index
  resMat <- create.matrix(100000,4);
  colnames(resMat) <- c("x","y","cor","ps")
  cnt <- 0;
  
  npeaks <- nrow(object@groupInfo);

  ncl <- sum(sapply(object@pspectra, length));
  npeaks.global <- 0; #Counter for % bar
  npspectra <- length(object@pspectra);

  #Check if object have been preprocessed with groupFWHM
  if(npspectra < 1){
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    cat('Calculating peak correlations in 1 group.\n % finished: '); 
    lp <- -1;
    pspectra_list     <- 1;
    object@psSamples  <- 1;
  }else{
    if(is.null(psg_list)){
      cat('\nCalculating peak correlations in',npspectra,'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating peak correlations in',length(psg_list),'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- psg_list;
      ncl <- sum(sapply(object@pspectra[psg_list],length));
    }
  }

  if(dim(EIC)[2] != npeaks){
      EIC <- t(EIC);
      #Second check, otherwise number of peaks != number of EIC curves
      if(dim(EIC)[2] != npeaks){
        stop(paste("Wrong dimension of EIC. It has ",dim(EIC)[1]," Rows for ",npeaks,"peaks",sep=""));
      }
  }
  lp <- -1;
  #Iterate over all PS-spectra
  for(j in 1:length(pspectra_list)){
    i  <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    npi <- length(pi);

    if( ((npi^2)/2 + cnt)  >= nrow(resMat)){
      #resize resMat
      resMat <- rbind(resMat,create.matrix(100000,4));
    }

    #percent output
    npeaks.global <- npeaks.global + length(pi);
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

    #Need at least two peaks for correlation
    if(npi < 2) {
      next;
    }

    #Suppress warnings
    options(warn = -1);
    res <- rcorr(EIC[, pi], type="pearson")
    options(warn = 0);

    #Set lower triangle to NA
    res$r[lower.tri(res$r,diag = TRUE)] <- NA;
    res$P[lower.tri(res$P,diag = TRUE)] <- NA;
    
    #Find peaks with have correlation higher corr_threshold and p <= 0.05
    index <- which(( res$r > corval) & (res$P <= pval))
    if(length(index) > 0){
     for( x in 1:(length(index))){
      col <- index[x] %/% npi + 1;
      row <- index[x] %%  npi;
      if(row == 0){
        row <- npi;
        col <- col - 1;
      }
      xi <- pi[row];
      yi <- pi[col];
      #y > x should be satisfied for every value, since we have set lower.tri to NA
      cnt<-cnt+1;
      resMat[cnt,] <- c(xi,yi,res$r[row, col],i);
      }
    }else{
      next;
    }
  }
  cat("\n");

  return(invisible(resMat[1:cnt,,drop=FALSE]))
})


setGeneric("calcCaS", function(object, ...) standardGeneric("calcCaS"))
setMethod("calcCaS", "xsAnnotate", function(object, corval=0.75, pval=0.05,
                                           intval="into") {
  #Calculate correlation across samples for a given xsAnnotate
  if (!sum(intval == c("into","intb","maxo"))){
       stop("unknown intensity value!")
  }
  npspectra <- length(object@pspectra);
  npeaks.global <- 0;
  ncl <- sum(sapply(object@pspectra, length));
  lp <- -1;
  #Columns Peak 1, Peak 2, correlation coefficienct, Pseudospectrum Index
  resMat <- create.matrix(100000,4);
  colnames(resMat) <- c("x","y","cor","ps")
  cnt <- 0;

  #Check that you have more than 1 sample
  #Precondition for correlation analysis
  if(length(object@xcmsSet@filepaths) < 3){
    cat('Calculation correlation across samples needs at least 3 samples\n');
    return(NULL);
  }else{
    #Check if object have been preprocessed with groupFWHM
    if(npspectra < 1){
      npspectra <- 1;
      object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
      cat('Calculating peak correlations across samples.\n % finished: '); 
      object@psSamples  <- 1;
    }else{
      cat('\nCalculating peak correlations across samples.\n % finished: '); 
    }
    npeaks <- 0;
    
    peaktable <- t(groupval(object@xcmsSet, value=intval));
    for(i in 1:npspectra){
      pi  <- object@pspectra[[i]];
      npi <- length(pi);
      if( ((npi^2)/2 + cnt)  >= nrow(resMat)){
        #resize resMat
        size <- max(100000, ((npi^2)/2 + 10000))
        resMat <- rbind(resMat,create.matrix(as.integer(size),4));
      }

      #percent output
      npeaks.global <- npeaks.global + length(pi);
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

      if(npi <2) {
        next;
      }

      #Suppress warnings
      options(warn = -1);
      res <- rcorr(peaktable[,pi], type="pearson");
      options(warn = 0);

      #Set lower triangle to NA
      res$r[lower.tri(res$r,diag = TRUE)] <- NA;
      res$P[lower.tri(res$P,diag = TRUE)] <- NA;
      
      #Find peaks with have correlation higher corr_threshold and p <= 0.05
      index <- which(( res$r > corval) & (res$P <= pval))
      if(length(index) > 0){
        for( x in 1:(length(index))){
          col <- index[x] %/% npi + 1;
          row <- index[x] %%  npi;
          if(row == 0){
            row <- npi;
            col <- col - 1;
          }
          xi <- pi[row];
          yi <- pi[col];
          #y > x should be satisfied for every value, since we have set lower.tri to NA
          cnt<-cnt+1;
          resMat[cnt,] <- c(xi,yi,res$r[row, col],i);
        }
      }else{
        next;
      }
    }
    cat("\n");
    return(invisible(resMat[1:cnt,,drop=FALSE]))
  }
})


setGeneric("calcIsotopes", function(object) standardGeneric("calcIsotopes"))
setMethod("calcIsotopes", "xsAnnotate", function(object){

  ncl <- sum(sapply(object@pspectra, length));
  npeaks.global <- 0; #Counter for % bar
  npspectra <- length(object@pspectra);
  #Columns Peak 1, Peak 2, correlation coefficienct, Pseudospectrum Index
  resMat <- matrix(nrow=0,ncol=4)
  colnames(resMat) <- c("x","y","cor","ps")

  if(length(object@isotopes) < 1){
    cat('Object contains no Isotope Information.\nRun findIsotopes before.\n');
    return(NULL);
  }else{
    if(nrow(object@isoID) > 0){
      cat('\nCalculating isotope assignments in',npspectra,'Groups... \n % finished: '); 
      lp <- -1;
      for(i in 1:npspectra){
          pi <- object@pspectra[[i]]
          npi <- length(pi)

          #percent output
          npeaks.global <- npeaks.global + length(pi);
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

          if(npi <2){
            next;
          }
          tmp <- matrix(nrow=0,ncol=2);
          invisible(lapply(pi,function(x) {
              if(!is.null(object@isotopes[[x]])) tmp <<- rbind(tmp,c(x,object@isotopes[[x]]$y))
            }))
          if(nrow(tmp)<1){
            next;
          }
          index.min <- min(tmp[,2]);
          index.max <- max(tmp[,2]);
          for(ii in index.min:index.max){
              tmp2 <- expand.grid(tmp[which(tmp[,2]==ii),1],tmp[which(tmp[,2]==ii),1],KEEP.OUT.ATTRS=FALSE)
              tmp2 <- tmp2[-which(tmp2[,1] >= tmp2[,2]),]
              resMat <- rbind(resMat,cbind(x=tmp2[,1],y=tmp2[,2],cor=1,ps=i));
          }
      }          
    }else{
      cat('Object contains no Isotopes!\n')
      return(resMat);
    }
  }
  return(resMat)  
})
##Methods for Combination of calc Results

setGeneric("combineCalc", function(object1,object2,...) standardGeneric("combineCalc"))

setMethod("combineCalc", signature("matrix","matrix"), function(object1,object2, method=getOption("BioC")$CAMERA$combineCalc.method,
                                       ...) {
    method <- match.arg(method, getOption("BioC")$CAMERA$combineCalc.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("combineCalc", method, sep=".")
    invisible(do.call(method, alist(object1,object2, ...)))
})


setGeneric("combineCalc.sum", function(object1,object2) standardGeneric("combineCalc.sum"))

setMethod("combineCalc.sum", signature("matrix","matrix"), function(object1,object2){
        
        if(ncol(object1) != 4){
          stop("first object is not a matrix with 4 columns");
        }
        if(ncol(object2) != 4){
          stop("second object is not a matrix with 4 columns");
        }
          
        combination = new.env(hash = TRUE)
        
        apply(object1,1,function(x){ 
          combination[[paste(x[c(1,2,4)],collapse=" ")]]<- x[3] 
        })

        apply(object2,1,function(x){
          if(is.null(combination[[paste(x[c(1,2,4)],collapse=" ")]])){
            combination[[paste(x[c(1,2,4)],collapse=" ")]]<- x[3]
          }else{
            combination[[paste(x[c(1,2,4)],collapse=" ")]]<- combination[[paste(x[c(1,2,4)],collapse=" ")]] + x[3];
          }
        })

        resMat <- matrix(ncol=4,nrow=length(ls(combination)));

        i<-1;y<-c();
        sapply(ls(combination), function(x) {
            y[c(1,2,4)] <- unlist(strsplit(x," "));
            y[3] <- combination[[x]];
            resMat[i,] <<- y; i<<-i+1;
        })
        resMat <- matrix(as.numeric(resMat),ncol=4);
        colnames(resMat) <- c("x","y","cor","ps")

        return(invisible(resMat));
})

##END Methods for Combination of calc Results


##Methods for Cluster Seperation of Pseudospectra

setGeneric("calcPC", function(object, method, ...) standardGeneric("calcPC"))

setMethod("calcPC", "xsAnnotate", function(object, method=getOption("BioC")$CAMERA$calcPC.method,
                                       ...) {

    method <- match.arg(method, getOption("BioC")$CAMERA$calcPC.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("calcPC", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

##label.propagation.community

setGeneric("calcPC.lpc", function(object, ...) standardGeneric("calcPC.lpc"))

setMethod("calcPC.lpc", "xsAnnotate", function(object, ajc=NULL,
                                               psg_list=NULL) {
  npspectra <- length(object@pspectra);
  pspectra  <- object@pspectra
  psSamples <- object@psSamples;
  npeaks.global <- 0
  ncl <- sum(sapply(object@pspectra, length));

  colnames(ajc)[3] <- c("weight") ##todo: Change to generell ajc interface

  #Information for % output
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
  if(is.null(ajc)){
    stop("Need ajc argument as weighted symbolic edge list. Example see manpage\n")
  }

  #peak counter
  npeaks <- 0;
  lp <- -1;
  for(j in 1:length(pspectra_list)){
    i  <- pspectra_list[j];#index of pseudospectrum
    pi <- object@pspectra[[i]]; #peak_id in pseudospectrum
    
          #percent output
          npeaks.global <- npeaks.global + length(pi);
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

    index <- which(ajc[,4] == i)
    if(length(index) < 1){
      g <- graph.data.frame(vertices=as.data.frame(pi),d=matrix(nrow=0,ncol=2), directed=FALSE)
    }else{
      g <- graph.data.frame(vertices=as.data.frame(pi),d=as.data.frame(ajc[index,1:3,drop=FALSE]), directed=FALSE)
    }
    lpc <- label.propagation.community(g,initial=0:(length(pi)-1))
    pspectra[[i]] <- pi[which(lpc==0)];
    if(max(lpc) > 0){
      for(ii in 1:max(lpc)){
        npspectra <- npspectra +1;
        pspectra[[npspectra]] <- pi[which(lpc==ii)];
        psSamples <- c(psSamples,psSamples[i])
      }
    }
  }
  object@pspectra  <- pspectra;
  object@psSamples <- psSamples;
  cat("\n");
  cat("New number of ps-groups: ",length(pspectra),"\n");
  return(object)
})

##highly connected subgraph
setGeneric("calcPC.hcs", function(object, ...) standardGeneric("calcPC.hcs"))

setMethod("calcPC.hcs", "xsAnnotate", function(object, ajc=NULL,
                                               psg_list=NULL) {

  npspectra <- length(object@pspectra);
  pspectra  <- object@pspectra
  psSamples <- object@psSamples;
  npeaks.global <- 0;
  ncl <- sum(sapply(object@pspectra, length));
  colnames(ajc)[3] <- c("weight") ##todo: Change to generell ajc interface

  #Information for % output
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

  #peak counter
  npeaks <- 0;
  lp <- -1;
  for(j in 1:length(pspectra_list)){
    i  <- pspectra_list[j];#index of pseudospectrum
    pi <- object@pspectra[[i]]; #peak_id in pseudospectrum
    
#percent output
          npeaks.global <- npeaks.global + length(pi);
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

    index <- which(ajc[,4] == i)
    if(length(index) < 1){
      g <- ftM2graphNEL(matrix(nrow=0,ncol=2),V=as.character(pi),edgemode="undirected")
    }else{
      g <- ftM2graphNEL(ajc[index,1:2,drop=FALSE],W=ajc[index,3,drop=FALSE], V=as.character(pi), edgemode="undirected");
    }
    hcs <- highlyConnSG(g);
    
    #order cluster after size
    cnts <- sapply(hcs$clusters,length);
    grps <- 1:length(hcs$clusters);     
    grps <- grps[order(cnts, decreasing = TRUE)]

    for (ii in 1:length(grps)){
      if(ii==1){
        #save old pspectra number
        pspectra[[i]] <- as.integer(hcs$cluster[[grps[ii]]])
      } else {
        npspectra <- npspectra +1
        pspectra[[npspectra]] <- as.integer(hcs$cluster[[grps[ii]]])
        psSamples[npspectra] <- psSamples[i]
      }
    }
  }
  object@pspectra  <- pspectra;
  object@psSamples <- psSamples;
  cat("\n");
  cat("New number of ps-groups: ",length(pspectra),"\n");
  return(object)
})

##END Methods for Cluster Seperation of Pseudospectra
calcCL3 <- function(object, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th, psg_list=psg_list){

  nrow <- nrow(object@groupInfo);
  CL   <- vector("list", nrow);
  CIL  <- list();
  CI   <- NULL;
  ncl  <- length(CL);
  
  npeaks <- 0; #Counter for % bar
  npspectra <- length(object@pspectra);
  
  #Check if object have been preprocessed with groupFWHM
  if(npspectra < 1){
    npspectra <- 1;
    object@pspectra[[1]] <- seq(1:nrow(object@groupInfo));
    cat('Calculating peak correlations for 1 big group.\nTry groupFWHM before, to reduce runtime. \n% finished: '); 
    lp <- -1;
    pspectra_list     <- 1;
    object@psSamples  <- 1;
  }else{
    if(is.null(psg_list)){
      cat('\nCalculating peak correlations in',npspectra,'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- 1:npspectra;
    }else{
      cat('\nCalculating peak correlations in',length(psg_list),'Groups... \n % finished: '); 
      lp <- -1;
      pspectra_list <- psg_list;
      ncl <- sum(sapply(object@pspectra[psg_list],length));
    }
  }

  psSamples <- object@psSamples;

  for(j in 1:length(pspectra_list)){
    i  <- pspectra_list[j];
    pi <- object@pspectra[[i]];
    f  <- psSamples[j]; #?????
    npi <- length(pi);

    #percent output
    npeaks <- npeaks + length(pi);
    perc   <- round((npeaks) / ncl * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { 
      cat(perc,' '); 
      lp <- perc;
    }
    if (.Platform$OS.type == "windows"){ 
      flush.console();
    }
    #end percent output

    if(npi < 2) {
      next;
    }

    #Suppress warnings
    options(warn = -1);
    res <- rcorr(EIC[, pi], type="pearson")
    options(warn = 0);

    #Set lower triangle to NA
    res$r[lower.tri(res$r,diag = TRUE)] <- NA;
    res$P[lower.tri(res$P,diag = TRUE)] <- NA;
    
    #Find peaks with have correlation higher corr_threshold and p <= 0.05
    index <- which(( res$r > cor_eic_th) & (res$P <= 0.05))
    if(length(index) > 0){
     for( x in 1:(length(index))){
      col <- index[x] %/% npi + 1;
      row <- index[x] %%  npi;
      if(row == 0){
        row <- npi;
        col <- col - 1;
      }
      xi <- pi[row];
      yi <- pi[col];
      CL[[xi]] <- c(CL[[xi]], yi)
      CL[[yi]] <- c(CL[[yi]], xi) ## keep the list symmetric
      CIL[[length(CIL)+1]] <- list(p=c(xi,yi), cor=res$r[row, col])
      }
    }else{
      next;
    }
  }

  if (length(CIL) >0){ 
    CI <- data.frame(t(sapply(CIL, function(x) x$p)), sapply(CIL, function(x) x$cor) );
  } else { 
    return(NULL);
  }
  colnames(CI) <- c('xi', 'yi', 'cors');
  cat("\n");
  return(invisible(list(CL=CL, CI=CI)))
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

getMaxScans <- function(object){

 if (!class(object) == "xsAnnotate"){
    stop ("no xsAnnotate object");
  }

  nfiles <- length(filepaths(object@xcmsSet))
  maxscans <- 0
  if(nfiles == 1){
    if (file.exists(filepaths(object@xcmsSet)[1])) { 
      xraw <- xcmsRaw(filepaths(object@xcmsSet)[1],profstep=0)
      maxscans <- length(xraw@scantime)     
    }else {
      stop('Raw data file:',filepaths(xs)[1],' not found ! \n');
    }
  }else {
    #Get scantime length for every xraw
    for (f in 1:nfiles){
      if(file.exists(filepaths(object@xcmsSet)[f])) { 
      xraw <- xcmsRaw(filepaths(object@xcmsSet)[f], profstep=0);
      maxscans <- max(maxscans, length(xraw@scantime));
      } else {
        stop('Raw data file:',filepaths(xs)[f],' not found ! \n');
      }
    }
  }
  return(maxscans)
}

setGeneric("getAllPeakEICs", function(object, index, maxscans) standardGeneric("getAllPeakEICs"))

setMethod("getAllPeakEICs", "xsAnnotate", function(object, index=NULL, maxscans=NULL){

  if(is.null(maxscans)){
    stop("maxscans is not set. Use getMaxScans beforehand.\n")
  }
  nfiles <- length(filepaths(object@xcmsSet))
  scantimes <- list()
#   maxscans <- maxscans;
  if(nfiles == 1){
    #Single File
    if (file.exists(filepaths(object@xcmsSet)[1])) { 
     xraw <- xcmsRaw(filepaths(object@xcmsSet)[1],profstep=0)
     maxscans <- length(xraw@scantime)
     scantimes[[1]] <- xraw@scantime
     pdata <- as.data.frame(object@xcmsSet@peaks) 
     EIC <- create.matrix(nrow(pdata),maxscans)
     EIC[,] <- getEIC4Peaks(xraw,pdata,maxscans)
    } else {
      stop('Raw data file:',filepaths(xs)[f],' not found ! \n');
    }
  } else {
    gval <- groupval(object@xcmsSet);
    if(is.null(index)){
      cat("Missing Index, generate all EICs from sample 1.\n");
      index <- rep(1, nrow(gval));
    }else if(length(index) != nrow(gval)){
      stop("Index length must equals number of peaks.\n");
    }

    cat('Generating EIC\'s .. \n')
  
    #Get scantime length for every xraw
#     for (f in 1:nfiles){
#       xraw <- xcmsRaw(filepaths(object@xcmsSet)[f], profstep=0);
#       maxscans <- max(maxscans, length(xraw@scantime));
#       scantimes[[f]] <- xraw@scantime;
#     }

    #generate EIC Matrix
    EIC <- create.matrix(nrow(gval),maxscans)

    #na flag, stores if sample contains NA peaks
    na.flag <- 0;

    for (f in 1:nfiles){
      if (file.exists(filepaths(object@xcmsSet)[f])) { 
        xraw <- xcmsRaw(filepaths(object@xcmsSet)[f], profstep=0);
        idx.peaks <- which(index == f);
        if(length(idx.peaks) > 0){
          pdata <- as.data.frame(object@xcmsSet@peaks[gval[idx.peaks,f],]) # data for peaks from file f
          if(length(which(is.na(pdata[,1]))) >0){
            na.flag <- 1;
          }
          if(length(idx.peaks)==1){
            pdata <- t(pdata);
          }
          EIC[idx.peaks,] <- getEIC4Peaks(xraw,pdata,maxscans)
        }
      } else {
        stop('Raw data file:',filepaths(xs)[f],' not found ! \n')
      }
    }
    if(na.flag ==1){
      cat("Warning: Found NA peaks in selected sample.\n");
    }
  }
  invisible(list(scantimes=scantimes,EIC=EIC)); 
})

# 
# getAllPeakEICs <- function(xs,index=NULL){
# 
#   peaki <- getPeaksIdxCol(xs,NULL)
#   nfiles <- length(filepaths(xs))
#   scantimes <- list()
#   maxscans <- 0
# 
#   if(is.null(index)){
#     cat("Missing Index, generate all EICs from sample 1.\n");
#     index <- rep(1,nrow(peaki));
#   }
#   cat('Generating EIC\'s .. \n') 
#   if (nfiles > 1) { 
#     for (f in 1:nfiles){
#       xraw <- xcmsRaw(filepaths(xs)[f], profstep=0);
#       maxscans <- max(maxscans, length(xraw@scantime));
#       scantimes[[f]] <- xraw@scantime;
#     }
#     EIC <- array(integer(0), c(nrow(peaki),maxscans))
#     na.flag <- 0;
#     for (f in 1:nfiles){
#       if (file.exists(filepaths(xs)[f])) { 
#         xraw <- xcmsRaw(filepaths(xs)[f], profstep=0);
#         idx.peaks <- which(index == f);
#         if(length(idx.peaks) > 0){
#             pdata <- as.data.frame(xs@peaks[peaki[idx.peaks,f],]) # data for peaks from file f
#             if(length(which(is.na(pdata[,1]))) >0){
#               na.flag <- 1;
#             }
#             if(length(idx.peaks)==1){
#               pdata <- t(pdata);
#             }
#             EIC[idx.peaks,] <- getEIC4Peaks(xraw,pdata,maxscans)
#           }
#         }
#         else stop('Raw data file:',filepaths(xs)[f],' not found ! \n')
#       }
#     if(na.flag ==1){
#       cat("Found NA peaks in selected samples. Those will be seperated in new pcgroups in each case.\nUse fillpeaks if not desired!\n");
#     }
#   }  else { ## create EIC's for single file
#        if (file.exists(filepaths(xs)[1])) { 
#           xraw <- xcmsRaw(filepaths(xs)[1],profstep=0)
#           maxscans <- length(xraw@scantime)
#           scantimes[[1]] <- xraw@scantime
#           pdata <- as.data.frame(xs@peaks[peaki,]) 
#           EIC <- array(NA,c(nrow(pdata),maxscans))   
#           EIC[,] <- getEIC4Peaks(xraw,pdata,maxscans)
#         }  else stop('Raw data file:',filepaths(xs)[f],' not found ! \n') 
#   } 
#   invisible(list(scantimes=scantimes,EIC=EIC)); 
# }
# 

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


fast_corr <- function(x){
  
  x[is.na(x)] <- 1e30; #same factor as in rcorr

  p <- as.integer(ncol(x))
  if(p<1)
    stop("must have >1 column")
  
  n <- as.integer(nrow(x))
  if(n<5)
    stop("must have >4 observations")
  
  x <- scale(x);
  r <- crossprod(x) / (n-1);
  
  r[r>1e29] <- NA; #sace factor as in rcorr

  npair <- matrix(rep(n,p*p),ncol=p)

  P <- matrix(2*(1-pt(abs(r)*sqrt(npair-2)/sqrt(1-r*r), npair-2)),ncol=p);
  
  P[abs(r)==1] <- 0;
  diag(P) <- NA;
  invisible(list(r=r, n=npair, P=P))
}
