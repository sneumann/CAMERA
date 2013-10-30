##############################################################
## private class to calculate (dynamic) rulesets
##

setClass("ruleSet",
         representation(ionlistfile="character",
                        ionlist="data.frame", 
                        neutrallossfile="character",
                        neutralloss="data.frame", 
                        neutraladditionfile="character",
                        neutraladdition="data.frame",
                        maxcharge="numeric",
                        mol="numeric",
                        nion="numeric",
                        nnloss="numeric",
                        nnadd="numeric",
                        nh="numeric",
                        polarity="character",
                        rules="data.frame",
                        lib.loc="character"),        
         contains=c("Versioned"),
         prototype=prototype(
           ionlistfile="",
           ionlist=data.frame(),
           neutrallossfile="",
           neutralloss=data.frame(),           
           neutraladditionfile="",
           neutraladdition=data.frame(),
           maxcharge=numeric(),
           mol=numeric(),
           nion=numeric(),
           nnloss=numeric(),
           nnadd=numeric(),
           nh=numeric(),
           polarity=NULL,
           rules=data.frame(),
           lib.loc=NULL,
           new("Versioned", versions=c(ruleSet="0.1.1"))),
         validity=function(object) {
           TRUE
         })

setGeneric("setDefaultLists", function(object, lib.loc=.libPaths()) standardGeneric("setDefaultLists"))
setGeneric("readLists", function(object) standardGeneric("readLists"))
setGeneric("setDefaultParams", function(object) standardGeneric("setDefaultParams"))
setGeneric("setParams", function(object, maxcharge, mol, nion, nnloss, nnadd, nh, 
                                 polarity, lib.loc) standardGeneric("setParams"))
setGeneric("generateRules", function(object) standardGeneric("generateRules"))
setGeneric("generateRules2", function(object) standardGeneric("generateRules2"))

setMethod("show",
          signature="ruleSet",
          function(object) {
            cat("Backing files: ", object@ionlistfile, object@neutrallossfile, object@neutraladditionfile,"\n")
            cat("Ion lists: ",
                "  ", nrow(object@ionlist), " ions",
                "  ", nrow(object@neutralloss), " neutral losses",
                "  ", nrow(object@neutraladdition), " neutral additions",
                "\n")
            object
          })


setMethod("setDefaultLists",
          signature="ruleSet",
          function(object, lib.loc=.libPaths()) {
            
            ##Read Tabellen
            object@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA",
                                              lib.loc=lib.loc)[1]
            if (!file.exists(object@ionlistfile)) stop('ions.csv not found.')
            
            object@neutrallossfile <- system.file('lists/neutralloss.csv', 
                                                  package = "CAMERA", lib.loc=lib.loc)[1]
            if (!file.exists(object@neutrallossfile)) stop('neutralloss.csv not found.')
            
            object@neutraladditionfile <- system.file('lists/neutraladdition.csv', 
                                                      package = "CAMERA", lib.loc=lib.loc)[1]
            if (!file.exists(object@neutraladditionfile)) stop('neutraladdition.csv not found.')
            object
          })

setMethod("readLists",
          signature="ruleSet",
          function(object) {
            
            object@ionlist <- read.table(object@ionlistfile, header=TRUE, dec=".", sep=",",
                                         as.is=TRUE, stringsAsFactors = FALSE);
            
            object@neutralloss <- read.table(object@neutrallossfile, header=TRUE, dec=".", sep=",",
                                             as.is=TRUE, stringsAsFactors = FALSE);
            
            object@neutraladdition <- read.table(object@neutraladditionfile,
                                                 header=TRUE, dec=".", sep=",",
                                                 as.is=TRUE, stringsAsFactors = FALSE);            
            object
          })

setMethod("setDefaultParams",
          signature="ruleSet",
          function(object) {
            object@maxcharge=3
            object@mol=3
            object@nion=2
            object@nnloss=1
            object@nnadd=1
            object@nh=2
            object@polarity="positive"
            object
          })

setMethod("setParams",
          signature=c("ruleSet", "numeric","numeric","numeric","numeric","numeric",
                      "numeric","character","character"), 
          function (object, maxcharge=3, mol=3, nion=2, nnloss=1, nnadd=1, nh=2, 
                    polarity=NULL, lib.loc=NULL) {
            object@maxcharge=maxcharge 
            object@mol=mol
            object@nion=nion
            object@nnloss=nnloss 
            object@nnadd=nnadd 
            object@nh=nh
            object@polarity=polarity
            object@lib.loc=lib.loc
            object
          })

setMethod("generateRules",
          signature="ruleSet", 
          function (object) {

            maxcharge=object@maxcharge
            mol=object@mol
            nion=object@nion
            nnloss=object@nnloss
            nnadd=object@nnadd 
            nh=object@nh
            polarity=object@polarity

            ionlist=object@ionlist
            neutralloss=object@neutralloss
            neutraladdition=object@neutraladdition              

            rules=object@rules
            
            name<-c();
            nmol<-c();
            charge<-c();
            massdiff<-c();
            oidscore<-c();
            quasi<-c();
            ips<-c();

            ##Erzeuge Regeln
            tmpname   <- c();
            tmpnmol   <- c();
            tmpcharge <- 0;
            tmpmass   <- 0;
            tmpips    <- 0;

            ## Molekülionen
            if(polarity=="positive"){
              ## Wasserstoff, hard codiert
              for(k in 1:mol){
                if(k == 1){
                  str    <- "";
                  tmpips <- 1.0;
                }else{
                  str    <-  k;
                  tmpips <- 0.5;
                };
                name     <- append(name, paste("[", str, "M+H]+", sep=""));
                charge   <- append(charge, 1);
                massdiff <- append(massdiff, 1.007276);
                nmol     <- append(nmol, k);
                if(k == 1) {
                  quasi  <- append(quasi, 1);
                } else { 
                  quasi  <- append(quasi, 0);
                };
                oidscore <- append(oidscore, 1);
                ips      <- append(ips, tmpips)
                name     <- append(name,paste("[",str,"M+2H]2+",sep=""));
                charge<-append(charge,2);
                massdiff<-append(massdiff,2.014552);
                nmol<-append(nmol,k);quasi<-append(quasi,0);
                oidscore<-append(oidscore,2);
                ips<-append(ips,tmpips)
                name<-append(name,paste("[",str,"M+3H]3+",sep=""));
                charge<-append(charge,3);
                massdiff<-append(massdiff,3.021828);
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
                  massdiff <- append(massdiff,ionlist[i,3]+1.007276);
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
                  ## Austausch erstmal nur einen pro Ion
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
                    ## Austausch hat stattgefunden, einfach bsp 1
                    tmpname<-paste(tmpname,"-H",sep="");
                    tmpcharge<-tmpcharge-1;
                    tmpmass<-tmpmass-1.007276;
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

              ## Erzeuge Neutral Addition
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
                massdiff<- append(massdiff,neutraladdition[i,2]+1.007276);
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
                  massdiff<-  append(massdiff,-neutralloss[i,2]+1.007276*ii);
                  nmol<-append(nmol,1);
                  quasi<-append(quasi,0)
                  oidscore<-append(oidscore,oid+1);oid<-oid+1;
                  ips<-append(ips,0.25);
                }
              }
              ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
              if(length(index<-which(ruleset[,"charge"]>maxcharge))>0){
                ruleset<- ruleset[-index,];
              }
            }else if(polarity=="negative"){
              ## Wasserstoff, hard codiert
              for(k in 1:mol){
                if(k==1){str<-"";tmpips<-1.0;}else{str<-k;tmpips<-0.5};
                name<-append(name,paste("[",str,"M-H]-",sep=""));
                charge<-append(charge,-1);massdiff<-append(massdiff,-1.007276);nmol<-append(nmol,k);if(k==1){quasi<-append(quasi,1);}else{quasi<-append(quasi,0);};oidscore<-append(oidscore,1);ips<-append(ips,tmpips)
                name<-append(name,paste("[",str,"M-2H]2-",sep=""));charge<-append(charge,-2);massdiff<-append(massdiff,-2.014552);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,2);ips<-append(ips,tmpips)
                name<-append(name,paste("[",str,"M-3H]3-",sep=""));charge<-append(charge,-3);massdiff<-append(massdiff,-3.021828);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,3);ips<-append(ips,tmpips)
                oid<-3;
                for(i in 1:nrow(ionlist)){
                  if(ionlist[i,2]>=0){
                    if(ionlist[i,2] == 1){
                      name<-append(name,paste("[",str,"M-2H+",ionlist[i,1],"]-",sep=""));
                      charge <- append(charge,ionlist[i,2]-2);
                      massdiff<- append(massdiff,ionlist[i,3]-(2*1.007276));
                      nmol <- append(nmol,k);
                      quasi <- append(quasi,0);
                      oidscore<-append(oidscore,oid+i);
                      ips<-append(ips,0.25);
                      next;
                    } else {
                      if(ionlist[i,2] > maxcharge) {next;}
                      localCharge <- ionlist[i,2]+1
                      name<-append(name,paste("[",str,"M-",localCharge,"H+",ionlist[i,1],"]-",sep=""));
                      charge <- append(charge,ionlist[i,2]-localCharge);
                      massdiff<- append(massdiff,ionlist[i,3]-(localCharge*1.007276));
                      nmol <- append(nmol,k);
                      quasi <- append(quasi,0);
                      oidscore<-append(oidscore,oid+i);
                      ips<-append(ips,0.25);
                      next;
                    }
                  }
                  if(ionlist[i,2]== -1){
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]2-",sep=""));
                  }else{
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]",ionlist[i,2]+1,"-",sep=""));
                  }
                  charge <- append(charge,ionlist[i,2]-1);
                  massdiff<- append(massdiff,ionlist[i,3]-1.007276);
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
                    tmpmass<-tmpmass-1.007276;
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
                massdiff<- append(massdiff,neutraladdition[i,2]-1.007276);
                nmol<-append(nmol,1);
                quasi<-append(quasi,0)
                oidscore<-append(oidscore,oid+1);oid<-oid+1;
                ips<-append(ips,0.5);
              }
              oid<-max(oidscore);

              ##Erzeuge Neutral loss
              index<-which(quasi==1)
              for(i in 1:nrow(neutralloss)){
                name<-append(name,paste("[M-H-",neutralloss[i,1],"]-",sep=""));
                charge<-append(charge,-1);
                massdiff<-  append(massdiff,-neutralloss[i,2]-1.007276);
                nmol<-append(nmol,1);
                quasi<-append(quasi,0)
                oidscore<-append(oidscore,oid+1);oid<-oid+1;
                ips<-append(ips,0.25);
              }
              ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
              if(length(index<-which(ruleset[,"charge"]< -maxcharge))>0){
                ruleset<- ruleset[-index,];
              }
            }else stop("Unknown error")

            ## Update object rules and return ruleset
            object@rules=ruleset
            object;
          })


setMethod("generateRules2", signature="ruleSet", function (object) {

  maxcharge=object@maxcharge
  mol=object@mol
  nion=object@nion
  nnloss=object@nnloss
  nnadd=object@nnadd 
  nh=object@nh
  polarity=object@polarity

  ionlist=object@ionlist
  neutralloss=object@neutralloss
  neutraladdition=object@neutraladdition              

  rules=object@rules
            
  ##Create Rules
  ruleset     <- matrix(nrow=0, ncol=8);
  colnames(ruleset) <- c("name", "nmol","charge", "massdiff", "typ","mandatory",
                         "score", "parent")

  tmpname   <- c();
  tmpcharge <- 0;
  tmpmass   <- 0;
  tmpionparent <- NA;
  massH <- 1.007276
  
  ##Positive Rule set
  if(polarity == "positive"){
    charge <- "+"
    chargeValue <- 1
  }else if(polarity == "negative"){
    charge <- "-"
    chargeValue <- -1
  }else{
    stop("Unknown polarity mode in rule set generation! Debug!\n")
  }

  #Hydrogen part is hard coded
  for(k in 1:mol) {
    if(k == 1){
      #For M+H
      str    <- "";
      tmpips <- 1.5;
      quasi  <- 1 
    }else{
      #For xM+H
      str    <-  k;
      tmpips <- 0;
      quasi  <- 0 
    }

    for(xh in seq.int(nh)){
      if(xh == 1){
        ruleset <- rbind(ruleset, cbind(paste("[", str, "M", charge, "H]", charge, sep=""), k, chargeValue, 
                                        massH*chargeValue, "A", 
                                        quasi, tmpips+0.5, tmpionparent));
       } else {
        ruleset <- rbind(ruleset, cbind(paste("[", str, "M", charge, xh, "H]", xh, charge, sep=""), 
                                        k,xh*chargeValue, massH*xh*chargeValue, "A", 0, 0.5, tmpionparent+1));
      }
    }
        
    for(i in 1:nrow(ionlist)){
      #Add change of kat- respectively anion against one hydrogen
      if(ionlist[i, 2] > 0 & charge == "-"){
        if(ionlist[i, 2] == 1){
          sumcharge <- "";
          hdiff <- 1;
        }else{
          sumcharge <- ionlist[i, 2];
          hdiff <- ionlist[i,2]-1;
        }
        ruleset <- rbind(ruleset, cbind(paste("[", str, "M", charge,"H+", ionlist[i,1], "-", sumcharge,"H]", "",
                                               charge, sep=""), k, ionlist[i, 2]*chargeValue, 
                                        ionlist[i, 3]-massH*(1+hdiff),
                                        "A", 0, 0.25, tmpionparent));  
      }
        
      #xM+H + additional Kation like Na or K      
      if(ionlist[i,2] <= 0 & chargeValue > 0) {
        #Ions with negative charge, when polarity is positive
        next;
      } else if(ionlist[i, 2] >= 0 & chargeValue < 0){
        #Ions with positive charge, when polarity is negative
        next;
      }

      #Add Rule to Ruleset
      if(abs(ionlist[i, 2]) == 1){
        ruleset <- rbind(ruleset, cbind(paste("[", str, "M", charge,"H+", ionlist[i,1], "]2", 
                                              charge, sep=""), k, ionlist[i, 2]+1*chargeValue, 
                                        ionlist[i, 3]+massH*chargeValue, "A", 0, 0.25, tmpionparent));
      }else{
        ruleset <- rbind(ruleset, cbind(paste("[", str, "M", charge,"H+", ionlist[i, 1], "]", 
                                              ionlist[i, 2]+(1*chargeValue),
                                        charge, sep=""), k, ionlist[i, 2]+1*chargeValue, 
                                        ionlist[i, 3]+massH*chargeValue, "A" ,0, 0.25, tmpionparent));
      }
          
    }#End for loop nrow(ionlist)
      
    ##Coeff - coefficient Matrix, for generating rules with
    #combination of kat- or anionsexchange ions like [M-2H+Na] (M-H and change of H against Na)
    coeff <- expand.grid(rep(list(0:nion), nrow(ionlist)))
    if(chargeValue > 0){
      index <- which(ionlist[, 2] <= 0);
    }else{
      index <- which(ionlist[, 2] >= 0);
    }

    if(length(index) > 0){
      coeff[, index] <- 0;
    }
    
    tmp <- NULL;
    for(i in 1:nrow(ionlist)){
      if((chargeValue > 0 & ionlist[i, 2] <= 0) | (chargeValue < 0 & ionlist[i,2] >=0)){
        next;
      }
      #Austausch erstmal nur einen pro Ion
      tmp <- rbind(tmp, t(apply(coeff, 1, function(x) {
                                            x[i] <- x[i]+1;
                                            x[nrow(ionlist)+1] <- -1;
                                            x }
                   )));
    }
    coeff <- cbind(coeff, rep(0, nrow(coeff)));
    colnames(coeff)[4] <- "Var4"
    colnames(tmp)[4] <- "Var4"
    coeff <- unique(rbind(coeff, tmp));
    
    for(i in 1:nrow(coeff)){
      if(sum(coeff[i, 1:nrow(ionlist)]) > 2 | sum(coeff[i, 1:nrow(ionlist)]) < 1){
        next;
      }
          
      tmpname   <- paste("[",str,"M",sep="");
      tmpcharge <- 0;
      tmpmass   <- 0;

      for(ii in 1:(ncol(coeff))){
        if(coeff[i,ii] > 0){
          if(coeff[i,ii] > 1){
            tmpname <- paste(tmpname, "+", coeff[i, ii], ionlist[ii, 1], sep="");
          }else{
            tmpname <- paste(tmpname, "+", ionlist[ii, 1], sep="");
          }
          tmpcharge <- tmpcharge + coeff[i, ii] * ionlist[ii, 2];
          tmpmass   <- tmpmass + coeff[i, ii] * ionlist[ii, 3];
        } else if (coeff[i,ii] < 0){
          tmpname <- paste(tmpname, "-H", sep="");
          tmpcharge <- tmpcharge - 1;
          tmpmass   <- tmpmass  - massH;
        }
      }

      if(abs(tmpcharge) > 1){
        tmpname <- paste(tmpname, "]", tmpcharge, charge, sep="")
      }else{
        tmpname <- paste(tmpname, "]", charge, sep="")
      }

      if(tmpcharge > maxcharge | tmpcharge == 0){
        next;
      }

      if(sum(coeff[i, ]) == 1 && k == 1 && coeff[i, 4] >=0){
        ruleset <- rbind(ruleset, cbind(tmpname, k, tmpcharge, tmpmass, "A", 1, 0.75, tmpionparent));
      }else{
        ruleset <- rbind(ruleset, cbind(tmpname, k, tmpcharge, tmpmass, "A", 0, 0.25, tmpionparent));
      }
    }#end for loop nrow(coeff)
  }#end for loop k

  # Create neutral addition to M+H from list
  for(i in 1:nrow(neutraladdition)){
    #Add neutral ion to only M+H
    ruleset <- rbind(ruleset, cbind(paste("[M", charge, "H+", neutraladdition[i, 1], "]", charge, sep="") , 1, chargeValue, 
                                    neutraladdition[i, 2]+(massH*chargeValue), "A", 0, 0.25, 1));
  }

  ## Add neutral loss from list to ruleset
  for(i in 1:nrow(neutralloss)){
    ruleset <- rbind(ruleset, cbind(paste("-", neutralloss[i, 1], sep=""), 1, 0, 
                                    -neutralloss[i, 2], "F", 0, 0.25, 1));
    #Eliminate rules with charge > maxcharge
    if(length(index <- which(ruleset[, "charge"] > maxcharge)) > 0){
      ruleset <- ruleset[-index, ];
    }
  }
  ruleset <- as.data.frame(ruleset, stringsAsFactors=FALSE)
  class(ruleset$nmol) <- "numeric"
  class(ruleset$charge) <- "numeric"
  class(ruleset$massdiff) <- "numeric"
  class(ruleset$mandatory) <- "numeric"
  class(ruleset$score) <- "numeric"
  class(ruleset$parent) <- "numeric"
  object@rules=ruleset
  object;

  return(object);
})


###############################################
##
## calcRules() convenience method with the old
## behaviour and signature
##

calcRules <- function (maxcharge=3, mol=3, nion=2, nnloss=1,
                       nnadd=1, nh=2, polarity=NULL, lib.loc = .libPaths(), newFragments=FALSE){

  r <- new("ruleSet")
  r <- setDefaultLists(r, lib.loc=lib.loc)
  r <- readLists(r)  
  r <- setParams(r, maxcharge=maxcharge, mol=mol, nion=nion,
                  nnloss=nnloss, nnadd=nnadd, nh=nh, polarity=polarity, lib.loc = lib.loc)  
  if(newFragments){
    r <- generateRules2(r)
  }else{
    r <- generateRules(r)
  }
  
  return(r@rules)
}
