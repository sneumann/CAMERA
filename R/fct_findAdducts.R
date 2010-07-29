##Functions for findAdducts


annotateGrpMPI <- function(params){
library(CAMERA);
result<-list();
  for(ii in 1:length(params$i)){
    result[[ii]]<-CAMERA:::annotateGrp(params$pspectra,params$i[[ii]],params$imz,params$rules,params$mzabs,params$devppm,params$isotopes,params$quasimolion);
  }
return(result);
}

annotateGrp <- function(ipeak,imz,rules,mzabs,devppm,isotopes,quasimolion) {
  mz     <- imz[ipeak];
  na_ini <- which(!is.na(mz))

  ML     <- massDiffMatrix(mz[na_ini],rules)
  hypothese <- createHypothese(ML,rules,devppm,mzabs,na_ini);

#   m      <- fastMatch(as.vector(ML), as.vector(ML), tol = max(2*devppm*mean(mz, na.rm=TRUE))+ mzabs)
#   c      <- sapply(m, length)
#   index  <- which(c >= 2)
#   if(length(index) == 0) {
#     return(NULL);
#   }
  
  #Erstelle Hypothesen
#   hypothese <- create_hypothese(m, index, ML, rules, na_ini)
  if(is.null(nrow(hypothese))){
    return(NULL);
  }
  
  #Entferne Hypothesen, welche gegen Isotopenladungen verstossen!
  if(length(isotopes) > 0){
      hypothese <- check_isotopes(hypothese, isotopes, ipeak)
  }
  if(nrow(hypothese) < 2){
    return(NULL);
  };
  
  #Test auf Quasi-Molekülionen
  hypothese <- check_quasimolion(hypothese, quasimolion)
  if(nrow(hypothese) < 2){
    return(NULL);
  };
  
  #Entferne Hypothesen, welche gegen OID-Score&Kausalität verstossen!
  hypothese <- check_oid_causality(hypothese, rules)
  if(nrow(hypothese) < 2){
    return(NULL);
  };
  
  #Prüfe IPS-Score
  hypothese <- check_ips(hypothese)
  if(nrow(hypothese)<2){return(NULL);};
  
  return(hypothese);
}


createHypothese <- function(ML,rules,devppm,mzabs,na_ini){
  ML.nrow <- nrow(ML);
  ML.vec <- as.vector(ML);
  max.value <- max(round(ML,0));
  hashmap <- vector(mode="list",length=max.value);
  for(i in 1:length(ML)){
    val <- trunc(ML[i],0);
    if(val>1){
      hashmap[[val]] <- c(hashmap[[val]],i);
    }
  }

  hypothese <- matrix(NA,ncol=9,nrow=0);
  colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "oidscore", "ips", "massgrp", "check");
  massgrp <- 1;

  for(i in 1:length(hashmap)){
    if(is.null(hashmap[[i]])){
      next;
    }
    candidates <- ML.vec[hashmap[[i]]];
    candidates.index <- hashmap[[i]];
    if(i != 1 && !is.null(hashmap[[i-1]]) && min(candidates) < i+(2*devppm*i+mzabs)){
      index <- which(ML.vec[hashmap[[i-1]]]> i-(2*devppm*i+mzabs))
      if(length(index)>0) {
        candidates <- c(candidates, ML.vec[hashmap[[i-1]]][index]);
        candidates.index <- c(candidates.index,hashmap[[i-1]][index]);
      }    
    }
    if(length(candidates) < 2){
      next;
    }
    tol <- max(2*devppm*mean(candidates, na.rm=TRUE))+ mzabs;
    result <- cutree(hclust(dist(candidates)),h=tol);
    index <- which(table(result) >= 2);
    if(length(index) == 0){
      next;
    }
    m <- lapply(index, function(x) which(result == x));
    for(ii in 1:length(m)){
        ini.adducts <- candidates.index[m[[ii]]];
        for( iii in 1:length(ini.adducts)){
          adduct <- ini.adducts[iii] %/% ML.nrow +1;
          mass   <- ini.adducts[iii] %% ML.nrow;
          if(mass == 0){
            mass <- ML.nrow;
            adduct <- adduct -1;
          }
          hypothese <- rbind(hypothese, c(na_ini[mass], adduct, rules[adduct, "nmol"], rules[adduct, "charge"], mean(candidates[m[[ii]]]), rules[adduct, "oidscore"], rules[adduct,"ips"],massgrp ,1));
        }
        massgrp <- massgrp +1;
    }
  }
  return(hypothese);
}



create_hypothese <- function(m, index, ML, rules, na_ini){
a <- m[index];
# a2 <- lapply(a, sort);
b <- a[which(duplicated(a) == FALSE)];
b_length <- sapply(b, length);
b_ini    <- 1:length(b);
b2       <- b[order(b_length)];
b_ini    <- b_ini[order(b_length)];
add      <- 1;
ini_new  <- c();

b2_length <- sapply(b2,length);
b.index <- vector("numeric",length=max(b_length)-1);

for( i in 2:max(b_length) - 1){
  b.index[i-1] <- which( b2_length > i)[1];
}
b.index[max(b_length)-1] <- which(b2_length >= max(b_length)-1)[1]

# time <- proc.time()
while(length(b2) > 0){
  ind <- b.index[length(b2[[1]])-1];
  ini <- which(sapply(b2[ind:length(b2)], function(x) {all((b2[[1]]) %in% x)}) == TRUE);
  if(length(ini) >= 1){
    ini_new <- append(ini_new, add);
    b2  <- b2[-1];
    add <- add+1;
  }else{
    b2 <- b2[-1];
    add<-add+1;
  }
}
# proc.time() - time




if(length(ini_new > 0)){
  ini_new <- b_ini[ini_new];
  b <- b[-ini_new];
}

nrow_b   <- length(b);
ncol_b   <- sapply(b, length)
nrow_ML  <- nrow(ML);
ncol_ML  <- ncol(ML);
ML.v     <- as.vector(ML);
hypomass <- sapply(b, function(x) {mean(ML.v[x])})
hypo     <- matrix(NA,ncol=8);
colnames(hypo) <- c("massID","ruleID","nmol","charge","mass","oidscore","ips","massgrp");

for(row in 1:nrow_b){
  for(col in 1:ncol_b[row]){
    adduct <- b[[row]][col] %/% nrow_ML + 1;
    mass   <- b[[row]][col] %% nrow_ML;
    if(mass == 0){
      mass   <- nrow_ML;
      adduct <- adduct - 1;
    }
    hypo   <- rbind(hypo, c(na_ini[mass], adduct, rules[adduct, "nmol"], rules[adduct, "charge"], hypomass[row], rules[adduct, "oidscore"], rules[adduct,"ips"], row));
  }
}
hypo  <- hypo[-1,];
check <- 1;
hypothese <- matrix(NA, ncol=9);
colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "oidscore", "ips", "massgrp", "check");
hypothese <- cbind(hypo, check);

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
    tmpname   <- c();
    tmpnmol   <- c();
    tmpcharge <- 0;
    tmpmass   <- 0;
    tmpips    <- 0;

    #Molekülionen
    if(polarity=="positive"){
        #Wasserstoff, hard codiert
        for(k in 1:mol){
            if(k == 1){
              str    <- "";
              tmpips <- 1;
            }else{
              str    <-  k;
              tmpips <- 0.5;
            };
            name     <- append(name, paste("[", str, "M+H]+", sep=""));
            charge   <- append(charge, 1);
            massdiff <- append(massdiff, 1.0076);
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
            name<-append(name,paste("[M-H-",neutralloss[i,1],"]-",sep=""));
            charge<-append(charge,-1);
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
