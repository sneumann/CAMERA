##Functions for findAdducts

annotateGrpMPI2 <- function(params){
  library(CAMERA);
  result <- vector(mode="list", length=length(params$i));
  names(result) <- params$i;
  for(ii in 1:length(params$i)){
    res <- CAMERA:::annotateGrp(papply_commondata$pspectra[[params$i[[ii]]]], papply_commondata$imz, 
                                papply_commondata$rules,papply_commondata$mzabs,papply_commondata$devppm,
                                papply_commondata$isotopes, papply_commondata$quasimolion,
                                papply_commondata$rules.idx);
    if(!is.null(res)){
      result[[ii]]<-res;
    }
  }
  return(result);
}

annotateGrpMPI <- function(params, globalParams){
  library(CAMERA);

  result <- vector(mode="list", length=length(params$i));
  names(result) <- params$i;
  for(ii in 1:length(params$i)){
    res <- CAMERA:::annotateGrp(globalParams$pspectra[[params$i[[ii]]]],
                                globalParams$imz,
                                globalParams$rules,
                                globalParams$mzabs,
                                globalParams$devppm,
                                globalParams$isotopes,
                                globalParams$quasimolion,
                                globalParams$rules.idx);
    if(!is.null(res)){
      result[[ii]]<-res;
    }
  }
  return(result);
}

annotateGrp <- function(ipeak, imz, rules, mzabs, devppm, isotopes, quasimolion, rules.idx) {
  #m/z vector for group i with peakindex ipeak
  mz     <- imz[ipeak];
  naIdx <- which(!is.na(mz))

  #Spectrum have only annotated isotope peaks, without monoisotopic peak
  #Give error or warning?
  if(length(na.omit(mz[naIdx])) < 1){
    return(NULL);
  }

  ML <- massDiffMatrix(mz[naIdx], rules[rules.idx,]);
  
  hypothese <- createHypothese(ML, rules[rules.idx, ], devppm, mzabs, naIdx);
  
  #create hypotheses
  if(is.null(nrow(hypothese)) || nrow(hypothese) < 2 ){
    return(NULL);
  }

  #remove hypotheses, which violates via isotope annotation discovered ion charge 
  if(length(isotopes) > 0){
    hypothese <- checkIsotopes(hypothese, isotopes, ipeak);
  }
  
  if(nrow(hypothese) < 2){
    return(NULL);
  }
  
  #Test if hypothese grps include mandatory ions
  #Filter Rules #2
  if(length(quasimolion) > 0){
    hypothese <- checkQuasimolion(hypothese, quasimolion);
  }
  
  if(nrow(hypothese) < 2){
    return(NULL);
  };
  
  #Entferne Hypothesen, welche gegen OID-Score&Kausalitaet verstossen!
  hypothese <- checkOidCausality(hypothese, rules[rules.idx, ]);
  if(nrow(hypothese) < 2){
    return(NULL);
  };
  
  #Pruefe IPS-Score
  hypothese <- checkIps(hypothese)
  if(nrow(hypothese) < 2){
    return(NULL)
  }
  
  #We have hypotheses and want to add neutral losses
  if("typ" %in% colnames(rules)){
    hypothese <- addFragments(hypothese, rules, mz)
  
    hypothese <- resolveFragmentConnections(hypothese)
  }
  return(hypothese);
}


createHypothese <- function(ML, rules, devppm, mzabs, naIdx){
  ML.nrow <- nrow(ML);
  ML.vec <- as.vector(ML);
  max.value <- max(round(ML, 0));
  hashmap <- vector(mode="list", length=max.value);
  for(i in 1:length(ML)){
    val <- trunc(ML[i],0);
    if(val>1){
      hashmap[[val]] <- c(hashmap[[val]],i);
    }
  }
  if("ips" %in% colnames(rules)){
    score <- "ips"
  }else{
    score <- "score"
  }
  hypothese <- matrix(NA,ncol=8,nrow=0);
  colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "score", "massgrp", "check");
  massgrp <- 1;

  for(i in seq(along=hashmap)){
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
    result <- cutree(hclust(dist(candidates)), h=tol);
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
          hypothese <- rbind(hypothese, c(naIdx[mass], adduct, rules[adduct, "nmol"], rules[adduct, "charge"], mean(candidates[m[[ii]]]),  rules[adduct,score],massgrp ,1));
        }
        if(length(unique(hypothese[which(hypothese[, "massgrp"] == massgrp), "massID"])) == 1){
          ##only one mass annotated
          hypothese <- hypothese[-(which(hypothese[,"massgrp"]==massgrp)),,drop=FALSE]
        }else{
          massgrp <- massgrp +1;
        }
    }
  }
  return(hypothese);
}

resolveFragmentConnections <- function(hypothese){
  #Order hypothese after mass
  hypothese <- hypothese[order(hypothese[, "mass"], decreasing=TRUE), ]
  
  for(massgrp in unique(hypothese[, "massgrp"])){
    index <- which(hypothese[, "massgrp"] == massgrp & !is.na(hypothese[, "parent"]))
    if(length(index) > 0) {
      index2 <- which(hypothese[, "massID"] %in% hypothese[index, "massID"] & hypothese[, "massgrp"] != massgrp)
      if(length(index2) > 0){
        massgrp2del <- which(hypothese[, "massgrp"] %in% unique(hypothese[index2, "massgrp"]))
        hypothese <- hypothese[-massgrp2del, ]
      }
    }
  }
  return(hypothese)
}

addFragments <- function(hypothese, rules, mz){
  #check every hypothese grp
  fragments <- rules[which(rules[, "typ"] == "F"), , drop=FALSE]
  hypothese <- cbind(hypothese, NA);
  colnames(hypothese)[ncol(hypothese)] <- "parent"
  if(nrow(fragments) < 1){
    #no fragment exists in rules
    return(hypothese)
  }
  
  orderMZ <- cbind(order(mz),order(order(mz)))
  sortMZ <- cbind(mz,1:length(mz))
  sortMZ <- sortMZ[order(sortMZ[,1]),]
  
  for(massgrp in unique(hypothese[, "massgrp"])){
    for(index in which(hypothese[, "ruleID"] %in% unique(fragments[, "parent"]) & 
      hypothese[, "massgrp"] == massgrp)){
      massID <- hypothese[index, "massID"]
      ruleID <- hypothese[index, "ruleID"]
      indexFrag <- which(fragments[, "parent"] == ruleID)
      
      while(length(massID) > 0){
        result <- fastMatch(sortMZ[1:orderMZ[massID[1],2],1], mz[massID[1]] + 
          fragments[indexFrag, "massdiff"], tol=0.05)
        invisible(sapply(1:orderMZ[massID[1],2], function(x){
          if(!is.null(result[[x]])){
            massID <<- c(massID, orderMZ[x,1]);        
            indexFrags <- indexFrag[result[[x]]];
            tmpRes <- cbind(orderMZ[x,1], as.numeric(rownames(fragments)[indexFrags]), fragments[indexFrags, c("nmol", "charge")],
                                             hypothese[index, "mass"], fragments[indexFrags, c("score")],
                                             massgrp, 1, massID[1], deparse.level=0)
            colnames(tmpRes) <- colnames(hypothese)
            hypothese <<- rbind(hypothese, tmpRes);
          }
        }))
        massID <- massID[-1];
      }    
    }
  }
  return(hypothese)
}


getderivativeIons <- function(annoID, annoGrp, rules, npeaks){
    derivativeIons <- vector("list", npeaks);
    charge <- 0;
    #check that we have annotations
    if(nrow(annoID) < 1){
      return(derivativeIons);
    }
    
    for(i in 1:nrow(annoID)){
        peakid  <-  annoID[i, 1];
        grpid   <-  annoID[i, 2];
        ruleid  <-  annoID[i, 3];
        parent  <-  annoID[i, 4];
#         if(is.null(derivativeIons[[peakid]])){
#           #Peak has no annotation
#           if(charge == 0 | rules[ruleid, "charge"] == charge){
#             derivativeIons[[peakid]][[1]] <- list(rule_id = ruleid, charge = rules[ruleid, "charge"], 
#                                                   nmol= rules[ruleid, "nmol"], name=paste(rules[ruleid, "name"]), 
#                                                   mass=annoGrp[which(annoGrp[, "id"] == grpid), 2],parent=parent)
#             }
#         }else{
#           #Peak has already annotations
          if(charge == 0 | rules[ruleid, "charge"] == charge){
            mass <- annoGrp[which(annoGrp[, "id"] == grpid), 2];
            if(is.na(parent)){
              name <- paste(rules[ruleid, "name"]);
            } else {
              #look for name
              name <- paste(rules[ruleid, "name"]);
              for(ii in seq(along=derivativeIons[[parent]])){
                if(derivativeIons[[parent]][[ii]]$mass == mass){
                  
                  break;
                }                  
              }
            }
            derivativeIons[[peakid]][[(length(derivativeIons[[peakid]])+1)]] <- 
              list(rule_id = ruleid, charge=rules[ruleid, "charge"], 
                   nmol=rules[ruleid, "nmol"], name=name, 
                   mass=mass, 
                   parent=parent)
          }
#         }
        charge=0;
    }
  return(derivativeIons);
}

checkIps <- function(hypothese){
  for(hyp in 1:nrow(hypothese)){
    if(length(which(hypothese[, "massgrp"] == hypothese[hyp, "massgrp"])) < 2){
      hypothese[hyp, "check"] = 0;
    }
  }
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ];
  if(is.null(nrow(hypothese))) {
    hypothese <- matrix(hypothese, byrow=F, ncol=9)
  }
  if(nrow(hypothese) < 1){
    colnames(hypothese)<-c("massID", "ruleID", "nmol", "charge", "mass", "oidscore", "ips","massgrp", "check")
    return(hypothese)
  }
  for(hyp in 1:nrow(hypothese)){
    if(length(id <- which(hypothese[, "massID"] == hypothese[hyp, "massID"] & hypothese[, "check"] != 0)) > 1){
      masses <- hypothese[id, "mass"]
      nmasses <- sapply(masses, function(x) { 
                                  sum(hypothese[which(hypothese[, "mass"] == x), "score"]) 
                                  })
      masses <- masses[-which(nmasses == max(nmasses))];
      if(length(masses) > 0){
        hypothese[unlist(sapply(masses, function(x) {which(hypothese[, "mass"]==x)})), "check"]=0;
      }
    }
  }
  
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
  return(hypothese)
}

checkOidCausality <- function(hypothese,rules){
  #check every hypothese grp
  for(hyp in unique(hypothese[,"massgrp"])){
    hyp.nmol <- which(hypothese[, "massgrp"] == hyp & hypothese[, "nmol"] > 1)

    for(hyp.nmol.idx in hyp.nmol){
      if(length(indi <- which(hypothese[, "mass"] == hypothese[hyp.nmol.idx, "mass"] & 
        abs(hypothese[, "charge"]) == hypothese[, "nmol"])) > 1){
        if(hyp.nmol.idx %in% indi){
          #check if [M+H] [2M+2H]... annotate the same molecule
          massdiff <- rules[hypothese[indi, "ruleID"], "massdiff"] / 
            rules[hypothese[indi, "ruleID"], "charge"]
          if(length(indi_new <- which(duplicated(massdiff))) > 0){
            hypothese[hyp.nmol.idx, "check"] <- 0;
          }
        }
      }
    }
  }
    
#     #check nmol
#     if(hypothese[hyp, "nmol"] > 1){
#       #nmol > 1;
#       checkSure <- TRUE;
#       if(hypothese[hyp, "charge"] == 1){
#         #nmol > 1 and charge = 1; e.g. [2M+H]+, ensure [M+H] is there
#         if(length(which(hypothese[, "mass"] == hypothese[hyp, "mass"] & hypothese[, "oidscore"] == hypothese[hyp, "oidscore"])) > 1){
#           #same oidscore is there, could also be [3M+H]; otherwise could not check
#           for(prof in (hypothese[hyp, "nmol"] - 1):1){
#           #check if [M+H] is there, for a [3M+H], [2M+H] and [M+H] has to be there
#             indi <- which(hypothese[,"mass"] == hypothese[hyp,"mass"] & hypothese[,"oidscore"] == hypothese[hyp,"oidscore"] & hypothese[,"nmol"] == prof)
#             if(length(indi) == 0){
#               checkSure <- FALSE;
#               hypothese[hyp,"check"] <- 0;
#               next;
#             }
#           }
#         }
#       }else if(abs(hypothese[hyp, "charge"]) == hypothese[hyp, "nmol"]){
#         #nmol > 1 and charge = nmol; e.g. [2M+2H]2+
#         if(length(which(hypothese[, "mass"] == hypothese[hyp, "mass"] & hypothese[, "oidscore"] == hypothese[hyp, "oidscore"])) > 1){
#           for(prof in (hypothese[hyp,"nmol"]-1):1){
#             indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
#             if(length(indi) == 0){
#               checkSure <- FALSE;
#               hypothese[hyp,"check"] <- 0;#next;
#             }
#           }
#         }
#         if(length(indi <- which(hypothese[, "mass"] == hypothese[hyp, "mass"] & abs(hypothese[, "charge"]) == hypothese[, "nmol"])) > 1){
#           #check if [M+H] [2M+2H]... annotate the same molecule
#           massdiff <- rules[hypothese[indi, "ruleID"], "massdiff"] / rules[hypothese[indi, "ruleID"], "charge"]
#           if(length(indi_new <- which(duplicated(massdiff))) > 0){
#             checkSure <- FALSE;
#             hypothese[hyp, "check"] <- 0;
#           }
#         }
#       }
#       if(checkSure){
#         hypothese[hyp, "check"] <- 1;
#       }
#     }
#   }
  hypothese <- hypothese[which(hypothese[, "check"] == TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
  return(hypothese)
}

checkQuasimolion <- function(hypothese, quasimolion){
  hypomass <- unique(hypothese[, "mass"])
  for(mass in 1:length(hypomass)){
    if(!any(quasimolion %in% hypothese[which(hypothese[, "mass"] == hypomass[mass]), "ruleID"])){
      hypothese[which(hypothese[, "mass"] == hypomass[mass]), "check"] = 0;
    }else if(is.null(nrow(hypothese[which(hypothese[, "mass"] == hypomass[mass]), ]))){
      hypothese[which(hypothese[, "mass"] == hypomass[mass]), "check"] = 0;
    }
  }
  
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), , drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
    
  return(hypothese)
}

checkIsotopes <- function(hypothese, isotopes, ipeak){
  for(hyp in 1:nrow(hypothese)){
    peakid <- ipeak[hypothese[hyp, 1]];
    if(!is.null(isotopes[[peakid]])){
      #Isotope da
      explainable <- FALSE;
      if(isotopes[[peakid]]$charge == abs(hypothese[hyp, "charge"])){
        explainable <- TRUE;
      }
      if(!explainable){
        #delete Rule
        hypothese[hyp,"check"]=0;
      }
    }
  }
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
  
  return(hypothese)
}

checkHypothese <- function(hypothese){
  if(is.null(nrow(hypothese))){
    hypothese <- matrix(hypothese, byrow=F, ncol=8)
  } 
  colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "score", "massgrp", "check")
  for(i in unique(hypothese[,"massgrp"])){
   if(length(unique(hypothese[which(hypothese[, "massgrp"] == i), "massID"])) == 1){
    ##only one mass annotated
    hypothese <- hypothese[-(which(hypothese[,"massgrp"]==i)), , drop=FALSE]
   } 
  }
  return(hypothese)
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


massDiffMatrix <- function(m, rules){
  #m - m/z vector
  #rules - annotation rules
  nRules <- nrow(rules);
  DM   <- matrix(NA, length(m), nRules)
  
  for (i in seq_along(m)){
    for (j in seq_len(nRules)){
      DM[i, j] <- (abs(rules[j, "charge"] * m[i]) - rules[j, "massdiff"]) / rules[j, "nmol"]    # ((z*m) - add) /n
    }
  }
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
