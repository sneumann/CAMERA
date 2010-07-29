##Functions for findIsotopes
calcIsotopes <- function(maxiso,maxcharge){
        M_C12_C13 = 1.0033

        ## Calculate ISO/CHARGE - Matrix
        IM <- matrix(NA,maxiso,maxcharge)
        for (i in 1:maxiso)
        for (j in 1:maxcharge) IM[i,j] <- (i* M_C12_C13) / j

        return(IM)

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