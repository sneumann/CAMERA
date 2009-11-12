###Constructor###
setClass("xsAnnotate",
    representation(
                    peaks = "matrix" ,
                    pspectra = "list",
                    isotopes="list",
                    derivativeIons="list",
                    formula="list",
                    sample="numeric",
                    grp_info="numeric",
                    category="character",
                    xcmsSet="xcmsSet",
                    ruleset="data.frame",
                    annoID="matrix",
                    annoGrp="matrix",
                    isoID="matrix",
                    polarity="character"),
    prototype(
                    peaks= matrix(ncol=0,nrow=0),
                    pspectra = list(),
                    isotopes=list(),
                    derivativeIons=list(),
                    formula=list(),
                    sample=NULL,
                    category="",
                    grp_info=NULL,
                    xcmsSet=NULL,
                    ruleset=NULL,
                    annoID=matrix(ncol=3,nrow=0),
                    annoGrp=matrix(ncol=3,nrow=0),
                    isoID=matrix(ncol=4,nrow=0),
                    polarity="")
            );

xsAnnotate <- function(xs=NULL,sample=NA,category=NA){

if(is.null(xs)) { stop("no argument was given"); }
else if(!class(xs)=="xcmsSet") stop("xs is no xcmsSet object ");

object  <- new("xsAnnotate");

if(length(levels(sampclass(xs))) > 1){
    #more than one sample
    if (!nrow(xs@groups) > 0) stop ('First argument must be a xcmsSet with group information or contain only one sample.') #checking alignment
    else if( any(is.na(c(sample,category))) | any(is.null(c(sample,category)))) stop ('For a grouped xcmsSet parameter sample and category must be set and not NULL.')
    else if(! category %in% sampclass(xs)) stop (paste('Paramter category:',category,'is not in',paste(levels(sampclass(xs)),collapse=" ")))
    else if(length(which(sampclass(xs) == category))< sample | sample < 1) stop("Parameter sample must be lower equal than number of samples for the specific category and greater than 0.")
    else {object@sample=sample; object@category=category;}
}else if(length(sampnames(xs)) > 1){
    #one sample category, more than one sample
    if (!nrow(xs@groups) > 0) stop ('First argument must be a xcmsSet with group information or contain only one sample.') #checking alignment
    else if(is.na(sample) | is.null(sample)) stop ('For a grouped xcmsSet parameter sample must be set.')
    else if(length(xs@filepaths) < sample | sample < 1) stop("Parameter sample must be lower equal than number of samples and greater than 0.")
    else { object@sample = sample;object@category = "";}
}else if(length(sampnames(xs)) == 1){
    #only one sample
    object@sample= -1;
    object@category="";
    index= -1;
}else stop("Unknown error with a grouped xcmsSet");

if(object@sample>0){
    #give Information about sample number
    if(object@category==""){
        cat(paste("xcmsSet contains more than one sample, using sample:",sample,"\n"));
        index<-sample;
    }else{
        cat(paste("xcmsSet contains more than one sample, using sample:",sample,"and category:",category,"\n"));
        index<-which(sampclass(xs)==category)[1]+sample-1;
    }
}
#generate Peaklist
tmp <- getPeaks(xs,index);
#remove N/A from peaklist, keep index
if(length(index<-which(is.na(tmp[,1])))>0){
    object@grp_info <-  which(!is.na(tmp[,1]))
    object@peaks    <-  tmp[-index,];
}else{object@peaks  <-  tmp}

#save xcmsSet in the xsAnnotate object
object@xcmsSet  <-  xs;
colnames(object@annoID) <-  c("id","grp_id","rule_id");
colnames(object@annoGrp)<-  c("id","mass","ips");
colnames(object@isoID)  <-  c("mpeak","isopeak","iso","charge")
return(object);
}

setMethod("show", "xsAnnotate", function(object){
#show main information
cat("An \"xsAnnotate\" object!\n");
cat("With",length(object@pspectra),"groups (pseudospectra)\n");
cat("With",nrow(object@peaks),"peaks\n");
if(object@sample>-1){
    if(object@category==""){
        cat(paste("Using sample:",object@sample,"\n"));
    } else { cat(paste("Using sample:",object@sample,"and category:",object@category,"\n"));}
}
memsize <- object.size(object)
cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
})
###End Constructor###

###xsAnnotate generic Methods###
setGeneric("groupFWHM", function(object,sigma=6,perfwhm=0.6) standardGeneric("groupFWHM"))
setMethod("groupFWHM","xsAnnotate", function(object,sigma=6,perfwhm=0.6) {
#Gruppierung nach fwhm
# sigma - number of standard deviation arround the mean (6 = 2 x 3 left and right)
# perfwhm - 0.3;
if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
names<-colnames(object@peaks)
object@pspectra <- list()
if(!"mz" %in% names) stop ("Corrupt peaktable!\n")
if(object@peaks[1,"rt"] == -1) { cat("Warning: no retention times avaiable. Do nothing"); }
else{
    maxo <- object@peaks[,'maxo']; #max intensities of all peaks
    while(!all(is.na(maxo)==TRUE)){
        iint<-which.max(maxo);rtmed<-object@peaks[iint,"rt"]; #highest peak in whole spectra
        rt_min <- object@peaks[iint,"rtmin"];rt_max <- object@peaks[iint,"rtmax"] #begin and end of the highest peak
        hwhm <- ((rt_max-rt_min)/sigma*2.35*perfwhm)/2; #fwhm of the highest peak
        irt<-which(object@peaks[,'rt']>(rtmed-hwhm)&object@peaks[,'rt']<(rtmed+hwhm)&!is.na(maxo)) #all other peaks whose retensiontimes are in the fwhm of the highest peak
        if(length(irt)>0){
            #if peaks are found
            object@pspectra[[length(object@pspectra)+1]]<-irt; #create groups
            maxo[irt]<-NA; #set itensities of peaks to NA, due to not to be found in the next cycle
        }
    }
}
return(invisible(object)); #return object
})

setGeneric("groupCorr",function(object,cor_eic_th=0.75) standardGeneric("groupCorr"));
setMethod("groupCorr","xsAnnotate", function(object,cor_eic_th=0.75) {
if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
#restore xcmsSet from xsAnnotate object
xs<-object@xcmsSet
if (!class(xs)=="xcmsSet")     stop ("xs is not an xcmsSet object")
#check if LC data is available
if(object@peaks[1,"rt"] == -1) {
    cat("Warning: no retention times avaiable. Do nothing\n")
}else {
    if(object@sample== -1){
        tmp <- getAllEICs(xs,1);
    }else if(object@category==""){
        tmp <- getAllEICs(xs,object@sample);
        if(!is.null(object@grp_info)) tmp$EIC <- tmp$EIC[object@grp_info,,,drop=FALSE];
    }else {
        index<-which(sampclass(xs)==object@category)[1]+object@sample-1;
        tmp <- getAllEICs(xs,index);if(!is.null(object@grp_info)) tmp$EIC <- tmp$EIC[object@grp_info,,,drop=FALSE];
    }
    EIC <- tmp$EIC
    scantimes <- tmp$scantimes
    tmp <- calcCL(object,xs, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th)
    object@pspectra <- calc_pc(object@peaks,tmp$CL,tmp$CI)

    ##Workarround: peaks without groups
    peaks<-vector("logical",nrow(object@peaks))
    npspectra<-length(object@pspectra)
    for(i in 1:npspectra){
        peaks[object@pspectra[[i]]]<-TRUE;
    }
    index<-which(peaks==FALSE);
    for(i in 1:length(index)){
        object@pspectra[npspectra+i]<-index[i];
    }
}

return(invisible(object));
})

setGeneric("findIsotopes", function(object,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01) standardGeneric("findIsotopes"));
setMethod("findIsotopes","xsAnnotate", function(object,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01){
if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
#berechne IsotopenMatrix
IM <- calcIsotopes(maxiso=maxiso,maxcharge=maxcharge);
#Normierung
devppm <- ppm / 1000000;
#get mz,rt,into from peaktable
imz <- object@peaks[,"mz"];irt <- object@peaks[,"rt"];mint <- object@peaks[,"into"];
isotope <- vector("list",length(imz));
npspectra <- length(object@pspectra);
isomatrix<-matrix(NA,ncol=5);

#wenn vorher nicht groupFWHM aufgerufen wurde, gruppiere alle Peaks in eine Gruppe
if(npspectra < 1) { npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@peaks));}

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
#             cat(j)
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
                                if (!ISO_RULE1) break
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
if(is.null(nrow(isomatrix))) isomatrix = matrix(isomatrix,byrow=F,ncol=length(isomatrix))
object@isoID<-rbind(object@isoID,isomatrix[,1:4]);
# Zähler für Isotopengruppen
globalcnt <- 0;oldnum<-0;
if(nrow(isomatrix)>0){
  for( i in 1:nrow(isomatrix)){
    if(!isomatrix[i,1]==oldnum){
        globalcnt<-globalcnt+1; isotope[[isomatrix[i,1]]]<-list(y=globalcnt, iso=0, charge= isomatrix[i,4], val=isomatrix[i,5]);oldnum<-isomatrix[i,1];
    }
    isotope[[isomatrix[i,2]]]<-list(y=globalcnt,iso=isomatrix[i,3],charge=isomatrix[i,4],val=isomatrix[i,5]);
  }
}
object@isotopes <- isotope;
return(object);
})

setGeneric("findAdducts",function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL,rules=NULL) standardGeneric("findAdducts"));
setMethod("findAdducts", "xsAnnotate", function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL,rules=NULL){
#Normierung
devppm = ppm / 1000000;
#hole die wichtigen Spalten aus der Peaktable
imz <- object@peaks[,"mz"];irt <- object@peaks[,"rt"];mint <- object@peaks[,"into"];
#anzahl peaks in gruppen für % Anzeige
sum_peaks<-sum(sapply(object@pspectra,length));
#Isotopen
isotopes <- object@isotopes;
#Adduktliste
derivativeIons<-vector("list",length(imz));
#Sonstige Variablen
oidscore<-c();index<-c();
annoID=matrix(ncol=3,nrow=0)
annoGrp=matrix(ncol=3,nrow=0)
colnames(object@annoID) <-  c("id","grp_id","rule_id");
colnames(object@annoGrp)<-  c("id","mass","ips");

if(!(object@polarity=="")){
  cat(paste("polarity is set in xsAnnotate:",object@polarity,"\n"));
  if(!is.null(object@ruleset)){
    rules<-object@ruleset;
  }else stop("ruleset could not read from object!\nFor recalculation set polarity = NULL!\n")
}else {

  #Erkenne polarität
  if(!is.null(polarity)){
      if(polarity %in% c("positive","negative")){
          rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=polarity);
          object@polarity=polarity;
      }else stop("polarity mode unknown, please choose between positive and negative.")
  }else if(length(object@xcmsSet@polarity)>0){
      index<-which(sampclass(object@xcmsSet)==object@category)[1]+object@sample-1
      if(object@xcmsSet@polarity[index] %in% c("positive","negative")){
          rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=object@xcmsSet@polarity[index]);
          object@polarity=polarity;
      }else stop("polarity mode in xcmsSet unknown, please define variable polarity.")
  }else stop("polarity mode could not be estimated from the xcmsSet, please define variable polarity!")
  #save ruleset
  object@ruleset<-rules;
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
if(npspectra < 1){ npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@peaks)); }
cat('\nCalculating possible adducts in',npspectra,'Groups... \n % finished: '); lp <- -1;
#zähler für % Anzeige
npeaks<-0;massgrp<-0;
#für alle Gruppen
for(i in 1:npspectra){
    #       cat(i);
    #Indizes der Peaks in einer Gruppe
    ipeak <- object@pspectra[[i]];
    #Zähler hochzählen und % ausgeben
    npeaks<-npeaks+length(ipeak);
    perc <- round((npeaks) / sum_peaks * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    if (.Platform$OS.type == "windows") flush.console()

    #wenn mehr als ein Peaks in einer Gruppe ist
    if(length(ipeak)>1){
        mz <- imz[ipeak];
        na_ini<-which(!is.na(mz))
        ML <- massDiffMatrix(mz[na_ini],rules)
        m <- fastMatch(as.vector(ML),as.vector(ML),tol = max(2*devppm*mean(mz,na.rm=TRUE))+ mzabs)
        c<-sapply(m,length)
        index<-which(c>=2)
        if(length(index)==0) next;

        #Erstelle Hypothesen
        hypothese<-create_hypothese(m,index,ML,rules,na_ini)
        if(is.null(nrow(hypothese)))next;

        #Entferne Hypothesen, welche gegen Isotopenladungen verstossen!
        if(length(isotopes)>0){
            hypothese <-check_isotopes(hypothese,isotopes,ipeak)
        }
        if(nrow(hypothese)<2){next};

        #Test auf Quasi-Molekülionen
        hypothese <-check_quasimolion(hypothese,quasimolion)
        if(nrow(hypothese)<2){next};

        #Entferne Hypothesen, welche gegen OID-Score&Kausalität verstossen!
        hypothese <- check_oid_causality(hypothese,rules)
        if(nrow(hypothese)<2){next};

        #Prüfe IPS-Score
        hypothese <- check_ips(hypothese)
        if(nrow(hypothese)<2){next};

        #Speichern
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
})


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
  fragment<-rep("",nrow(object@peaks))
  for(i in 1:npspectra){
   print (i);
    index<-object@pspectra[[i]];
    peaktable<-object@peaks[index,];
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
index<-object@pspectra[[grp]];
peaktable<-object@peaks[index,]
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
colnames(peaktable)<-colnames(object@peaks)

return(invisible(data.frame(peaktable,isotopes,adduct,grp,stringsAsFactors=FALSE)));
}

getPeaklist<-function(object){
    xs<-object@xcmsSet
    if (nrow(xs@groups) > 0 && length(sampnames(xs))> 1){
        groupmat <- groups(xs)
        peaklist<- cbind(groupmat,groupval(xs, "medret", "into"))
        cnames <- colnames(peaklist)
        if (cnames[1] == 'mzmed') cnames[1] <- 'mz' else stop ('Peak information ?!?')
        if (cnames[4] == 'rtmed') cnames[4] <- 'rt' else stop ('Peak information ?!?')
        colnames(peaklist) <- cnames
    }else{
        peaklist<-object@peaks;
    }
    adduct<-vector("character",nrow(object@peaks));
    isotopes<-vector("character",nrow(object@peaks));
    pcgroup<-vector("character",nrow(object@peaks));
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
    if(!is.null(object@grp_info)){
        isotopes_tmp<-vector("character",nrow(peaklist));adduct_tmp<-vector("character",nrow(peaklist));pcgroup_tmp<-vector("character",nrow(peaklist));
        for(i in 1:length(isotopes)){
            isotopes_tmp[object@grp_info[i]]<-isotopes[i];
            adduct_tmp[object@grp_info[i]]<-adduct[i];
            pcgroup_tmp[object@grp_info[i]]<-pcgroup[i];
        }
        isotopes<-isotopes_tmp;
        adduct<-adduct_tmp;
        pcgroup<-pcgroup_tmp;
    }
    rownames(peaklist)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
    return(invisible(data.frame(peaklist,isotopes,adduct,pcgroup,stringsAsFactors=FALSE,row.names=NULL)));
}

annotate<-function(xs,sigma=6, perfwhm=0.6,cor_eic_th=0.75,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01,multiplier=3,sample=1,category=NA,polarity="positive"){
    if (!class(xs)=="xcmsSet")     stop ("xs is not an xcmsSet object")
    xs_anno  <- xsAnnotate(xs,sample=sample,category=category);
    xs_anno2 <- groupFWHM(xs_anno,sigma=sigma,perfwhm=perfwhm);
    xs_anno3 <- groupCorr(xs_anno2,cor_eic_th=cor_eic_th);
    xs_anno4 <- findIsotopes(xs_anno3,maxcharge=maxcharge,maxiso=maxiso,ppm=ppm,mzabs=mzabs);
    xs_anno5 <- findAdducts(xs_anno4,multiplier=multiplier,ppm=ppm,mzabs=mzabs,polarity=polarity);
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
            charge<-append(charge,1);massdiff<-append(massdiff,1.0076);nmol<-append(nmol,k);if(k==1){quasi<-append(quasi,1);}else{quasi<-append(quasi,0);};oidscore<-append(oidscore,1);ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M+2H]2+",sep=""));charge<-append(charge,2);massdiff<-append(massdiff,2.0152);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,2);ips<-append(ips,tmpips)
            name<-append(name,paste("[",str,"M+3H]3+",sep=""));charge<-append(charge,3);massdiff<-append(massdiff,3.0228);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,3);ips<-append(ips,tmpips)
            oid<-3;
            for(i in 1:nrow(ionlist)){
                if(ionlist[i,2]<=0)next;
                if(ionlist[i,2]==1){
                    name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]2+",sep=""));
                }else{
                    name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]",ionlist[i,2]+1,"+",sep=""));
                }
                charge	<-	append(charge,ionlist[i,2]+1);
                massdiff<-	append(massdiff,ionlist[i,3]+1.0076);
                nmol	<-	append(nmol,k);
                quasi	<-	append(quasi,0);
                oidscore<-append(oidscore,oid+i);
                if(tmpips>0.75){
                    ips<-append(ips,0.5)
                }else{
                    ips<-append(ips,tmpips);
                }
                #Austausch
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
                nmol	<-append(nmol,k);
                oidscore<-append(oidscore,oid+i)
                if(sum(coeff[i,])==1&& k==1){
                    quasi   <-append(quasi,1);
                    ips     <-append(ips,tmpips);
                }else{
                    quasi	<-append(quasi,0);
                    if(tmpips>0.75){
                        ips	<-append(ips,0.75);
                    }else{
                        ips	<-append(ips,tmpips);
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
                    charge	<-	append(charge,ionlist[index2[ii],2]);
                    massdiff<-	append(massdiff,neutraladdition[i,2]+ionlist[index2[ii],3]);
                    nmol	<-	append(nmol,1);
                    quasi	<-	append(quasi,0);
                    oidscore<-append(oidscore,oid+1);oid<-oid+1;
                    ips<-append(ips,0.5);
                }
            }
            if(neutraladdition[i,1]=="NaCOOH"){next;}
            name<-append(name,paste("[M+H+",neutraladdition[i,1],"]+",sep=""));
            charge<-append(charge,+1);
            massdiff<-	append(massdiff,neutraladdition[i,2]+1.0076);
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
                    if(ionlist[i,2]>1)next;
                    name<-append(name,paste("[",str,"M-2H+",ionlist[i,1],"]-",sep=""));
                    charge	<-	append(charge,ionlist[i,2]-2);
                    massdiff<-	append(massdiff,ionlist[i,3]-(2*1.0076));
                    nmol	<-	append(nmol,k);
                    quasi	<-	append(quasi,0);
                    oidscore<-append(oidscore,oid+i);
                    ips<-append(ips,0.25);
                    next;
                }
                if(ionlist[i,2]== -1){
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]2-",sep=""));
                }else{
                    name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]",ionlist[i,2]+1,"-",sep=""));
                }
                charge	<-	append(charge,ionlist[i,2]-1);
                massdiff<-	append(massdiff,ionlist[i,3]-1.0076);
                nmol	<-	append(nmol,k);
                quasi	<-	append(quasi,0);
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
                nmol	<-append(nmol,k);
                oidscore<-append(oidscore,oid+i)
                if(sum(coeff[i,])==1&& k==1){
                    quasi   <-append(quasi,1);
                }else{
                    quasi	<-append(quasi,0);
                }
                ips	<-append(ips,tmpips);
            }
        }
        oid<-max(oidscore);
        ##Erzeuge Neutral Addition
        index<-which(quasi==1)
        for(i in 1:nrow(neutraladdition)){
            if(length(index2<-which(ionlist[,2]<0))>0){
                for(ii in 1:length(index2)){
                    if(ionlist[index2[ii],2]< -1){
                        name	<-	append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]",abs(ionlist[index2[ii],2]),"-",sep=""));
                    }else{
                        name	<-	append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]-",sep=""));
                    }
                    charge	<-	append(charge,ionlist[index2[ii],2]);
                    massdiff<-	append(massdiff,neutraladdition[i,2]+ionlist[index2[ii],3]);
                    nmol	<-	append(nmol,1);
                    quasi	<-	append(quasi,0);
                    oidscore<-append(oidscore,oid+1);oid<-oid+1;
                    ips<-append(ips,0.5);
                }
            }
            name<-append(name,paste("[M-H+",neutraladdition[i,1],"]-",sep=""));
            charge<-append(charge,-1);
            massdiff<-	append(massdiff,neutraladdition[i,2]-1.0076);
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

combine_xsanno <- function(xsa.pos,xsa.neg,pos=TRUE,tol=2){
# zwei Annotationsobject (pos,neg)
# soll pos zurückgegeben werden? standard= TRUE

##Test Objekte
if (!class(xsa.pos)=="xsAnnotate") stop ("xsa.pos is no xsAnnotate object")
if (!class(xsa.neg)=="xsAnnotate") stop ("xsa.neg is no xsAnnotate object")
if (xsa.pos@polarity!= "positive" & xsa.neg@polarity!= "negative") stop ("xsAnnotate object have wrong polarities.\nOnly pos/neg is allowed!")

##1.Step
rt1 <- sapply(xsa.pos@pspectra, function(x) {mean(xsa.pos@peaks[x,"rt"])})
rt2 <- sapply(xsa.neg@pspectra, function(x) {mean(xsa.neg@peaks[x,"rt"])})
m<-fastMatch(rt1,rt2,tol=tol)
rule_hh<-data.frame("[M-H+H]",1,1,2.0152,1,1,1)
colnames(rule_hh)<-c("name","nmol","charge","massdiff","oidscore","quasi","ips")
endresult<-matrix(ncol=6,nrow=0);
colnames(endresult)<-c("grp_pos","peak_pos","grp_neg","peak_neg","mass","check")

rules.pos<-xsa.pos@ruleset;
rules.neg<-xsa.neg@ruleset;
id_h_pos <- which(rules.pos[,"oidscore"] == 1 & rules.pos[,"quasi"] == 1);
id_h_neg <- which(rules.neg[,"oidscore"] == 1 & rules.neg[,"quasi"] == 1);
annoID.pos<-xsa.pos@annoID;annoID.neg<-xsa.neg@annoID;
annoGrp.pos<-xsa.pos@annoGrp;annoGrp.neg<-xsa.neg@annoGrp;
grp2del<-c();
for(i in 1:length(m)){
#     print(i);
    if(is.null(m[[i]])) next;
    for(j in 1:length(m[[i]]))
    {
#         print(j);
        grp.pos<-getpspectra(xsa.pos,i);
        grp.neg<-getpspectra(xsa.neg,m[[i]][j]);

        #Entferne Isotope aus dem Intensitätsvector, sollen nicht mitannotiert werden
        m.pos<-NA
        for(k in 1:nrow(grp.pos))
        {
         if(!is.null(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]]))
         {
             if(xsa.pos@isotopes[[xsa.pos@pspectra[[i]][k]]]$iso>0){m.pos<-append(m.pos,NA);next}
         }
         m.pos <-append(m.pos,grp.pos[k,1])
        }
        m.pos<-m.pos[-1];
        m.neg<-NA
        for(k in 1:nrow(grp.neg))
        {
         if(!is.null(xsa.neg@isotopes[[xsa.neg@pspectra[[m[[i]][j]]][k]]]))
         {
             if(xsa.neg@isotopes[[ xsa.neg@pspectra[[m[[i]][j]]][k]]]$iso>0){m.neg<-append(m.neg,NA);next}
         }
         m.neg<-append(m.neg,grp.neg[k,1])
        }
        m.neg<-m.neg[-1];
        hypothese.pos<-grp.pos[,"adduct"]
        if(length(index<-which(hypothese.pos != ""))>0)
        {
           liste<-unlist(strsplit(hypothese.pos[index]," "))
	   masslist.pos<-as.numeric(unique(liste[-grep('\\+',liste)]));
	   masslist.pos<-naOmit(masslist.pos)
        }else{masslist.pos<-NULL}
        hypothese.neg<-grp.neg[,"adduct"]
        if(length(index<-which(hypothese.neg != ""))>0)
        {
           liste<-unlist(strsplit(hypothese.neg[index]," "))
	   masslist.neg<-as.numeric(unique(liste[-grep('\\-',liste)]));
	   masslist.neg<-naOmit(masslist.neg)
        }else{masslist.neg<-NULL}
        na_ini1<-which(!is.na(m.pos));
        na_ini2<-which(!is.na(m.neg));
        results<-combine_hypothese(naOmit(m.pos),naOmit(m.neg),rule_hh,na_ini1,na_ini2);
        if(length(results)>0)
        {
            #M+H,M-H gefunden, Pfad: 1
            anno_ids.pos<-xsa.pos@pspectra[[i]][sapply(results,function(x) {x[1]})]
            anno_ids.neg<-xsa.neg@pspectra[[m[[i]][j]]][sapply(results,function(x) {x[2]})]
            for(ii in 1:length(results))
            {
                mass<-mean(c(grp.pos[results[[ii]][1],1],grp.neg[results[[ii]][2],1]));
                if(pos){
                  id.pos<-xsa.pos@pspectra[[i]][results[[ii]][1]];
                  if(length(xsa.pos@derivativeIons[[id.pos]]) >0){ #Pfad: 2
                      if(length(index.pos <- which(annoID.pos[,1]==id.pos & annoID.pos[,3] == id_h_pos))>0){ #Pfad: 3
                          #3.1
                          if(length(index.pos)>1) cat("index.pos are greater than allowed. please debug!");
                          if(!length(del_hypo<-which(annoID.pos[,1]==id.pos & annoID.pos[,3] != id_h_pos))>0) { #Pfad 3.1
                              #Pfad: 3.1.1
                              #bestätige Hypothese
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,2));
                          }else{
                              #Pfad: 3.1.2
                              #delete Hypothese
                              ##TODO: DEL
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,2));
                              for(iii in 1:length(del_hypo)){
                                del_grp<-annoID.pos[del_hypo[iii],2];
                                index2.pos<-which(annoID.pos[,2]==del_grp & annoID.pos[,3] == id_h_pos)
                                if(length(index2.pos)>0){
                                  tmp.pos<-annoID.pos[index2.pos,1]
                                  if(!tmp.pos %in% anno_ids.pos){
                                     grp2del<-append(grp2del,del_grp);
                                    }
                                  }else{ grp2del<-append(grp2del,del_grp);}
                              }
                          }
                      }else{
                        grp_hyp.pos <- annoID.pos[which(annoID.pos[,1]==id.pos),2];
                        index_pfad3_2<-annoID.pos[which(annoID.pos[,2] %in% grp_hyp.pos & annoID.pos[,3] == id_h_pos),1]
                        if(any(index_pfad3_2 %in% anno_ids.pos)) { #Pfad: 3.2
                            #Pfad: 3.2.1
                            next;
                        }else{
                            #Pfad: 3.2.2
                            ##Test auf Grp.
                            if(is.null(masslist.neg)){
                            
                            next;}
                            if(length(unlist(fastMatch(test_mass<-annoGrp.pos[annoID.pos[which(annoID.pos[,1]==id.pos),2],"mass"],masslist.neg,tol=0.05)))>0){
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,4));
                            }else{
                              del_hypo<-which(annoID.pos[,1]==id.pos & annoID.pos[,3] != id_h_pos);
                              del_grp<-annoID.pos[del_hypo,2];
                              grp2del<-append(grp2del,del_grp);
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,3));
                            }
                        }
                      }
                  }else{
                    #2.1
                    endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,1));
                  }
                }else if(!pos){
                  #Join Results to negative List
                  id.neg<-xsa.neg@pspectra[[m[[i]][j]]][results[[ii]][2]];
                  if(length(xsa.neg@derivativeIons[[id.neg]]) >0){ #Pfad: 2
                      if(length(index.neg <- which(annoID.neg[,1]==id.neg & annoID.neg[,3] == id_h_neg))>0){ #Pfad: 3
                          #3.1
                          if(length(index.neg)>1) cat("index.pos are greater than allowed. please debug!");
                          if(!length(del_hypo<-which(annoID.neg[,1]==id.neg & annoID.neg[,3] != id_h_neg))>0) { #Pfad 3.1
                              #Pfad: 3.1.1
                              #bestätige Hypothese
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,2));
                          }else{
                              #Pfad: 3.1.2
                              #delete Hypothese
                              ##TODO: DEL
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,2));
                              for(iii in 1:length(del_hypo)){
                                del_grp<-annoID.neg[del_hypo[iii],2];
                                index2.neg<-which(annoID.neg[,2]==del_grp & annoID.neg[,3] == id_h_neg)
                                if(length(index2.neg)>0){
                                  tmp.neg<-annoID.neg[index2.neg,1]
                                  if(!tmp.neg %in% anno_ids.neg){
                                     grp2del<-append(grp2del,del_grp);
                                    }
                                  }else{ grp2del<-append(grp2del,del_grp);}
                              }
                          }
                      }else{
                        grp_hyp.neg <- annoID.neg[which(annoID.neg[,1]==id.neg),2];
                        index_pfad3_2<-annoID.neg[which(annoID.neg[,2] %in% grp_hyp.neg & annoID.neg[,3] == id_h_neg),1]
                        if(any(index_pfad3_2 %in% anno_ids.neg)) { #Pfad: 3.2
                            #Pfad: 3.2.1
                            next;
                        }else{
                            #Pfad: 3.2.2
                            ##Test auf Grp.
                            if(length(unlist(fastMatch(test_mass<-annoGrp.neg[annoID.neg[which(annoID.neg[,1]==id.neg),2],"mass"],masslist.pos,tol=0.05)))>0){
                              endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,4));
                            }else{
                            del_hypo<-which(annoID.neg[,1]==id.neg & annoID.neg[,3] != id_h_neg);
                            del_grp<-annoID.neg[del_hypo,2];
                            grp2del<-append(grp2del,del_grp);
                            endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,3));
                            }
                        }
                      }
                  }else{
                    #2.1
                    endresult<-rbind(endresult,c(i,results[[ii]][1],m[[i]][j],results[[ii]][2],mass,1));
                  }
                }else cat("Error: pos have wrong value. please debug!\n");
            }
        }
    }
}
if(pos){
index<-which(annoID.pos[,2] %in% grp2del);
if(length(index)>0){
annoID.pos<-annoID.pos[-index,];
}
add_adducts<-vector("character",length(xsa.pos@isotopes));
old_grpid<-max(annoGrp.pos[,1]);
if(nrow(endresult)>0)
{
  for(i in 1:nrow(endresult))
  {
      if(!endresult[i,"check"] %in% c(2,4)){
      old_grpid<-old_grpid+1;
      peakid<-xsa.pos@pspectra[[endresult[i,"grp_pos"]]][endresult[i,"peak_pos"]];
      annoID.pos<-rbind(annoID.pos,c(peakid,old_grpid,id_h_pos))
      annoGrp.pos<-rbind(annoGrp.pos,c(old_grpid,endresult[i,"mass"],2))
      if(endresult[i,"check"]==1){
      add_adducts[peakid]<-paste("Found [M-H]");
      }else{
      add_adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");}
      }else{
          peakid<-xsa.pos@pspectra[[endresult[i,"grp_pos"]]][endresult[i,"peak_pos"]];
          add_adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");}
  }
}
xsa.pos@derivativeIons<-getderivativeIons(annoID.pos,annoGrp.pos,rules.pos,length(xsa.pos@isotopes))
peaklist<-getPeaklist(xsa.pos);
index<-ncol(peaklist)
peaklist<-cbind(peaklist,add_adducts);
colnames(peaklist)[index+1]<-"neg. Mode"
}else if(!pos){
index<-which(annoID.neg[,2] %in% grp2del);
if(length(index)>0){
annoID.neg<-annoID.neg[-index,];
}
add_adducts<-vector("character",length(xsa.neg@isotopes));
old_grpid<-max(annoGrp.neg[,1]);
if(nrow(endresult)>0)
{
  for(i in 1:nrow(endresult))
  {
      if(!endresult[i,"check"] %in% c(2,4)){
      old_grpid<-old_grpid+1;
      peakid<-xsa.neg@pspectra[[endresult[i,"grp_neg"]]][endresult[i,"peak_neg"]];
      annoID.neg<-rbind(annoID.neg,c(peakid,old_grpid,id_h_neg))
      annoGrp.neg<-rbind(annoGrp.neg,c(old_grpid,endresult[i,"mass"],2))
      if(endresult[i,"check"]==1){
      add_adducts[peakid]<-paste("Found [M+H]");
      }else{
      add_adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");}
      }else{
          peakid<-xsa.neg@pspectra[[endresult[i,"grp_neg"]]][endresult[i,"peak_neg"]];
          add_adducts[peakid]<-paste("Verified (",endresult[i,"check"],")",sep="");}
  }
}
xsa.neg@derivativeIons<-getderivativeIons(annoID.neg,annoGrp.neg,rules.neg,length(xsa.neg@isotopes))
peaklist<-getPeaklist(xsa.neg);
index<-ncol(peaklist)
peaklist<-cbind(peaklist,add_adducts);
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

calc_pc <-function(tsa,CL,CI) {

require(RBGL)

gm <- matrix(-1,1,2); colnames(gm) <- c('fromID','toID')
li <- sapply(CL,function(x) length(x) > 0); l <- which(li)
if (!any(li)) return(tsa)

## make connected node matrix
for (j in 1:length(l))
{  gm <- rbind(gm,cbind(l[j],CL[[ l[j]  ]])); }
gm <- gm[-1,]

## create list of connected components
V <- as.character(unique(as.vector(gm)))
OG <- new("graphNEL",nodes=V)
edgemode(OG) <- "undirected"
OG <- addEdge(from=as.character(gm[,"fromID"]),to=as.character(gm[,"toID"]),graph=OG,weights=1)
cc <- connComp(OG)

pci <- matrix(-1,1,2); colnames(pci) <- c('peakID','pcGroupID')

if (length(cc) > 0) {
  NG <- list();NG2<-list();
  cat('Calculating graph cross linking... \n% finished: '); lp <- -1;
  for (cci in 1:length(cc)) {
    perc <- round((cci) / length(cc) * 100)
    if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    if (.Platform$OS.type == "windows") flush.console()
    ## verify all correlation graphs
    G <- subGraph(cc[[cci]],OG)
    if (length(nodes(G)) > 2) {
      ## decomposition might be necessary
# 			G <- removeSelfLoops(G)
        hcs <-  highlyConnSG(G)
        lsg <- sapply(hcs$clusters,function(x) length(x))
        lsg.i <- which(lsg > 1)
        if (length(lsg.i)<1) next;
        for (i in 1:length(lsg.i)) {
          NG[[length(NG)+1]] <- subGraph(hcs$clusters[[lsg.i[i]]],G)
        }
    } else {  NG[[length(NG)+1]] <- G; }
  }
  ## calculate all new pspectra
  for (j in 1:length(NG)){
    NG2[[j]] <- as.numeric(nodes(NG[[j]]))
  }
}
return(NG2)
}

getAllEICs <- function(xs,index=1) {
peaki <- getPeaksIdxCol(xs,NULL);
if(is.matrix(peaki)) peaki<-peaki[,index]
#   nfiles <- length(xs@filepaths)
scantimes <- list()
maxscans <- 0
#   if (nfiles > 1) {
#       cat('Searching maxima .. \n')
#       for (f in 1:nfiles){
#         cat('\Reading raw data file:',xs@filepaths[f])
#         xraw <- xcmsRaw(xs@filepaths[f],profstep=0)
#         cat(',', length(xraw@scantime),'scans. \n')
#         maxscans <- max(maxscans,length(xraw@scantime))
#         scantimes[[f]] <- xraw@scantime
#       }
#
#       for (f in 1:nfiles){
#         if (file.exists(xs@filepaths[f])) {
#           cat('Reading raw data file:',xs@filepaths[f],'\n')
#           xraw <- xcmsRaw(xs@filepaths[f],profstep=0)
#           cat('Generating EIC\'s .. \n')
#           pdata <- as.data.frame(xs@peaks[peaki[,f],]) # data for peaks from file f
#           if (f==1) EIC <- array(integer(0),c(nrow(pdata),maxscans,length(xs@filepaths)))
#           EIC[,,f] <- getEICs(xraw,pdata,maxscans)
#         }
#         else stop('Raw data file:',xs@filepaths[f],' not found ! \n')
#       }
#   }  else { ## create EIC's for single file
    if (file.exists(xs@filepaths[index])) {
        cat('Reading raw data file:',xs@filepaths[index],'\n')
        xraw <- xcmsRaw(xs@filepaths[index],profstep=0)
        cat('Generating EIC\'s .. \n')
        maxscans <- length(xraw@scantime)
        scantimes[[1]] <- xraw@scantime
#           pdata <- as.data.frame(xs@peaks[peaki[,index],])
        pdata <- as.data.frame(xs@peaks[peaki,])
        EIC <- array(integer(0),c(nrow(pdata),maxscans,1))
        EIC[,,1] <- getEICs(xraw,pdata,maxscans)
        }  else stop('Raw data file:',xs@filepaths[index],' not found ! \n')
#   }

#   if (!is.null(file)) save(EIC,scantimes,file=file,compress=TRUE)
#     else
    invisible(list(scantimes=scantimes,EIC=EIC));
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
return(ts)
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

calcCL <-function(object,xs, EIC, scantimes, cor_eic_th){

        CL <- vector("list",nrow(object@peaks))
        CIL <- list()
        ncl<-length(CL);npeaks=0;
        npspectra <- length(object@pspectra);
        cat('Calculating peak correlations... \n% finished: '); lp <- -1;

        #Wenn groupFWHM nicht vorher aufgerufen wurde!
        if(npspectra<1){
          npspectra<-1;object@pspectra[[1]]<-seq(1:nrow(object@peaks));
        }
        for(i in 1:npspectra){
                pi <- object@pspectra[[i]];
                #percent output
                npeaks<-npeaks+length(pi);
                perc <- round((npeaks) / ncl * 100)
                if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
                if (.Platform$OS.type == "windows") flush.console()
                #end percent output
                if(length(pi)>1)
                {
                        for(x in 2:length(pi))
                        {
                                xi <- pi[x];
                                for (y in 1:(x-1))
                                {
                                        yi <- pi[y];
                                        if ( ! (yi %in% CL[[xi]] || yi == xi))
                                        {
# 						cors <- rep(0.0,1)
                                                cors <-0;
                                                ##debug
                                                f <- object@sample; #sample 1
# 						p<-1;
# 						Nf<-1;
                                                ##end debug
                                                eicx <-  EIC[xi,,1]
                                                eicy <-  EIC[yi,,1]
                                                px <- object@peaks[xi,]
                                                py <- object@peaks[yi,]
                                                #No RT workarround
# 						if(px["rtmin"]== -1)
# 						{ rti=1;}
# 						else
# 						{
                                                        crt <- range(px["rtmin"],px["rtmax"],py["rtmin"],py["rtmax"])
                                                        rti <- which(scantimes[[1]] >= crt[1] & scantimes[[1]] <= crt[2])
# 						}
                                                if (length(rti)>1)
                                                {
                                                        dx <- eicx[rti]; dy <- eicy[rti]
                                                        dx[dx==0] <- NA; dy[dy==0] <- NA;
                                                        if (length(which(!is.na(dx) & !is.na(dy))) >= 4)
                                                        {
                                                                ct <- NULL
                                                                options(show.error.messages = FALSE)
                                                                try(ct <- cor.test(dx,dy,method='pearson',use='complete'))
                                                                options(show.error.messages = TRUE)
                                                                if (!is.null(ct) && !is.na(ct))
                                                                {
                                                                        if (ct$p.value <= 0.05) cors <- ct$estimate else cors <- 0;
                                                                }
                                                                else cors <- 0;
                                                        }
                                                        else cors <- 0;
                                                }
                                                else cors <- 0;
# 	       					cot <- which(cors >= cor_eic_th)
# 	       					if (length(cot) / Nf >= gcmfrac)
                                                if(cors>cor_eic_th)
                                                {
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
        if (length(CIL) >0) CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) )
        else return(NULL)
        colnames(CI) <- c('xi','yi','cors')
        return(invisible(list(CL=CL,CI=CI)))
}
###End xsAnnotate intern Function###
