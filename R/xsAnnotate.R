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
				xcmsSet="xcmsSet"),
	
	 prototype(
				peaks= matrix(ncol=0,nrow=0),
				pspectra = list(),
				isotopes=list(),
				derivativeIons=list(),
				formula=list(),
				sample=NULL,
				category="",
				grp_info=NULL,
				xcmsSet=NULL)
	);

xsAnnotate <- function(xs=NULL,sample=NA,category=NA){
 
 if(is.null(xs)) { stop("no argument was given"); }
 else if(!class(xs)=="xcmsSet") stop("xs is no xcmsSet object ");

 object		<- new("xsAnnotate");
if(length(levels(sampclass(xs))) > 1)
{
	if (!nrow(xs@groups) > 0) stop ('First argument must be a xcmsSet with group information or contain only one sample.') #checking alignment
	else if( any(is.na(c(sample,category))) | any(is.null(c(sample,category)))) stop ('For a grouped xcmsSet parameter sample and category must be set and not NULL.')
	else if(! category %in% sampclass(xs)) stop (paste('Paramter category:',category,'is not in',paste(levels(sampclass(xs)),collapse=" ")))
	else if(length(which(sampclass(xs) == category))< sample | sample < 1) stop("Parameter sample must be lower equal than number of samples for the specific category and greater than 0.")
	else {object@sample=sample; object@category=category;}
}else if(length(sampnames(xs)) > 1)
{
	#one sample category, more than one sample
	if (!nrow(xs@groups) > 0) stop ('First argument must be a xcmsSet with group information or contain only one sample.') #checking alignment
	else if(is.na(sample) | is.null(sample)) stop ('For a grouped xcmsSet parameter sample must be set.')
	else if(length(xs@filepaths) < sample | sample < 1) stop("Parameter sample must be lower equal than number of samples and greater than 0.")
	else { object@sample = sample;object@category = "";}
}else if(length(sampnames(xs)) == 1)
{
	#juhu einfacher Fall
	object@sample= -1;
	object@category="";
	index= -1;
}else stop("Unknown error with a grouped xcmsSet");


if(object@sample>0)
{
 if(object@category=="")
 {
	cat(paste("xcmsSet contains more than one sample, using sample:",sample,"\n"));
	index<-sample;
 }else{
	cat(paste("xcmsSet contains more than one sample, using sample:",sample,"and category:",category,"\n"));
	index<-which(sampclass(xs)==category)[1]+sample-1;	
 }
}

tmp <- getPeaks(xs,index);
if(length(index<-which(is.na(tmp[,1])))>0)
{
	object@grp_info<-which(!is.na(tmp[,1]))
	object@peaks	<- tmp[-index,];
}else{object@peaks	<- tmp}
object@xcmsSet	<- xs;

##Negative Warning.



return(object);
}

setMethod("show", "xsAnnotate", function(object) {

	cat("An \"xsAnnotate\" object!\n");
	cat("With",length(object@pspectra),"groups (pseudospectra)\n");
	cat("With",nrow(object@peaks),"peaks\n");
	if(object@sample>-1)
	{
   		if(object@category=="")
 		{ cat(paste("Using sample:",object@sample,"\n")); }
 		else { cat(paste("Using sample:",object@sample,"and category:",object@category,"\n")); }
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
if(!"mz" %in% names) stop ("Corrupt peaktable!\n")
if(object@peaks[1,"rt"] == -1) { cat("Warning: no retention times avaiable. Do nothing"); }
else
{
	maxo <- object@peaks[,'maxo']; #max intensities of all peaks
	while(!all(is.na(maxo)==TRUE))
	{
		iint<-which.max(maxo);rtmed<-object@peaks[iint,"rt"]; #highest peak in whole spectra
		rt_min <- object@peaks[iint,"rtmin"];rt_max <- object@peaks[iint,"rtmax"] #begin and end of the highest peak
# 		sd_est<- (rt_max-rt_min)/sigma;
# 		hwhm<-sd_est*sqrt( 2*log(2) - 2*log(1/(sd_est*sqrt(2*pi)))*object@peaks[iint,"maxo"])
		hwhm <- ((rt_max-rt_min)/sigma*2.35*perfwhm)/2; #fwhm of the highest peak
		irt<-which(object@peaks[,'rt']>(rtmed-hwhm)&object@peaks[,'rt']<(rtmed+hwhm)&!is.na(maxo)) #all other peaks whose retensiontimes are in the fwhm of the highest peak
		if(length(irt)>0) #if peaks are found
		{
			object@pspectra[[length(object@pspectra)+1]]<-irt; #create groups
			maxo[irt]<-NA; #set itensities of peaks to NA, due to not to be found in the next cycle
		}
	}
}
return(invisible(object)); #return object
})

setGeneric("groupCorr",function(object,cor_eic_th=0.75) standardGeneric("groupCorr"));
setMethod("groupCorr","xsAnnotate", function(object,cor_eic_th=0.75) {
xs<-object@xcmsSet
#Peak Correlation
if (!class(object)=="xsAnnotate") stop ("no xsAnnotate object")
# if(!exists(as.character(substitute(xs)))) stop ("xcmsSet object cannot be found")
if (!class(xs)=="xcmsSet")     stop ("xs is not an xcmsSet object")
if(object@peaks[1,"rt"] == -1) { cat("Warning: no retention times avaiable. Do nothing") }
else
{
  if(object@sample== -1)
  {
	tmp <- getAllEICs(xs,1);
  }else if(object@category=="")
  {
	tmp <- getAllEICs(xs,object@sample);
	if(!is.null(object@grp_info))tmp$EIC <- tmp$EIC[object@grp_info,,,drop=FALSE];
  }else
  {
  	index<-which(sampclass(xs)==object@category)[1]+object@sample-1;
	tmp <- getAllEICs(xs,index);if(!is.null(object@grp_info)) tmp$EIC <- tmp$EIC[object@grp_info,,,drop=FALSE];
  }
  EIC <- tmp$EIC
  scantimes <- tmp$scantimes
  rm(tmp);gc();

  tmp <- calcCL(object,xs, EIC=EIC, scantimes=scantimes, cor_eic_th=cor_eic_th)
  object@pspectra <- calc_pc(object@peaks,tmp$CL,tmp$CI)

  rm(tmp);gc();
}
  return(invisible(object));
})

setGeneric("findIsotopes", function(object,...) standardGeneric("findIsotopes"));
setMethod("findIsotopes","xsAnnotate", function(object,maxcharge=3,maxiso=4,ppm=5,mzabs=0.01){

#lade library
require(Hmisc)

#berechne IsotopenMatrix
IM <- calcIsotopes(maxiso=maxiso,maxcharge=maxcharge);

#Normierung
devppm <- ppm / 1000000;

#hole die wichtigen Spalten aus der Peaktable
imz <- object@peaks[,"mz"];irt <- object@peaks[,"rt"];mint <- object@peaks[,"into"];
isotope <- vector("list",length(imz)); 

npspectra <- length(object@pspectra);
isomatrix<-matrix(NA,ncol=5);

#wenn vorher nicht groupFWHM aufgerufen wurde, gruppiere alle Peaks in eine Gruppe
if(npspectra < 1)
{ npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@peaks)); }

#Suche Isotope in jeder Gruppe
for(i in 1:npspectra)
{
	#indizes der peaks aus der gruppe in der peaktable
	ipeak <- object@pspectra[[i]];
	
	#hat gruppe mehr als einen Peak, sonst mach nichts
	if(length(ipeak)>1)
	{
		#masse und intensität der Peaks
		mz <- imz[ipeak];int <- mint[ipeak];

		#matrix der peaks mit allen wichtigen Informationen
		spectra<-matrix(c(mz,int,ipeak),ncol=3)
		spectra<-spectra[order(spectra[,1]),];
		cnt<-nrow(spectra);

		#für jeden Peak
		for ( j in 1:(length(mz)-1))
		{
			#erzeuge Differenzmatrix
			MI <- spectra[j:cnt,1] - spectra[j,1];

			#für alle erlaubte Ladungen
			for( charge in maxcharge:1)
			{
				#Suche Übereinstimmungen der Isotopenabständen mit der Differenzmatrix
				m <- find.matches(MI,IM[,charge],tol= max(2*devppm*mz)+ mzabs,scale=1)
				if (!is.null(dim(m))) stop('--------??') #Unbekannter Fehler

				#Für jeden Match, teste welches Isotope gefunden wurde
				if(any(m$matches>0))
				{
					#für alle erlaubten Isotopenpeaks
					for( iso in 1:maxiso)
					{
						#wurde der iso-Isotopenpeak gefunden?
						pos <- which(m$matches == iso)
						if (isTRUE(pos > 0))
						{ # Isotop Nr. iso scheint zu existieren
                    					dev <- (devppm * spectra[pos+j-1,1]) + (devppm + spectra[j,1])
                    					if (isTRUE(all.equal(spectra[pos+j-1,1],spectra[j,1] + IM[iso,charge] ,tolerance=(dev + mzabs),scale=1)))
							{ # Isotop Nr. iso existiert
                       						int.available <- all(!is.na(c(spectra[pos+j-1,2],spectra[j,2])))
                       						if (iso == 1)
								{ #wenn der erste Isotopenpeak gefunden wurde
                          						if (int.available) 
                            						ISO_RULE1 <- (spectra[pos+j-1,2] < spectra[j,2] ) ## isotopic rule 
                              						else ISO_RULE1 <- TRUE
                          						ISO_RULE <- TRUE
                          					}
                         					else
								{
									# Sind alle anderen isotopen Peaks da?
									test<-match(apply(isomatrix[,c(1,3,4)],1,function(x) {paste(x,collapse=" ")}),apply(matrix(cbind(spectra[j,3],1:(iso-1),charge),ncol=3),1, function(x) {paste(x,collapse=" ")}))
									if(length(na.omit(test))==(iso-1))
									{
                          							ISO_RULE1 <- TRUE 
	                          						if (int.available) ISO_RULE <- (spectra[pos+j-1,2] < spectra[j,2]) 
        	                    						else ISO_RULE <- TRUE
									}
									else
									{
										ISO_RULE1 <- FALSE
									}
                         					}
                       						if (!ISO_RULE1) break 
                       						if (ISO_RULE1 && ISO_RULE)
								{
									#Neues Isotope gefunden
									val=FALSE; #TODO: Was macht das hier?
									#TODO: Intrinsische Ladungen betrachten
									if(!length(which(isomatrix[,1]==spectra[j,3] & isomatrix[,2]==spectra[pos+j-1,3]))>0)
									{
										if(!length(which(isomatrix[,2]==spectra[j,3])>0))
										{
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

# Zähler für Isotopengruppen
globalcnt <- 0;oldnum<-0;
for( i in 1:nrow(isomatrix))
{
	if(!isomatrix[i,1]==oldnum) {
	globalcnt<-globalcnt+1; isotope[[isomatrix[i,1]]]<-list(y=globalcnt, iso=0, charge= isomatrix[i,4], val=isomatrix[i,5]);oldnum<-isomatrix[i,1];}
	isotope[[isomatrix[i,2]]]<-list(y=globalcnt,iso=isomatrix[i,3],charge=isomatrix[i,4],val=isomatrix[i,5]);
}
# for(i in 1:nrow(isomatrix))
# {
# 	if(is.null(isotope[[isomatrix[i,1]]]))
# 	{
# 		isotope[[isomatrix[i,1]]][[1]]<-list( y = isomatrix[i,2], iso=isomatrix[i,3], charge= isomatrix[i,4], val=isomatrix[i,5])
# 	}
# 	else
# 	{
# 		isotope[[isomatrix[i,1]]][[(length(isotope[[isomatrix[i,1]]])+1)]] <- list(y = isomatrix[i,2], iso=isomatrix[i,3], charge= isomatrix[i,4], val=isomatrix[i,5]) 
# 	}
# }
object@isotopes <- isotope;
return(object);
})

setGeneric("findAdducts",function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL) standardGeneric("findAdducts"));
setMethod("findAdducts", "xsAnnotate", function(object,ppm=5,mzabs=0.015,multiplier=3,polarity=NULL){
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

#Erkenne polarität
if(!is.null(polarity))
{
	if(polarity %in% c("positive","negative"))
	{
		rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=polarity);
	}else stop("polarity mode unknown, please choose between positive and negative.")
}else if(length(object@xcmsSet@polarity)>0){
	index<-which(sampclass(object@xcmsSet)==object@category)[1]+object@sample-1
	if(object@xcmsSet@polarity[index] %in% c("positive","negative"))
	{
		rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=object@xcmsSet@polarity[index]);
	}else stop("polarity mode in xcmsSet unknown, please define variable polarity.")
}else stop("polarity mode couldn't be estimated from the xcmsSet, please define variable polarity!")



#Einlesen der Regeln
# rules<-calcRules(maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2)
quasimolion<-which(rules[,"quasi"]==1)

#Entferne Isotope aus dem Intensitätsvector, sollen nicht mitannotiert werden
if(length(isotopes)>0){

	for(x in 1:length(isotopes))
	{
		if(!is.null(isotopes[[x]]))
		{
			if(isotopes[[x]]$iso!=0)imz[x]=NA;
		}
	}
}

#Anzahl Gruppen
npspectra <- length(object@pspectra);

#Wenn vorher nicht gruppiert wurde, alle Peaks in eine Gruppe stecken
if(npspectra < 1)
{ npspectra <- 1;object@pspectra[[1]]<-seq(1:nrow(object@peaks)); }

cat('\n Calculating possible adducts in',npspectra,'Groups... \n % finished: '); lp <- -1;
#zähler für % Anzeige
npeaks<-0;

#für alle Gruppen
for(i in 1:npspectra)
{
# 	cat(i);
	#Indizes der Peaks in einer Gruppe
	ipeak <- object@pspectra[[i]];
	#Zähler hochzählen und % ausgeben
	npeaks<-npeaks+length(ipeak);
	perc <- round((npeaks) / sum_peaks * 100)
    	if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    	if (.Platform$OS.type == "windows") flush.console()

	#wenn mehr als ein Peaks in einer Gruppe ist
	if(length(ipeak)>1)
	{
		mz <- imz[ipeak];
		ML <- massDiffMatrix(mz,rules)
		m <- find.matches(as.vector(ML),as.vector(ML),tol= max(2*devppm*mean(mz,na.rm=TRUE))+ mzabs,scale=1,maxmatch=20)
		if (is.null(dim(m$matches)))next;
# 		m <- find.matches(as.vector(ML),asm.vector(ML),tol= 0.005 + devgen,scale=1,maxmatch=20)
# 		index<-which(m$matches[,ncol(m$matches)]>0)
		c<-apply(m$matches,1,function(x) {length(which(na.omit(x)>0))})
		#Kausalität, man geht danach wo mehr Hypothesen gefunden werden
 		index<-which(c>=4)
		if(length(index)==0) next;
 		a<-m$matches[index,];
		a[a==0]<-NA;
		if(is.null(nrow(a))) a <- matrix(a,byrow=F,ncol=20)
		a<-t(apply(a,1,function(x) {sort(x,na.last=TRUE)}));
		b<-a[which(duplicated(a)==FALSE),];
		if(is.null(nrow(b))) b = matrix(b,byrow=F,ncol=ncol(a)) 
# 		ini<-which.max(apply(b,1,function(x){length(na.omit(x))}))
		ini_old<-c();
		##TODO: Nötig?
		for(z in 1:nrow(b))
		{
			index<-which(apply(b,1,function(x) {all(na.omit(b[z,]) %in% x)})==TRUE)
			index<-index[which(!index==z)];
			if(length(index)>0)
			{
				if(length(na.omit(b[index,]))>length(na.omit(b[z,])))
				{
					ini_old<-append(ini_old,z);
				}
			}
		}
		if(length(ini_old>0))b<-b[-ini_old,]
		if(is.null(nrow(b))) b = matrix(b,byrow=F,ncol=ncol(a))

		nrow_b<-nrow(b);
#  		ncol_b<-ncol(length(b));
		ncol_b<-apply(b,1,function(x) {length(na.omit(x))})
		nrow_ML<-nrow(ML);
 		ncol_ML<-ncol(ML);
		c<-as.vector(ML);
		hypomass<-apply(b,1,function(x) {mean(c[x],na.rm = TRUE)})
		hypo<-matrix(NA,ncol=7)
		colnames(hypo)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips")
		for(row in 1:nrow_b)
		{
			for(col in 1:ncol_b[row])
			{	
				adduct<-b[row,col]%/%nrow_ML+1;
				mass <- b[row,col]%%nrow_ML;if(mass==0){mass<-nrow_ML;adduct<-adduct-1;}
				hypo <- rbind(hypo,c(mass,adduct,rules[adduct,"nmol"],rules[adduct,"charge"],hypomass[row],rules[adduct,"oidscore"],rules[adduct,"ips"]));
			}
		}
		hypo<-hypo[-1,]
		check<-1
		hypothese<-matrix(NA,ncol=8)
		colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","check")
		hypothese<-cbind(hypo,check);
		if(is.null(nrow(hypothese)))next;

		#Entferne Hypothesen, welche gegen Isotopenladungen verstossen!
		for(hyp in 1:nrow(hypothese))
		{
			peakid<-ipeak[hypothese[hyp,1]];
			if(!is.null(isotopes[[peakid]]))
			{
				#Isotope da
				explainable<-FALSE;
				if(isotopes[[peakid]]$charge==abs(hypothese[hyp,"charge"]))
				{
					explainable<-TRUE;
				}
				if(!explainable)
				{
						#delete Rule
# 						object@derivativeIons[[ii]]<-NULL;
					hypothese[hyp,"check"]=0;
					##Änderung weil, sonst Regeln verloren gehen
# 					ihypo <- which(hypothese[,"massID"]==hypothese[hyp,"massID"] & hypothese[,"check"]!=0 & hypothese[,"mass"]==hypothese[hyp,"mass"])
# 					if(length(ihypo)<1)
# 					{
# 						ihypo<-which(hypothese[,"mass"]==hypothese[hyp,"mass"])
# 						if(!(is.null(ihypo)))
# 						{
# 							hypothese[ihypo,"check"]=0;
# 						}
# 					}
				}
			}
		}
		hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
		if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=8)
		colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","check")
		if(nrow(hypothese)<2){next};

		##Test auf Quasi-Molekülionen
		hypomass<-unique(hypothese[,"mass"])
		for(mass in 1:length(hypomass))
		{
			if(!any(quasimolion %in% hypothese[which(hypothese[,"mass"]==hypomass[mass]),"ruleID"]))
			{
				hypothese[which(hypothese[,"mass"]==hypomass[mass]),"check"]=0;
			}
			else if(is.null(nrow(hypothese[which(hypothese[,"mass"]==hypomass[mass]),])))
			{
				hypothese[which(hypothese[,"mass"]==hypomass[mass]),"check"]=0;
			}
		}
		hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
		if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=8)
		colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","check")
		if(nrow(hypothese)<2){next};

		#Entferne Hypothesen, welche gegen OID-Score&Kausalität verstossen!
		for(hyp in 1:nrow(hypothese))
# 		for(hyp in 12:12)
		{
			#check nmol
			if(hypothese[hyp,"nmol"]>1)
			{
				check_sure<-TRUE;
				if(hypothese[hyp,"charge"]==1)
				{
					if(length(which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"]))>1){
					for(prof in (hypothese[hyp,"nmol"]-1):1)
					{
						indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
						if(length(indi)==0)
						{
							check_sure<-FALSE;
							hypothese[hyp,"check"]<-0;next;
						}
					}
					}
				}
				else if(hypothese[hyp,"charge"]==hypothese[hyp,"nmol"])
				{
					if(length(which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"]))>1)
					{
						for(prof in (hypo[hyp,"nmol"]-1):1)
						{
							indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
# 							indi<-which(rules[,"massdiff"]==rules[hypothese[hyp,"ruleID"],"massdiff"] && rules[,"nmol"]==prof);
							if(length(indi)==0)
							{
								check_sure<-FALSE;
								hypothese[hyp,"check"]<-0;#next;
							}
						}
					}
					if(length(indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"charge"]==hypothese[,"nmol"]))>1)
					{
						massdiff<-rules[hypothese[indi,"ruleID"],"massdiff"]/rules[hypothese[indi,"ruleID"],"charge"]
						if(length(indi_new<-which(duplicated(massdiff)))>0)
						{
							check_sure<-FALSE;
							hypothese[indi[indi_new],"check"]<-0;
						}
					}
				}
				if(check_sure){hypothese[hyp,"check"]<-1;}
			}else
			{
				if(hypo[hyp,"charge"]>1)
				{
					##todo
# 					check[blub]<-TRUE;
				}else
				{
					#nothing to say
# 					check[blub]<-TRUE;
				}
			}
		}
		hypothese<-hypothese[which(hypothese[,"check"]==TRUE),];
		charge=0;
		if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=8)
		colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","check");
		if(nrow(hypothese)<2){next};

		#Lösche xM+xIon
		

		##Schritt?? Score proof
		for(hyp in 1:nrow(hypothese))
# 		for(hyp in 18:24)
		{
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
		if(is.null(nrow(hypothese))) hypothese = matrix(hypothese,byrow=F,ncol=8)
		colnames(hypothese)<-c("massID","ruleID","nmol","charge","mass","oidscore","ips","check")
		if(nrow(hypothese)<2){next};

		for(hyp in 1:nrow(hypothese))
		{
			peakid<-ipeak[hypothese[hyp,"massID"]];
			if(is.null(derivativeIons[[peakid]]))
			{
				if(charge==0 || rules[hypothese[hyp,3],"charge"]==charge)
				{
				 derivativeIons[[peakid]][[1]]<- list( rule_id = hypothese[hyp,2], charge=rules[hypothese[hyp,3],"charge"], nmol= rules[hypothese[hyp,3],"nmol"], name=paste(rules[hypothese[hyp,2],"name"]),mass=hypothese[hyp,"mass"])
				}
			}
			else
			{
				if(charge==0 || rules[hypothese[hyp,3],"charge"]==charge)
				{
				derivativeIons[[peakid]][[(length(derivativeIons[[peakid]])+1)]] <- list( rule_id = hypothese[hyp,2], charge=rules[hypothese[hyp,3],"charge"], nmol= rules[hypothese[hyp,3],"nmol"], name=paste(rules[hypothese[hyp,2],"name"]),mass=hypothese[hyp,"mass"])
				}
			}
			charge=0;
		}
	}
}
cat("\n");
object@derivativeIons <- derivativeIons;
return(object)
})
###End xsAnnotate generic Methods###

###xsAnnotate exported Methods###
getGroup <- function(object,grp){
	index<-object@pspectra[[grp]];
	peaktable<-object@peaks[index,]
	adduct<-vector("character",length(index));
	isotopes<-vector("character",length(index));
	ions<-object@derivativeIons;
	iso<-object@isotopes;
	lions <- ions[index];
	liso  <- iso[index];
 
	for(i in 1:length(lions))
	{
		if(!is.null(lions[[i]]))
		{
			if(length(lions[[i]])>1)
			{	
				names<-c();
				for(ii in 1:length(lions[[i]]))
				{
					names<-paste(names,lions[[i]][[ii]]$name,lions[[i]][[ii]]$mass);
				}
				adduct[i]<-names;
			}
			else
			{
				adduct[i]<-paste(lions[[i]][[1]]$name,lions[[i]][[1]]$mass);
			}
		}
		if(!is.null(liso[[i]]))
		{
# 			if(length(liso[[i]])>1)
# 			{	
# 				names<-c();
# 				for(ii in 1:length(liso[[i]]))
# 				{
# 					names<-paste(names,"[",liso[[i]][[ii]]$y,"] ",liso[[i]][[ii]]$iso," ",liso[[i]][[ii]]$charge,"+ ",sep="");
# 				}
# 				isotopes[i]<-names;
# 			}
# 			else
# 			{
				isotopes[i]<-paste("[",liso[[i]]$y,"] ",liso[[i]]$iso," ",liso[[i]]$charge,"+",sep="");
# 			}		
		}
	}
	return(invisible(data.frame(peaktable,isotopes,adduct,grp,stringsAsFactors=FALSE)));
}

getPeaklist<-function(object){
	xs<-object@xcmsSet
	if (nrow(xs@groups) > 0 && length(sampnames(xs))> 1) {
 		groupmat <- groups(xs)
# 		peaklist <- data.frame(cbind(groupmat,groupval(xs, "medret", "into")),row.names = NULL)
# 		index<-which(sampclass(xs)==category)[1]+sample-1;
		peaklist<- cbind(groupmat,groupval(xs, "medret", "into"))
		cnames <- colnames(peaklist)
		if (cnames[1] == 'mzmed') cnames[1] <- 'mz' else stop ('Peak information ?!?')
		if (cnames[4] == 'rtmed') cnames[4] <- 'rt' else stop ('Peak information ?!?')
      		colnames(peaklist) <- cnames
# 		if(!is.null(object@grp_info))peaklist <- peaklist[object@grp_info,,drop=FALSE];
# 		ts <- xs@peaks[which(xs@peaks[,"sample"]==index),]
  	}else
	{
  		peaklist<-object@peaks;
	}
	adduct<-vector("character",nrow(object@peaks));
	isotopes<-vector("character",nrow(object@peaks));
	pcgroup<-vector("character",nrow(object@peaks));
	for(i in 1:length(isotopes))
	{
		if(length(object@derivativeIons)>0 && !(is.null(object@derivativeIons[[i]])))
		{
			if(length(object@derivativeIons[[i]])>1)
			{	
				names<-c();
				for(ii in 1:length(object@derivativeIons[[i]]))
				{
					names<-paste(names,object@derivativeIons[[i]][[ii]]$name,object@derivativeIons[[i]][[ii]]$mass);
				}
				adduct[i]<-names;
			}
			else
			{
				adduct[i]<-paste(object@derivativeIons[[i]][[1]]$name,object@derivativeIons[[i]][[1]]$mass);
			}
		}else adduct[i]<-"";
		if(length(object@isotopes)>0 && !is.null(object@isotopes[[i]]))
		{
				num_iso<-object@isotopes[[i]]$iso;
 				if(num_iso==0)	str_iso<-"[M]"
 				else {str_iso<-paste("[M+",num_iso,"]",sep="")}
				isotopes[i]<-paste("[",object@isotopes[[i]]$y,"]",str_iso,object@isotopes[[i]]$charge,"+",sep="");
# 			}		
		}else isotopes[i]<-"";
# 		if(!is.null(object@isotopes[[i]])) isotopes[i]<-object@isotopes[[i]];
# 		if(!is.null(object@derivativeIons[[i]])) adduct[i]<-object@derivativeIons[[i]];
	}
	for(i in 1:length(object@pspectra))
	{
		index<-object@pspectra[[i]];
		pcgroup[index]<-i;
	}
# 	pcgroup[which(pcgroup==FALSE)]<-0;
# 	pcgroup<-as.numeric(pcgroup);
	if(!is.null(object@grp_info))
	{
		isotopes_tmp<-vector("character",nrow(peaklist));adduct_tmp<-vector("character",nrow(peaklist));pcgroup_tmp<-vector("character",nrow(peaklist));
		for(i in 1:length(isotopes))
		{
			isotopes_tmp[object@grp_info[i]]<-isotopes[i];
			adduct_tmp[object@grp_info[i]]<-adduct[i];
			pcgroup_tmp[object@grp_info[i]]<-pcgroup[i];
		}
		isotopes<-isotopes_tmp;
		adduct<-adduct_tmp;
		pcgroup<-pcgroup_tmp;
# 		pcgroup[which(pcgroup==0)]<-"";
	}
	return(invisible(data.frame(peaklist,isotopes,adduct,pcgroup,stringsAsFactors=FALSE)));
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
calcFragments <- function (){
	name<-c();massdiff<-c();
	name 		<- c(
		"-CH2",
		"-CH3",
		"-NH3",
		"-H20",
		"-C2H2",
		"-CO",
		"-CO2"
		);
	massdiff	<- c(
		14.01565,
		15.02347,
		17.02655,
		18.01056,
		26.01565,
		27.99491,
		43.98982
		);
	return(data.frame(name,massdiff));
}

calcRules <- function (maxcharge=3,mol=3,nion=2,nnloss=1,nnadd=1,nh=2,polarity=NULL)
{
	name<-c();nmol<-c();charge<-c();massdiff<-c();oidscore<-c();quasi<-c();ips<-c();
	##Read Tabellen
# 	ionlist<-"/home/ckuhl/src/R/CAMERA/inst/lists/ions.csv"
	ionlist <- system.file('lists/ions.csv', package = "CAMERA")[1]
  	if (!file.exists(ionlist)) stop('ionlist.csv not found.')
	ionlist<-read.table(ionlist,header=TRUE,dec=".",sep=",",as.is=TRUE);

# 	neutralloss<-"/home/ckuhl/src/R/CAMERA/inst/lists/neutralloss.csv"
	neutralloss <- system.file('lists/neutralloss.csv', package = "CAMERA")[1]
  	if (!file.exists(neutralloss)) stop('neutralloss.csv not found.')
	neutralloss <- read.table(neutralloss,header=TRUE,dec=".",sep=",",as.is=TRUE);
# 
	neutraladdition <- system.file('lists/neutraladdition.csv', package = "CAMERA")[1]
# 	neutraladdition<-"/home/ckuhl/src/R/CAMERA/inst/lists/neutraladdition.csv"
  	if (!file.exists(neutraladdition)) stop('neutraladdition.csv not found.')
	neutraladdition <- read.table(neutraladdition,header=TRUE,dec=".",sep=",",as.is=TRUE);
	##End Read Tabellen


	##Erzeuge Regeln
	tmpname<-c();tmpnmol<-c();tmpcharge<-0;tmpmass<-0;tmpips<-0;
	#Molekülionen
	if(polarity=="positive")
	{
        	#Wasserstoff, hard codiert
		for(k in 1:mol)
		{
			if(k==1){str<-"";tmpips<-1;}else{str<-k;tmpips<-0.5};
			name<-append(name,paste("[",str,"M+H]+",sep=""));  charge<-append(charge,1);massdiff<-append(massdiff,1.0076);nmol<-append(nmol,k);if(k==1){quasi<-append(quasi,1);}else{quasi<-append(quasi,0);};oidscore<-append(oidscore,1);ips<-append(ips,tmpips)
			name<-append(name,paste("[",str,"M+2H]2+",sep=""));charge<-append(charge,2);massdiff<-append(massdiff,2.0152);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,2);ips<-append(ips,tmpips)
			name<-append(name,paste("[",str,"M+3H]3+",sep=""));charge<-append(charge,3);massdiff<-append(massdiff,3.0228);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,3);ips<-append(ips,tmpips)
			oid<-3;
			for(i in 1:nrow(ionlist))
			{
				if(ionlist[i,2]<=0)next;
				if(ionlist[i,2]==1)
				{
					name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]2+",sep=""));
				}else
				{
					name<-append(name,paste("[",str,"M+H+",ionlist[i,1],"]",ionlist[i,2]+1,"+",sep=""));
				}
				charge	<-	append(charge,ionlist[i,2]+1);
				massdiff<-	append(massdiff,ionlist[i,3]+1.0076);
				nmol	<-	append(nmol,k);
				quasi	<-	append(quasi,0);
				oidscore<-append(oidscore,oid+i);
				ips<-append(ips,tmpips);
				#Austausch
			}
			oid<-oid+nrow(ionlist);
			coeff<-expand.grid(rep(list(0:nion),nrow(ionlist)))
			if(length(list<-which(ionlist[,2]<=0))>0)
			{
				coeff[,list]<-0;
			}
			coeff<-unique(coeff);
			coeff<-cbind(coeff,rep(0,nrow(coeff)));
			coeff<-coeff[-1,]
			tmp<-NULL;
			for(i in 1:nrow(ionlist))
			{
				if(ionlist[i,2]<=0)next;
				#Austausch erstmal nur einen pro Ion
				tmp<-rbind(tmp,t(apply(coeff,1,function(x) {x[i]<-x[i]+1;x[nrow(ionlist)+1]<-1;x})));
			}
			coeff<-unique(rbind(coeff,tmp));
			for(i in 1:nrow(coeff))
			{
				tmpname<-paste("[",str,"M",sep="");
				tmpcharge<-0;tmpmass<-0;
				for(ii in 1:(ncol(coeff)-1))
				{
					if(coeff[i,ii]>0)
					{
						if(coeff[i,ii]>1)
						{
							tmpname<-paste(tmpname,"+",coeff[i,ii],ionlist[ii,1],sep="");
						}else
						{
							tmpname<-paste(tmpname,"+",ionlist[ii,1],sep="");
						}
						tmpcharge<-tmpcharge+coeff[i,ii]*ionlist[ii,2];
						tmpmass<-tmpmass+coeff[i,ii]*ionlist[ii,3]
					}
				}
				if(coeff[i,ncol(coeff)]>0)
				{#Austausch hat stattgefunden, einfach bsp 1
					tmpname<-paste(tmpname,"-H",sep="");
					tmpcharge<-tmpcharge-1;
					tmpmass<-tmpmass-1.0076;
					tmpips<-0.25;
				}
				if(tmpcharge>1)
				{
					tmpname<-paste(tmpname,"]",tmpcharge,"+",sep="")
				}else
				{
					tmpname<-paste(tmpname,"]+",sep="")
				}
				name<-append(name,tmpname)
				charge<-append(charge,tmpcharge)
				massdiff<-append(massdiff,tmpmass)
				nmol	<-append(nmol,k);
				oidscore<-append(oidscore,oid+i)
				if(sum(coeff[i,])==1&& k==1)
				{
					quasi   <-append(quasi,1);
				}else
				{
				 	quasi	<-append(quasi,0);
				}
				ips	<-append(ips,tmpips);
			}
		}
		oid<-max(oidscore);
		##Erzeuge Neutral Addition
		index<-which(quasi==1)
		for(i in 1:nrow(neutraladdition))
		{
			if(length(index2<-which(ionlist[,2]>0))>0)
			{
				for(ii in 1:length(index2))
				{
					if(ionlist[index2[ii],2] > 1)
					{
						name	<-	append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]",abs(ionlist[index2[ii],2]),"+",sep=""));
					}else
					{
						name	<-	append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]+",sep=""));
					}
					charge	<-	append(charge,ionlist[index2[ii],2]);
					massdiff<-	append(massdiff,neutraladdition[i,2]+ionlist[index2[ii],3]);
					nmol	<-	append(nmol,1);
					quasi	<-	append(quasi,0);
					oidscore<-append(oidscore,oid+1);oid<-oid+1;
					ips<-append(ips,0.5);
				}
			}
			name<-append(name,paste("[M+H+",neutraladdition[i,1],"]+",sep=""));
			charge<-append(charge,+1);
			massdiff<-	append(massdiff,neutraladdition[i,2]+1.0076);
			nmol<-append(nmol,1);
			quasi<-append(quasi,0)
			oidscore<-append(oidscore,oid+1);oid<-oid+1;
			ips<-append(ips,0.5);
		}
		ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
		if(length(index<-which(ruleset[,"charge"]>maxcharge))>0)
		{
			ruleset<- ruleset[-index,];
		}
	}else if(polarity=="negative")
	{
        	#Wasserstoff, hard codiert
		for(k in 1:mol)
		{
			if(k==1){str<-"";tmpips<-1;}else{str<-k;tmpips<-0.5};
			name<-append(name,paste("[",str,"M-H]-",sep=""));  charge<-append(charge,-1);massdiff<-append(massdiff,-1.0076);nmol<-append(nmol,k);if(k==1){quasi<-append(quasi,1);}else{quasi<-append(quasi,0);};oidscore<-append(oidscore,1);ips<-append(ips,tmpips)
			name<-append(name,paste("[",str,"M-2H]2-",sep=""));charge<-append(charge,-2);massdiff<-append(massdiff,-2.0152);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,2);ips<-append(ips,tmpips)
			name<-append(name,paste("[",str,"M-3H]3-",sep=""));charge<-append(charge,-3);massdiff<-append(massdiff,-3.0228);nmol<-append(nmol,k);quasi<-append(quasi,0);oidscore<-append(oidscore,3);ips<-append(ips,tmpips)
			oid<-3;
			for(i in 1:nrow(ionlist))
			{
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
				if(ionlist[i,2]== -1)
				{
					name<-append(name,paste("[",str,"M-H+",ionlist[i,1],"]2-",sep=""));
				}else
				{
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
			if(length(list<-which(ionlist[,2]>=0))>0)
			{
				coeff[,list]<-0;
			}
			coeff<-unique(coeff);
			coeff<-cbind(coeff,rep(0,nrow(coeff)));
			coeff<-coeff[-1,]
#  			tmp<-NULL;
#  			for(i in 1:nrow(ionlist))
#  			{
#  				if(ionlist[i,2]<=0 | k>1)next;
#  				#Austausch erstmal nur einen pro Ion
# 				tmp<-rbind(tmp,t(apply(coeff,1,function(x) {x[i]<-x[i]+1;x[nrow(ionlist)+1]<-1;x})));
#  			}
#  			coeff<-unique(rbind(coeff,tmp));
			for(i in 1:nrow(coeff))
			{
				tmpname<-paste("[",str,"M",sep="");
				tmpcharge<-0;tmpmass<-0;
				for(ii in 1:(ncol(coeff)-1))
				{
					if(coeff[i,ii]>0)
					{
						if(coeff[i,ii]>1)
						{
							tmpname<-paste(tmpname,"+",coeff[i,ii],ionlist[ii,1],sep="");
						}else
						{
							tmpname<-paste(tmpname,"+",ionlist[ii,1],sep="");
						}
						tmpcharge<-tmpcharge+coeff[i,ii]*ionlist[ii,2];
						tmpmass<-tmpmass+coeff[i,ii]*ionlist[ii,3]
					}
				}
				if(coeff[i,ncol(coeff)]>0)
				{#Austausch hat stattgefunden, einfach bsp 1
					tmpname<-paste(tmpname,"-H",sep="");
					tmpcharge<-tmpcharge-1;
					tmpmass<-tmpmass-1.0076;
					tmpips<-0.5;
				}
				if(tmpcharge< -1)
				{
					tmpname<-paste(tmpname,"]",abs(tmpcharge),"-",sep="")
				}else
				{
					tmpname<-paste(tmpname,"]-",sep="")
				}
				name<-append(name,tmpname)
				charge<-append(charge,tmpcharge)
				massdiff<-append(massdiff,tmpmass)
				nmol	<-append(nmol,k);
				oidscore<-append(oidscore,oid+i)
				if(sum(coeff[i,])==1&& k==1)
				{
					quasi   <-append(quasi,1);
				}else
				{
				 	quasi	<-append(quasi,0);
				}
				ips	<-append(ips,tmpips);
			}
		}
		oid<-max(oidscore);
		##Erzeuge Neutral Addition
		index<-which(quasi==1)
		for(i in 1:nrow(neutraladdition))
		{
			if(length(index2<-which(ionlist[,2]<0))>0)
			{
				for(ii in 1:length(index2))
				{
					if(ionlist[index2[ii],2]< -1)
					{
						name	<-	append(name,paste("[M+",ionlist[index2[ii],1],"+",neutraladdition[i,1],"]",abs(ionlist[index2[ii],2]),"-",sep=""));
					}else
					{
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
		ruleset <- data.frame(name,nmol,charge,massdiff,oidscore,quasi,ips)
		if(length(index<-which(ruleset[,"charge"]< -maxcharge))>0)
		{
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
  
  if (length(cc) > 0)
  {
   	NG <- list();NG2<-list();	
	cat('\n Calculating graph cross linking... \n % finished: '); lp <- -1;
   	for (cci in 1:length(cc))
	{
# 		cat(cci," ");
		perc <- round((cci) / length(cc) * 100)
    		if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
    		if (.Platform$OS.type == "windows") flush.console()
		## verify all correlation graphs
       		G <- subGraph(cc[[cci]],OG) 
       		if (length(nodes(G)) > 2)
		{
			## decomposition might be necessary 
# 			G <- removeSelfLoops(G)
            		hcs <-  highlyConnSG(G)
            		lsg <- sapply(hcs$clusters,function(x) length(x))
            		lsg.i <- which(lsg > 1)
            		if (length(lsg.i)<1) next; 
            		for (i in 1:length(lsg.i))
			{
                		NG[[length(NG)+1]] <- subGraph(hcs$clusters[[lsg.i[i]]],G)
            		}
       		}
		else	{  NG[[length(NG)+1]] <- G; }
   	}
	
   	## calculate all new pspectra
   	for (j in 1:length(NG))
	{ 
       		NG2[[j]] <- as.numeric(nodes(NG[[j]]))    
   	}
  }
  return(NG2)
}

getAllEICs <- function(xs,index=1) {

    peaki <- getPeaksIdxCol(xs,NULL);
   if(is.matrix(peaki))peaki<-peaki[,index]
#   nfiles <- length(xs@filepaths)
	scantimes <- list()
	maxscans <- 0
#   if (nfiles > 1) { 
#       cat('Searching maxima .. \n') 
#       for (f in 1:nfiles){
#         cat('Reading raw data file:',xs@filepaths[f]) 
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
 if (!is.double(xraw@env$mz) || !is.double(xraw@env$intensity) || !is.integer(xraw@scanindex)  ) stop('mz/int not double.')
  
 if (length(timerange) >= 2) {
      timerange <- range(timerange)
      tidx <- which((xraw@scantime >= timerange[1]) & (xraw@scantime <= timerange[2]))
      scanrange <- range(tidx)
  } else if (length(scanrange) < 2)
      scanrange <- c(1, length(xraw@scantime))
  else
      scanrange <- range(scanrange)
  
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

calcCL <-function(object,xs, EIC, scantimes, cor_eic_th){

	CL <- vector("list",nrow(object@peaks))
  	CIL <- list()
	ncl<-length(CL);npeaks=0;
	npspectra <- length(object@pspectra);
	cat('\n Calculating peak correlations... \n % finished: '); lp <- -1;

	#Wenn groupsFWHM nicht vorher aufgerufen wurde!
	if(npspectra<1)
	{
		npspectra<-1;object@pspectra[[1]]<-seq(1:nrow(object@peaks));
	}

	for(i in 1:npspectra)
	{
		pi <- object@pspectra[[i]];
		
		#percent output
# 		cat(i);
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
# 	cat("\n");
	if (length(CIL) >0) CI <- data.frame(t(sapply(CIL,function(x) x$p)),sapply(CIL,function(x) x$cor) )   
	else return(NULL)
	
	colnames(CI) <- c('xi','yi','cors')
	return(invisible(list(CL=CL,CI=CI)))
}
###End xsAnnotate intern Function###