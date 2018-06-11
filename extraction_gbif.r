library(rgbif)
library(lubridate)
library(stringi)

#SET THE SPECIES LIST TO DOWNLOAD:
liste=data.frame(species=c("Thecla betulae","Episyrphus balteatus"))
#be sure that you have the right GBIF ID:
liste$ID=NA
liste$canon=NA
liste$confi=NA
liste$rang=NA
for(i in 1:nrow(liste)){
obj=name_backbone(name=as.character(liste[i,"species"]),rank='species', kingdom='animals')
if(length(obj$usageKey)>0){
liste$ID[i]=obj$usageKey
liste$family[i]=obj$family
liste$order[i]=obj$order
liste$canon[i]=obj$canonicalName
liste$rang[i]=obj$rank
liste$confi[i]=obj$confidence
}
}
liste
#define a list of countries that we do not need (avoid to keep too many unnecessary records):
Pays_nondesire=c("Thailand","Mexico","Panama","Svalbard and Jan Mayen","Malaysia","Korea, Democratic People's Republic of","Taiwan, Province of China","Pakistan","Indonesia","Swaziland",
"American Samoa","Syrian Arab Republic","Tajikistan","Uzbekistan","Kyrgyzstan","United Arab Emirates","Grenada","Bolivia, Plurinational State of","Uganda",
"Tanzania, United Republic of","Malawi","Puerto Rico","Taiwan","Kenya","Trinidad and Tobago","Afghanistan","Armenia","Iran, Islamic Republic of","Panama","Botswana",
"New Zealand","Japan","Australia","United States","Canada","Brazil","Argentina","Madagascar","Suriname","China","Nigeria","Albania","Turkey","Russian Federation","Colombia","Korea, Republic of","Mongolia","Morocco","Israel","Venezuela, Bolivarian Republic of",
"Libya","Algeria","Benin","Jordan","Greenland","French Guiana","Ghana","Chile","Guyana","Togo","Lebanon","Tunisia","Iraq","India","Kyrgyzstan","Taiwan",
"New Zealand","Japan","Australia","United States","Canada","Brazil","Argentina","Madagascar","Suriname","China","Nigeria","Turkey","Russian Federation","Colombia",
"Korea, Republic of","Mongolia","Morocco","Israel","Venezuela, Bolivarian Republic of","Libya","Algeria","Benin","Jordan","Greenland","French Guiana","Ghana","Chile","Guyana",
"Togo","Lebanon","Tunisia","Iraq")

#DATA EXTRACTION:

for(i in 1:nrow(liste)){

#CHECK that the number of records is under 200 000, if it is not hte case, go further:
if(occ_count(taxonKey=liste[i,"ID"])<200000){
#Download records
b=occ_search(taxonKey=liste[i,"ID"],limit=200000)
#check that the tabe is not empty and there is a column day:
if(length(b$data)>0 && ("day" %in% names(b$data))){
tab=b$data
#keep only precise records
tab=subset(tab,!is.na(year) & !is.na(month) & !is.na(day))
#remove larvae records
if("lifeStage" %in% names(tab)){tab$lifeStage[which(is.na(tab$lifeStage))]="unkown"
tab=subset(tab,lifeStage!="LARVA" & lifeStage!="EMRYO")}
#keep only records from europe (approximately) or not georeferenced
if(!("decimalLatitude" %in% names(tab))){tab$decimalLatitude=NA
tab$decimalLongitude=NA}
tab$decimalLatitude[which(tab$decimalLatitude==0)]=NA
tab=subset(tab,decimalLatitude>34 & decimalLatitude<72 | is.na(decimalLatitude))
tab$decimalLongitude[which(is.na(tab$decimalLatitude))]=NA
tab=subset(tab,decimalLongitude>-15 & decimalLongitude<32 | is.na(decimalLongitude))
#if the table is not empty yet after this selection:
if(dim(tab)[1]>0){
#creates columns if they do not exist (download tables are moving (not always including same columns)
if(!("identifiedBy" %in% names(tab))){tab$identifiedBy=NA}
if(!("institutionCode" %in% names(tab))){tab$institutionCode=NA}
if(!("county" %in% names(tab))){tab$county=NA}
if(!("country" %in% names(tab))){tab$country=NA}
if(!("municipality" %in% names(tab))){tab$municipality=NA}
if(!("locality" %in% names(tab))){tab$locality=NA}
if(!("sex" %in% names(tab))){tab$sex=NA}
if(!("lifeStage" %in% names(tab))){tab$lifeStage=NA}
#Organize columns
tab$species=paste(liste[i,"species"])
tab$key=paste(liste[i,"ID"])
tab=tab[,c("gbifID","species","day","month","year","country","county","municipality","locality","decimalLatitude","decimalLongitude","identifiedBy","institutionCode","key","sex","lifeStage")]
names(tab)=c("Reference","species","day","month","year","country","county","municipality","locality","Latitude","Longitude","collector","Source","ID","sex","lifeStage")
#add taxonomic informations	
tab$order=paste(liste[i,"order"])
tab$family=paste(liste[i,"family"])
#add julian day
tab$julian.day=yday(as.Date(paste(tab$year,tab$month,tab$day,sep="-")))
dattot=tab}}
}else{
#Species with more than 200 000
#DOWNLOAD DATA YEAR BY YEAR on the studied period
for(j in seq(1960,2016,1)){
b=occ_search(taxonKey=liste[i,"ID"],year=j,limit=200000)
#check that the tabe is not empty and there is a column day:
if(length(b$data)>0 && ("day" %in% names(b$data))){
tab=b$data
#keep only precise records
tab=subset(tab,!is.na(year) & !is.na(month) & !is.na(day))
#remove larvae records
if("lifeStage" %in% names(tab)){tab$lifeStage[which(is.na(tab$lifeStage))]="unkown"
tab=subset(tab,lifeStage!="LARVA" & lifeStage!="EMRYO")}
#keep only records from europe (approximately) or not georeferenced
if(!("decimalLatitude" %in% names(tab))){tab$decimalLatitude=NA
tab$decimalLongitude=NA}
tab$decimalLatitude[which(tab$decimalLatitude==0)]=NA
tab=subset(tab,decimalLatitude>34 & decimalLatitude<72 | is.na(decimalLatitude))
tab$decimalLongitude[which(is.na(tab$decimalLatitude))]=NA
tab=subset(tab,decimalLongitude>-15 & decimalLongitude<32 | is.na(decimalLongitude))
#if the table is not empty yet after this selection:
if(dim(tab)[1]>0){
#creates columns if they do not exist (download tables are moving (not always including same columns)
if(!("identifiedBy" %in% names(tab))){tab$identifiedBy=NA}
if(!("institutionCode" %in% names(tab))){tab$institutionCode=NA}
if(!("county" %in% names(tab))){tab$county=NA}
if(!("country" %in% names(tab))){tab$country=NA}
if(!("municipality" %in% names(tab))){tab$municipality=NA}
if(!("locality" %in% names(tab))){tab$locality=NA}
if(!("sex" %in% names(tab))){tab$sex=NA}
if(!("lifeStage" %in% names(tab))){tab$lifeStage=NA}
#Organize columns
tab$species=paste(liste[i,"species"])
tab$key=paste(liste[i,"ID"])
tab=tab[,c("gbif","species","day","month","year","country","county","municipality","locality","decimalLatitude","decimalLongitude","identifiedBy","institutionCode","key","sex","lifeStage")]
names(tab)=c("Reference","species","day","month","year","country","county","municipality","locality","Latitude","Longitude","collector","Source","ID","sex","lifeStage")
#add taxonomic informations	
tab$order=paste(liste[i,"order"])
tab$family=paste(liste[i,"family"])
#add julian day
tab$julian.day=yday(as.Date(paste(tab$year,tab$month,tab$day,sep="-")))
if(j==1){dattot=tab}else{dattot=rbind(dattot,tab)}
}}}}

#UNIFORMIZATION OF THE DATASET
dattot$country=as.character(dattot$country)
dattot$municipality=as.character(dattot$municipality)
dattot$county=as.character(dattot$county)
dattot$country[which(is.na(dattot$country) | dattot$country=="")]="unkown"
dattot$municipality[which(is.na(dattot$municipality) | dattot$municipality=="")]="unkown"
dattot$county[which(is.na(dattot$county) | dattot$county=="")]="unkown"
dattot$day=as.numeric(as.character(dattot$day))
dattot$month=as.numeric(as.character(dattot$month))
dattot$year=as.numeric(as.character(dattot$year))
dattot$municipality=gsub("[:punct:]","_",stri_trans_general(dattot$municipality, "latin-ascii"),fixed=T)
dattot$Source=gsub("[:punct:]","_",stri_trans_general(dattot$Source, "latin-ascii"),fixed=T)
dattot$collector=gsub("[:punct:]","_",stri_trans_general(dattot$collector, "latin-ascii"),fixed=T)
dattot$locality=gsub("[:punct:]","_",stri_trans_general(dattot$locality, "latin-ascii"),fixed=T)
dattot$county=gsub("[:punct:]","_",stri_trans_general(dattot$county, "latin-ascii"),fixed=T)
#GEOGRAPHIC SELECTION:
dattot=subset(dattot,municipality!="unkown" | county!="unkown" | country=="Denmark" | country=="Liechtenstein" | country=="Belgium" | country=="Netherlands" | country=="Luxembourg" | country=="Switzerland" |
country=="Andorra" | country=="Kosovo" | !is.na(Latitude))
dattot=dattot[which(!(dattot$country %in% Pays_nondesire)),]

#save extraction
if(i==1){dataf=dattot}else{dataf=rbind(dataf,dattot)}
#DISPLAY STATE:
print(c(i,as.character(liste[i,1]),nrow(dataf)))
}
