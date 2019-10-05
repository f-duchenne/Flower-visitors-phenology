library(lubridate)
library(data.table)
library(dispmod)

donnf=as.data.frame(fread("data_set_to_download.txt",sep="\t",header=T)) #see the onlie version of the article to get the link to download the dataser

donnf$Year=as.numeric(as.character(donnf$Year)) #convert year to a numeric variable if it is not

liste=unique(donnf$Species_gen) # define the speciesand mode liste

res=as.data.frame(matrix(NA,1,7))
names(res)=c("Species_gen","Year_effect","Year_error","Year_pval","Year_var_effect","Year_var_error","Year_var_pval")

for (i in 1:length(liste)){
tab=subset(donnf,Species_gen==liste[i,1])

#centered spatial variables to make interaction terms more readable 
tab$Latitude=tab$Latitude-mean(tab$Latitude)
tab$Longitude=tab$Longitude-mean(tab$Longitude)
tab$Altitude=tab$Altitude-mean(tab$Altitude)

#AIC Selection on the model for MFD, for some species, omptimize it manually by removing interactions:
if(liste[i,1]=="Endochironomus albipennis_1" | liste[i,1]=="Endochironomus albipennis_2" | liste[i,1]=="Psodos quadrifaria_NA" | liste[i,1]=="Rhyacia helvetina_NA" |
liste[i,1]=="Bubas bison_NA" | liste[i,1]=="Colias tyche_NA" | liste[i,1]=="Rheumaptera cervinalis_NA" | liste[i,1]=="Eudonia angustea_2" |
liste[i,1]=="Onthophagus emarginatus_NA" | liste[i,1]=="Tanypus kraatzi_NA" | liste[i,1]=="Lithophane lamda_NA" | liste[i,1]=="Euonthophagus amyntas_NA"){
model=lm(Jday~Latitude*Longitude+Altitude+Year,data=tab)
model=step(model,scope=list(lower=~Latitude+Longitude+Year,upper=model))
}else{
model=lm(Jday~Latitude*Longitude+I(Latitude^2)*I(Longitude^2)+
I(Latitude^3)*I(Longitude^3)+Altitude+Year*(Latitude+Longitude),data=tab)
model=step(model,scope=list(lower=~Latitude+Longitude+Year,upper=model))}

#Optimization of the model on the variance, catching error in order to avoid to break the loop. If there is an error,
#we removed terms in order to get a model that converges without order
model2="error"
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Year+Latitude+Altitude+Longitude,data=tab)},
error = function(e) {"error"})
manuel="non"
if(class(model2)!="dispmod"){
manuel="oui"
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Year+Latitude+Longitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Latitude+Longitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Year+Latitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Year+Longitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Latitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Longitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Altitude,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~Year,data=tab)},
error = function(e) {"error"})}
if(class(model2)!="dispmod"){
model2=tryCatch({lm.disp(model$call$formula,var.formula=~1,data=tab)},
error = function(e) {"error"})}

#Optimize the model on the variance based on the AIC
a=model2$var$aic-1
while(model2$var$aic>a){
vec=rownames(summary(model2$var)$coef)[-c(1+which.max(summary(model2$var)$coef[-1,4]))][-1]
if(length(vec)==0){vec=1}
modelb=tryCatch({lm.disp(model$call$formula,var.formula=as.formula(paste("~",paste(vec,collapse="+"),sep="")),data=tab,maxit=45)},
error = function(e) {"error"})
if(class(modelb)=="dispmod"){if(modelb$var$aic<model2$var$aic){model2=modelb}}
a=model2$var$aic
}

model2$mean$terms=model$terms
#Keep coefficients in a dataframe defined above
res[i,1]=as.character(liste[i,1])
res[i,2]=model2$mean$coef["Year"]
res[i,3]=summary(model2)$mean$coef["Year",2]
res[i,4]=summary(model2)$mean$coef["Year",4]
res[i,5]=model2$var$coef["Year"]
res[i,6]=if("Year" %in% rownames(summa2$coefficients)){summary(model2)$var$coef["Year",2]}else{NA}
res[i,7]=if("Year" %in% rownames(summa2$coefficients)){summary(model2)$var$coef["Year",4]}else{NA}
}