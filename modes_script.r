library(rgbif)
library(lubridate)
library(pastecs)
library(fitdistrplus)


#DOWLOAD EUROPEAN RECORDS FROM GBIF FOR EXAMPLE DATASET:
b=occ_search(taxonKey=name_backbone(name='Araschnia levana',rank='species', kingdom='animals')$usageKey,decimalLatitude='32,72',
decimalLongitude='-15,36',limit=200000)
tab=b$data

#set the smoothing parameter (higher is the sp, lower is your ability to detect close modes):
sp=1.65

#define the new column
tab$gen=NA

#add the julian day column:
tab$julian.day=yday(as.Date(paste(tab$year,tab$month,tab$day,sep="-")))
#remove non precise dated records accumulated on 1st and last day of the year:
tab=subset(tab,julian.day>1 & julian.day<365 & !is.na(decimalLatitude))

#Linear model to take in account spatial variation of the phenology+phenological shifts with time
model=lm(julian.day~decimalLatitude*decimalLongitude+as.numeric(year),data=tab) #you can add the altitude
model=step(model,scope=list(lower=~decimalLatitude*decimalLatitude,upper=model))

#stock residuals and fitted
tab$resi=residuals(model)
tab$fitt=fitted(model)
#plot the distribution of residuals
plot(density(tab$resi,adjust=sp))

#detect multimodality in residual distribution
y=density(tab$resi,kernel="gaussian",from=min(tab$resi),to=max(tab$resi),adjust=sp)$y
ae=turnpoints(y)
func=function(x){rep(ae$peaks[x],ae$exaequos[x]+1)}
func1=function(x){rep(ae$pits[x],ae$exaequos[x]+1)}
density.tab=data.frame(x=density(tab$resi,kernel="gaussian",from=min(tab$resi),to=max(tab$resi),adjust=sp)$x,y=y,
peaks=unlist(lapply(1:length(ae$peaks),func)),pits=unlist(lapply(1:length(ae$pits),func1)))
#stock local maximums:
maxi=data.frame(x=subset(density.tab,peaks==TRUE & y>=0.05*max(density.tab$y))$x,y=subset(density.tab,peaks==TRUE & y>=0.05*max(density.tab$y))$y)
#stock local minimums which are between two local maximums,removing unnecessary local minimums:
mini.func=function(j){
bide=subset(density.tab,pits==TRUE & x>maxi[j,"x"] & x<maxi[(j+1),"x"])
return(bide[which.min(bide$y),c("x","y")])}
minis=unlist(sapply(1:(nrow(maxi)-1),mini.func))
mini=data.frame(x=minis[seq(1,length(minis),2)],y=minis[seq(2,length(minis),2)])

#calculate distance between and elevation of two successive maximums and remove too close maximums
maxi$dist=0
maxi$elev=0
for(i1 in 1:(nrow(maxi)-1)){
maxi[i1,"dist"]=abs(maxi[i1,"x"]-maxi[i1+1,"x"])
maxi[i1,"elev"]=1+(min(maxi[c(i1,i1+1),"y"])-mini[[i1,"y"]])/max(maxi[,"y"])
if(maxi[i1,"dist"]<(40/maxi[i1,"elev"])){
maxi=maxi[-(i1-1+which.min(maxi[c(i1,i1+1),"y"])),]
mini=mini[-i1,]}
}

# if distribution is multimodal:
if(length(mini$x)>0){
borne1=c(min(tab$resi),mini$x)
borne2=c(mini$x,max(tab$resi))
#Fit a gaussian on each mode
for(j in 1:length(borne1)){
temp=subset(tab,resi>borne1[j] & resi<=borne2[j])
fit=fitdist(temp$resi,"norm",method="mle",start=list(mean = mean(temp$resi), sd=sd(temp$resi)))
#define the density of each residual for this mode
tab[,paste("V",j,sep="")]=dnorm(tab$resi,fit$estimate[1],fit$estimate[2])*(nrow(temp)/nrow(tab))
}s

#calculate the ratio of the density of each mode on the total density (all modes included):
tab[,paste("V",1:length(borne1),sep="")]=tab[,paste("V",1:length(borne1),sep="")]/
apply(tab[,paste("V",1:length(borne1),sep="")],1,sum)

#attribuate a random value for each record
tab$random=runif(nrow(tab),0,1)
cumul=cumsum(as.data.frame(t(tab[,paste("V",1:length(borne1),sep="")])))
gen.func=function(x){
vec=c(tab[x,"random"]-cumul[,x])
which(vec==max(vec[which(vec<0)]))
}
tab$gen=sapply(1:nrow(tab),gen.func)

#plot the result
library(ggmap)
library(ggplot2)
library(gridExtra)

France_map=get_googlemap(center="Europe", maptype = "satellite",zoom=3,scale=2) #this step can failed, repeat this row again until it works 
France <- ggmap(France_map)+
 scale_y_continuous(limits = c(32,72), expand = c(0, 0))+
       scale_x_continuous(limits = c(-15,36), expand = c(0, 0))
pl1=ggplot(data=tab,aes(x=julian.day,fill=as.factor(gen)))+geom_histogram()+theme_bw()+theme(legend.position="none")+
ggtitle(unique(tab$name))+xlab("Observation day (julian day)")
pl2=France+geom_point(data=tab,aes(y=decimalLatitude,x=decimalLongitude,col=as.factor(gen)),size=0.4,alpha=0.4)+labs(col="Mode")
grid.arrange(pl1,pl2,ncol=2)