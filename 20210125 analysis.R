library(ggplot2)
library (plyr)
library(sp)
library(rgeos)
library(raster)



#PART 0: DATA PREP 

#load files. ID connects the different datasets. 
#locations with serial number, ID, date, observer initials, perch height and diameter, X and Y coordinates, Day of observation (relative to first day of Observation which is coded as 1)
loc=read.csv('locations.csv')
#all individuals included once, with sex and SVL and date and day of first observation
ind=read.csv('allind.csv')
#offspring ID and IDs of assigned parents
off=read.csv('offspring.csv')
#display behavior data for a subset of sampled males
behav=read.csv("behavior.csv")


#removing all rows without location data, three females collected for eggs but with no coordinates (173, 174, 175).

loc=subset(loc, !is.na(loc$Xcoord))

#summarizing the number of observations per individual

obs.sum=ddply(loc, .(ID), summarize, numberobs=length(ID))
min(obs.sum$numberobs)
max(obs.sum$numberobs)
mean(obs.sum$numberobs)
sd(obs.sum$numberobs)
median(obs.sum$numberobs)


#make new data frame of all mothers (subset of all females) and all potential fathers (all sampled males)

mothers=as.vector(levels(as.factor(off$Mother)))
males=as.vector(ind$ID[ind$Sex=="m"])

#create a dataset of all possible pairs of mothers and males in the site. 
#all questions asked here are about this set of pairs--which of them produce offspring, and which do not?

pairs=data.frame(rep(0, length(mothers)*length(males)))
names(pairs)="mothers" 

#matching males to mothers
for (i in 1:length(mothers)){
  pairs$mothers[(((i-1)*length(males))+1):(i*length(males))]=rep(mothers[i], length(males))
}
pairs$males=rep(males, length(mothers))

#incorporating empirical data into pairs dataframe on which pairs produce offspring and which do not.

for (i in 1:(length(mothers)*length(males))){
  pairs$Sire[i]=
    #if(paste(pairs$mothers, pairs$males)[i] %in% paste(off$Mother, off$Sire)) {1} else {0}
    if(pairs$mothers[i] %in% off$Mother & pairs$males[i] %in% off$Sire[which(off$Mother==pairs$mothers[i])]) {1} else {0}
  pairs$Count[i]=length(which((paste(off$Mother, off$Sire))==(paste(pairs$mothers, pairs$males)[i])))
}

pairs$mothers=as.factor(pairs$mothers)
pairs$males=as.factor(pairs$males)


#males who sired offspring.
sires=count(pairs$males[pairs$Sire==1])
names(sires)=c("ID", "NumFemales")

#marking individuals with offspring data (mothers+ sires)
for (i in 1:nrow(ind)){
  ind$offspring[i]=if (ind$ID[i] %in% mothers | ind$ID[i] %in% sires$ID) {"Yes"} else {"No"}
}

#removing rows for which there are no observations.
#these three appear to be females who were caught towards the end of the season but not included in ecology/behavior
#this removes 26 offspring from consideration, and highlights that these results cannot be used to have an absolute sense of male reproductive success
pairs=pairs[pairs$mothers!="C173"&pairs$mothers!="C174"&pairs$mothers!="C175", ]

#this leaves us with a dataset of 24 mothers, 33 potential sires (27 within plot plus 6 outside), and therefore 792 pairs 
#of these, a total of 146 offspring were produced from 42 pairs, and the remaining 750 pairs did not produce offspring


#prepping the behavioral dataset and merging with the SVL dataset.

behav$elements=behav$Dewlap+behav$Pushups+behav$Headbobs
#to look at correlations between behavioral variables
pairs(behav[,5:12])

behav.ind=ddply(behav, .(ID), summarize, time=sum(Total.time), disptime=sum(Total.display), dewlap=sum(Dewlap), pushup=sum(Pushups), headbob=sum(Headbobs), elements=sum(elements))

behav.ind$disp.rate=behav.ind$disptime/behav.ind$time
behav.ind$element.rate=behav.ind$elements/behav.ind$disptime

behav.svl=merge(ind, behav.ind, by="ID", all=TRUE)
for (i in 1:nrow(behav.svl)){
  behav.svl$behav[i]=if (is.na(behav.svl$time[i])==TRUE) {"No"} else {"Yes"}
}
behav.svl$behav=as.factor(behav.svl$behav)



#want to analyse not just counts of offspring but also the proportion of each females offspring sired by a particular male
totalcounts=ddply(pairs, .(mothers), summarize, total=sum(Count))

for (i in (1:nrow(pairs))){
  pairs$prop[i]=pairs$Count[i]/totalcounts$total[totalcounts$mothers == pairs$mothers[i]]
}

#to ask if number of eggs and number of sires are correlated
#offspring numbers for males, and num sires 
off.num=ddply(off, .(Mother), summarize, eggs=length(Mother), sires= length(levels(as.factor(Sire))))


mean(off.num$sires)
min(off.num$sires)
max(off.num$sires)

sum(off.num$sires>1)/nrow(off.num)

#PART II: TERRITORIAL ANALYSES

#redoing michele's old minimum convex polygon analyses
#locations and minimum convex polygons for individuals with three or four or more observations
#will compare regressions of number of observations vs polygon areas for the two datasets

loc2=merge(loc, ind, by="ID", all.x=TRUE)

byind=ddply(loc2, .(ID), summarize, numobs=length(ID), sex=Sex[1], sn=SerialNumber.y[1], off=offspring[1], first=min(Day.x), last=max(Day.x))
byindsub3=byind[byind$numobs>2,]


#area enclosing all sampled points
#removing males selected in a 15m buffer around plot who were potential sires only

loc3=subset(loc2, loc2$InPlot=="Yes")
totalpoly=loc3[chull(loc3$Xcoord, loc3$Ycoord), ]
Polygon(totalpoly[,7:8])@area


#individual minimum convex polygons
hulls3=list()
polypoints3=list()
polys3=list()


for(i in 1:nrow(byindsub3)){
  temp=subset(loc3, loc3$ID==byindsub3$ID[i])
  
  hulls3[[i]]=chull(temp$Xcoord, temp$Ycoord)
  
  tempcoord=data.frame(X=rep(0, length(hulls3[[i]])), Y=rep(0, length(hulls3[[i]])))
  tempcoord$ID=rep(byindsub3$ID[i], length(hulls3[[i]]))
  tempcoord$Sex=rep(byindsub3$sex[i], length(hulls3[[i]]))
  tempcoord$off=rep(byindsub3$off[i], length(hulls3[[i]]))
  for (j in 1:length(hulls3[[i]])){
    tempcoord$X[j]=temp$Xcoord[hulls3[[i]][j]]
    tempcoord$Y[j]=temp$Ycoord[hulls3[[i]][j]]
  }
  polypoints3[[i]] = tempcoord
  polys3[[i]]= SpatialPolygons(list(Polygons(list(Polygon(tempcoord[,1:2])), ID=i))) 
  byindsub3$area[i]=Polygon(tempcoord[,1:2])@area
}

polysdf3=ldply(polypoints3)


#was concerned momentarily by polygon warnings, but it's the mcp that has 3 corners, not the number of observations!
#not sure why it warns us about triangles...

#checking if there's a relationship between the number of observations and MCP area, to see if some individuals with lower obs numbers should be excluded
plot(byindsub3$area~byindsub3$numobs)
abline(lm((byindsub3$area~byindsub3$numobs)))
anova(lm(sqrt(byindsub3$area)~byindsub3$numobs))

#not significant, so including all individuals with 3 or more observations
#michele's results cut off at 4, but I can't recreate her results (they include a weirdly low degrees of freedom, unsure why)



#do the MCPs intersect or not?
#this is where I expect filtering to happen to only individuals with 3 or more obs, but it hasn't happened. 
for (i in 1:nrow(pairs)){
 pairs$overlap[i]= if ((pairs$mothers[i] %in% byindsub3$ID)&(pairs$males[i] %in% byindsub3$ID)) as.numeric(gIntersects(polys3[[which(byindsub3$ID==pairs$mothers[i])]],polys3[[which(byindsub3$ID==pairs$males[i])]])) else "NA"
}  

pairs$overlap=as.numeric(pairs$overlap)


#to calculate proportion of pairs with offspring that do not overlap

pairsbyind3=pairs[pairs$mothers%in%byindsub3$ID,]
pairsbyind3$overlap.num=pairsbyind3$overlap*pairsbyind3$Count
sum(pairsbyind3$overlap.num, na.rm=TRUE)/sum(pairsbyind3$Count)

#by female
mother.summary4=ddply(pairsbyind3, .(mothers), summarize, prop=sum(overlap.num, na.rm=TRUE)/sum(Count, na.rm=TRUE))
mean(mother.summary4$prop)


#by mother, counting number of overlaps
mother.summary=ddply(pairs, .(mothers), summarize, numsires=sum(Sire), numoverlap=sum(overlap, na.rm=TRUE))
mother.summary2=merge(mother.summary, byindsub3, by.x="mothers", by.y="ID")
#WOW there's only 17 females that make it through

mean(mother.summary2$numoverlap)
min(mother.summary2$numoverlap)
max(mother.summary2$numoverlap)

mean(mother.summary2$numsires)
min(mother.summary2$numsires)
max(mother.summary2$numsires)


#finally, are number of sires and number of overlapping males correlated?



#next step: males who overlap the terr.mothers females, who do and don't sire offspring. 

pairs.overlap=subset(pairs, pairs$overlap==1)
males.overlap=subset(behav.svl, behav.svl$ID %in% pairs.overlap$males)

#only two males that overlap territories of females but do not sire offspring, so can't do comparative analyses. 

whichmales=list()
for (i in 1:nrow(mother.summary)){
  whichmales[[i]]=(pairs.overlap$males[which((pairs.overlap$mothers==mother.summary$mothers[i])&(pairs.overlap$overlap==1))])
}

anova(lm(mother.summary2$numsires~mother.summary2$numoverlap))


#PART II: TIME SCALED DISTANCE ANALYSES

#split loc by ID, arrange into list, 
locbyID=dlply(loc, .(ID), subset)

#for calculating max velocity from data; adding 1 to time difference to get, essentially, max distance traveled in a single day, thus avoiding division by zero

velocityfunc=function(x1, y1, t1, x2, y2, t2){
  return(sqrt(((x1-x2)^2)+((y1-y2)^2))/(t2-t1+1))
}

distfunc=function(x1, y1, x2, y2){
  return(sqrt(((x1-x2)^2)+((y1-y2)^2)))
}

n=length(locbyID)
dmax=NULL
deltatmax=NULL
vmax=NULL
tmin=NULL
tmax=NULL
for (i in 1:n){
  k=nrow(locbyID[[i]])
  vvec=NULL
  dvec=NULL
  tvec=NULL
  for (j in 1:(k-1)){
    dvec[j]=distfunc(locbyID[[i]]$Xcoord[j], locbyID[[i]]$Ycoord[j], locbyID[[i]]$Xcoord[j+1], locbyID[[i]]$Ycoord[j+1])
    vvec[j]=velocityfunc(locbyID[[i]]$Xcoord[j], locbyID[[i]]$Ycoord[j], locbyID[[i]]$Day[j], locbyID[[i]]$Xcoord[j+1], locbyID[[i]]$Ycoord[j+1], locbyID[[i]]$Day[j+1])
    tvec[j]=abs(locbyID[[i]]$Day[j+1]-locbyID[[i]]$Day[j])
  }
  vmax[i]=max(vvec, na.rm=TRUE)
  tmin[i]=min(tvec, na.rm=TRUE)
  tmax[i]=max(tvec, na.rm=TRUE)
  dmax[i]=if (length(dvec)==0|(length(dvec)==1&is.na(dvec[1]))) {"NA"} else {dvec[which.max(vvec)]}
  deltatmax[i]=if (length(dvec)==0|(length(dvec)==1&is.na(dvec[1]))) {"NA"} else {tvec[which.max(vvec)]}
}

v=max(vmax)
veldf=data.frame(dmax, deltatmax, vmax)

#function for time-scaled distance from Lyons et al. 2013.
#went back and forth on adding 1 to the time difference 
#but ultimately decided against it because of weird behavior in aic values around s=0 when adding the 1.

tsd=function (x1, y1, t1, x2, y2, t2, s, v){
  return(sqrt(((x1-x2)^2)+((y1-y2)^2)+(s*v*(t1-t2)^2)))
}

#merging with other datasets
pairs.svl=merge(pairs, ind[,c(2,4)], by.x='males', by.y="ID")
pairs.svl$males=as.factor(pairs.svl$males)


pairs.behav=merge(pairs, behav.svl, by.x='males', by.y="ID", all.x=TRUE)
pairs.behav$males=as.factor(pairs.svl$males)

#okay, pairs.behav has all the pairs plus behavioral and SVL data for all males



#varying s to choose based on predictive power.  
svec=0.005*c(0:200)

#first, finding the minimum tsd for all male-female pairs, i.e. how close did they get in space and time?


tsdminpairvars=matrix(nrow=nrow(pairs.behav), ncol=length(svec))
tsdwhichmin1=matrix(nrow=nrow(pairs.behav), ncol=length(svec))
tsdwhichmin2=matrix(nrow=nrow(pairs.behav), ncol=length(svec))


for (i in 1:nrow(pairs.behav)){
  for (q in (1:length(svec))){
    df1=subset(loc, loc$ID==as.character(pairs.behav$mothers[i]))
    df2=subset(loc, loc$ID==as.character(pairs.behav$males[i]))
    int=NULL
    int.t=NULL
    for (j in 1:nrow(df1)){
      vec=NULL
      
      for (k in 1:nrow(df2)){
        vec[k]=tsd(df1$Xcoord[j], df1$Ycoord[j], df1$Day[j],df2$Xcoord[k], df2$Ycoord[k], df2$Day[k], svec[q], v )  
      }
      int[j]=min(na.omit(vec))
      int.t[j]=which.min(vec)
    }
    
    tsdminpairvars[i, q]=min(na.omit(int))
    tsdwhichmin1[i, q]=which.min(int)
    tsdwhichmin2[i, q]=int.t[which.min(int)]
  }
}

pairstsd=cbind(pairs.behav, as.data.frame(tsdminpairvars))

#write.csv(pairstsd, file="pairstsd.csv")
#this was the finer scale division of the 0-1 interval; to recreate figure, can load file directly

#pairstsd=read.csv("pairstsd.csv")
#pairstsd$Sire=as.factor(pairstsd$sire)

#change columns if changing model

params=NULL
aicvec=NULL

for (i in 1:length(svec)){
  temp=pairstsd[,c("males","Sire", "SVL", "disp.rate", "element.rate", paste("V",i, sep=""))]
  names(temp)=c("males", "Sire", "SVL","disp.rate", "element.rate",  "tsd")
  mod=glm(Sire~tsd, data=temp, family=binomial)  
  params[i]=mod$coefficients[-1]
  aicvec[i]=mod$aic  
}

#calibrating by the number of variables being different for s=0 and s>0
#aicvec2=NULL
#aicvec2[1]=aicvec[1]+2
#aicvec2[2:length(svec)]=aicvec[2:length(svec)]+4


s=svec[which.min(aicvec)]

which(aicvec<((min(aicvec))+2))

svec[which(aicvec<((min(aicvec))+2))]

pardf=as.data.frame(cbind(svec, params))
names(pardf)=c("s", "tsd")



#SAVE BEST FIT MODEL COEFS

pairs.behav$tsd=tsdminpairvars[,which.min(aicvec)]

pairs.behav$minobsday1=tsdwhichmin1[,which.min(aicvec)]
pairs.behav$minobsday2=tsdwhichmin2[,which.min(aicvec)]

pairs.behav$daydiff=abs(pairs.behav$minobsday1-pairs.behav$minobsday2)


pairs.behav$tsd0=tsdminpairvars[,1]

bestfitmod=glm(Sire~tsd+SVL+disp.rate+element.rate, data=pairs.behav, family=binomial)

bestfitcoefs=bestfitmod$coefficients[-1]

#now, resampling the response variable

pairs.random=pairs.behav[,c("Sire", "tsd", "SVL", "disp.rate", "element.rate")]


nrand=10000
randcoefs=matrix(nrow=nrand, ncol=4)
p=NULL

upper=NULL
lower=NULL

set.seed(456)

for (i in 1:nrand){
      temp=pairs.random
      temp$Sire=sample(temp$Sire)
    mod.rand=glm(Sire~tsd+SVL+disp.rate+element.rate, data=temp, family=binomial)
    randcoefs[i,]=mod.rand$coefficients[-1]
  }
  

for (j in 1:4){
    p[j]=sum(as.vector(abs(randcoefs[,j]-mean(randcoefs[,j])))>=(abs(bestfitcoefs[j]-mean(randcoefs[,j]))))/nrand
    upper[j]=quantile(randcoefs[,(j)], 0.975)
    lower[j]=quantile(randcoefs[,(j)], 0.025)
}

#for one-tailed p-values
for (j in 1:4){
  p[j]=if (bestfitcoefs[j]<0) (sum(as.vector(randcoefs[,j]<bestfitcoefs[j]))/nrand) else (sum(as.vector(randcoefs[,j]>bestfitcoefs[j]))/nrand)
  upper[j]=quantile(randcoefs[,(j)], 0.975)
  lower[j]=quantile(randcoefs[,(j)], 0.025)
}

names(p)=c("tsd", "SVL", "disp.rate", "element.rate")
names(upper)=c("tsd", "SVL", "disp.rate", "element.rate")
names(lower)=c("tsd", "SVL", "disp.rate", "element.rate")


#next step: subset to only with offspring, calculate proportion of offspring of each female sired by a particular male
#repeat above analysis, find appropriate error distribution (i.e. not binomial anymore), may need a transformation
#do the same analysis again, with same predictor variables

pairs.off=subset(pairs.behav[pairs.behav$Count>0,])

eggcount=ddply(pairs, ("mothers"), summarize, totaloff=sum(Count))

for (i in 1:nrow(pairs.off)){
  pairs.off$prop.off[i]=pairs.off$Count[i]/eggcount$totaloff[which(eggcount$mothers==pairs.off$mothers[i])]
}

#proportion of offspring
bestfitmod2=glm(prop.off~tsd+SVL+disp.rate+element.rate, data=pairs.off, family="gaussian")
bestfitcoefs2=bestfitmod2$coefficients[-1]
pairs.off.random=pairs.off


randcoefs2=matrix(nrow=nrand, ncol=4)
p2=NULL
upper2=NULL
lower2=NULL

set.seed(123)

for (i in 1:nrand){
  temp2=pairs.off.random
  temp2$prop.off=sample(temp2$prop.off)
  mod.rand2=glm(prop.off~tsd+SVL+disp.rate+element.rate, data=temp2, family="gaussian")
  randcoefs2[i,]=mod.rand2$coefficients[-1]
}


for (j in 1:4){
  p2[j]=sum(as.vector(abs(randcoefs2[,j]-mean(randcoefs2[,j])))>=(abs(bestfitcoefs2[j]-mean(randcoefs2[,j]))))/nrand
  upper2[j]=quantile(randcoefs2[,(j)], 0.975)
  lower2[j]=quantile(randcoefs2[,(j)], 0.025)
}

#for one-tailed p-values
for (j in 1:4){
  p2[j]=if (bestfitcoefs2[j]<0) (sum(as.vector(randcoefs2[,j]<bestfitcoefs2[j]))/nrand) else (sum(as.vector(randcoefs2[,j]>bestfitcoefs2[j]))/nrand)
  upper[j]=quantile(randcoefs[,(j)], 0.975)
  lower[j]=quantile(randcoefs[,(j)], 0.025)
}
names(p2)=c("tsd", "SVL", "disp.rate", "element.rate")
names(upper2)=c("tsd", "SVL", "disp.rate", "element.rate")
names(lower2)=c("tsd", "SVL", "disp.rate", "element.rate")


#number of offspring
bestfitmod3=glm(Count~tsd+SVL+disp.rate+element.rate, data=pairs.off, family="poisson")
bestfitcoefs3=bestfitmod3$coefficients[-1]



randcoefs3=matrix(nrow=nrand, ncol=4)
p3=NULL
upper3=NULL
lower3=NULL

set.seed(456)

for (i in 1:nrand){
  temp3=pairs.off.random
  temp3$Count=sample(temp3$Count)
  mod.rand3=glm(Count~tsd+SVL+disp.rate+element.rate, data=temp3, family="poisson")
  randcoefs3[i,]=mod.rand3$coefficients[-1]
}


for (j in 1:4){
  p3[j]=sum(as.vector(abs(randcoefs3[,j]-mean(randcoefs3[,j])))>=(abs(bestfitcoefs3[j]-mean(randcoefs3[,j]))))/nrand
  upper3[j]=quantile(randcoefs3[,(j)], 0.975)
  lower3[j]=quantile(randcoefs3[,(j)], 0.025)
}

#for one-tailed p-values
for (j in 1:4){
  p3[j]=if (bestfitcoefs3[j]<0) (sum(as.vector(randcoefs3[,j]<bestfitcoefs3[j]))/nrand) else (sum(as.vector(randcoefs3[,j]>bestfitcoefs3[j]))/nrand)
  upper[j]=quantile(randcoefs[,(j)], 0.975)
  lower[j]=quantile(randcoefs[,(j)], 0.025)
}
names(p3)=c("tsd", "SVL", "disp.rate", "element.rate")
names(upper3)=c("tsd", "SVL", "disp.rate", "element.rate")
names(lower3)=c("tsd", "SVL", "disp.rate", "element.rate")




#PART III: FIGURES
plot1=ggplot()

#Figure 1. Sampling Durations

plot1+theme_light(15)+
  geom_rect(aes(ymin=(2*byind$sn)-0.25, ymax=(2*byind$sn)+0.25,  xmin=byind$first, xmax=byind$last, fill=byind$sex, alpha=byind$off))+
  geom_point(aes(x=loc2$Day.x, y=(2*loc2$SerialNumber.y),color=loc2$Sex, alpha=loc2$offspring))+
  xlab("Day of Sampling")+ylab("Individual")+
  scale_fill_manual(values=c("#E5BA52", "#991B37"), name="Sex", labels=c("Female", "Male"))+
  scale_color_manual(values=c("#E5BA52", "#991B37"))+guides(color=FALSE)+
  scale_alpha_discrete(range=c(0.35, 1), name="Sire/Mother?")+
  theme(axis.text.y=element_blank())


#Figure 2. Spatial locations

plot1+geom_polygon(aes(x=totalpoly$Xcoord, y=totalpoly$Ycoord), fill=NA, color="black", linetype="dashed")+
  geom_point(aes(x=loc2$Xcoord, y=loc2$Ycoord, color=loc2$Sex), size=3, alpha=0.5)+scale_color_manual(values=c("#E5BA52", "#991B37"), name="Sex", labels=c("Female", "Male"))+theme_light(15)+xlab("X-coordinate")+ylab("Y-coordinate")+
  geom_polygon(aes(x=polysdf3$X, y=polysdf3$Y, by=polysdf3$ID, color=polysdf3$Sex), size=1, fill="white", alpha=0)+
  scale_size_discrete(range=c(1,2))+guides(size=FALSE)+theme(legend.position = c(0.9, 0.15))


  #geom_polygon(aes(x=polysdf3$X[polysdf3$ID%in%"C165"],  y=polysdf3$Y[polysdf3$ID=="C165"]))+
  #geom_polygon(aes(x=polysdf3$X[polysdf3$ID%in%"C133"],  y=polysdf3$Y[polysdf3$ID=="C133"]))



#geom_polygon(aes(x=polysdf3$X[polysdf3$ID %in%males.overlap$ID],  y=polysdf3$Y[polysdf3$ID %in%males.overlap$ID]))

#so weird that it says it's ignoring the 'by' aesthetic, since it influences whether polygons get made by ID or by sex!!

#figure to illustrate 

polys133=polysdf3[polysdf3$ID%in%whichmales[[9]],]


plot1+ theme_light(15)+xlab("X-coordinate")+ylab("Y-coordinate")+
  geom_polygon(aes(x=polysdf3$X[polysdf3$ID%in%"C133"],  y=polysdf3$Y[polysdf3$ID=="C133"]), fill="#E5BA52", color="black", size=2, linetype="dashed")+
  geom_polygon(aes(x=polys133$X, y=polys133$Y, by=polys133$ID), fill="#991B37", alpha=0.25, color="black")+
  geom_polygon(aes(x=polysdf3$X[polysdf3$ID%in%"C112"],  y=polysdf3$Y[polysdf3$ID=="C112"]), fill="#991B37", alpha=0.8, color="black", size=2)+
  geom_polygon(aes(x=polysdf3$X[polysdf3$ID%in%"C165"],  y=polysdf3$Y[polysdf3$ID=="C165"]), fill="#991B37", alpha=0.8, color="black", size=2)
  
C133tsd=pairs.behav[pairs.behav$mothers=="C133",]
C133tsd[C133tsd$Sire==1,]
C133tsdoverlap=pairs.behav[pairs.behav$mothers=="C133"&pairs.behav$males%in%whichmales[[9]],]
  
plot1+theme_light()+geom_histogram(aes(x=C133tsd$tsd), color="#7D5494", fill="#7D5494", alpha=0.5)+
  geom_histogram(aes(x=C133tsdoverlap$tsd), color="#7D5494", fill="#7D5494")+
  geom_point(aes(x=C133tsd[C133tsd$Sire==1,]$tsd, y=c(1.5,4.5)), size=3)+
  xlab("Time-Scaled Distance")+ylab("Count")


pairs.behav[(pairs.behav$males%in%("C165")&pairs.behav$mothers%in%("C133")),]
pairs.behav[(pairs.behav$males%in%("C112")&pairs.behav$mothers%in%("C133")),]

library("viridisLite")

locotherm=loc2[loc2$Sex=="m"&loc2$ID!="C165"&loc2$ID!="C112",]

plot1+theme_light()+geom_point(aes(x=loc$Xcoord[loc$ID=="C133"], y=loc$Ycoord[loc$ID=="C133"], color=loc$Day[loc$ID=="C133"]), size=6, alpha=0.7,shape=15)+
  geom_point(aes(x=loc$Xcoord[loc$ID=="C165"], y=loc$Ycoord[loc$ID=="C165"], color=loc$Day[loc$ID=="C165"]), size=6, alpha=0.7,shape=17)+
  geom_point(aes(x=loc$Xcoord[loc$ID=="C112"], y=loc$Ycoord[loc$ID=="C112"], color=loc$Day[loc$ID=="C112"]), size=6,shape=18)+
  scale_color_viridis_c(name="Day Observed")+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=8, ymax=12.5), alpha=0, color="black")+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=8, ymax=12.5), alpha=0, color="black")+
  geom_rect(aes(xmin=8, xmax=9.5, ymin=-6.25, ymax=-4), alpha=0, color="black")+
  xlab("X (in m)")+ylab("Y (in m)")+
  #geom_point(aes(x=locotherm$Xcoord, y=locotherm$Ycoord))+
  annotate(geom="label", x=6.6, y = 6.6, label = "C133 and C165 \n tsd = 3.31", label.size=NA)+
annotate(geom="label", x=11, y = -2., label = "C133 and C112 \n tsd = 0.63", label.size=NA)
  
  
  

loc[loc$ID%in%c("C133", "C112", "C165"),]


#Figure 3. Distances to Centroids

dtclist=list()

for (i in 1:n){
  k=nrow(locbyID[[i]])
  dtcvec=NULL
  centx=mean(locbyID[[i]]$Xcoord, na.rm=TRUE)
  centy=mean(locbyID[[i]]$Ycoord, na.rm=TRUE)
  for (j in 1:k){
    dtcvec[j]=tsd(locbyID[[i]]$Xcoord[j],locbyID[[i]]$Ycoord[j], 0, centx, centy, 0, 0, 0)
  }
  dtclist[[i]]=dtcvec
}

loc2$dtc=unlist(dtclist)

plot2=ggplot(loc2[loc2$ID=="C133",])
plot2+theme_light(15)+ geom_point(aes(x=Day.x, y=dtc))+geom_line(aes(x=Day.x, y=dtc, group=ID), size=1.3)+xlab("Day of Sampling")+ylab("Distance to MCP Centroid (in m)")+scale_color_manual(values=c("#65BADA", "#7D5494"), labels=c("Females", "Males"))

plot2+theme_light(15)+ geom_point(aes(x=Day.x, y=dtcvec))+geom_line(aes(x=Day.x, y=dtcvec), size=1.3)+xlab("Day of Sampling")+ylab("Distance to MCP Centroid (in m)")

+scale_color_manual(values=c("#65BADA", "#7D5494"), labels=c("Females", "Males"))






#Figure 4.AIC score of Sire ~ ts + SVL for s from 0 to 1.

plot1+theme_light(15)+geom_point(aes(x=svec, y=aicvec), color="grey50", size=2, alpha=0.5)+
  geom_hline(aes(yintercept=min(aicvec)+2))+
  geom_point(aes(x=svec[which.min(aicvec)], y=aicvec[which.min(aicvec)]), size=7, color="orange", alpha=0.7)+
  xlab("s (parameter for weighting time in calculating time scaled distance)")+
  ylab("Model AIC")+
  geom_label(aes(x=svec[which.min(aicvec)], y=aicvec[which.min(aicvec)]-0.5, label=paste("s =", svec[which.min(aicvec)])))


#Figure 5.Density plots for time-scaled distance between males and females for pairs that do and don't produce offspring

ggplot()+geom_density(aes(x=pairs.behav$tsd[pairs.behav$Sire=="0"], y=..count../sum(..count..), fill="No"))+geom_density(aes(x=pairs.behav$tsd[pairs.behav$Sire=="1"], y=..count../sum(..count..), fill="Yes"), alpha=0.85)+xlab("Minimum Time Scaled Distance between males and females")+ylab("Frequency")+theme_light(15)+scale_fill_manual("Offspring \nProduced", values = c("#E7C3F3", "#7D5494"))


#Figure 6. tsd and svl for sires

ggplot()+geom_point(aes(x=pairs.behav$tsd[pairs.behav$Sire=="0"], y=pairs.behav$SVL[pairs.behav$Sire=="0"], shape="No"), size=2, alpha=0.3)+geom_point(aes(x=pairs.behav$tsd[pairs.behav$Sire=="1"], y=pairs.behav$SVL[pairs.behav$Sire=="1"], shape="Yes", size=pairs.behav$Count[pairs.behav$Sire=="1"]), alpha=0.75)+scale_size("Number Offspring", range=c(5,15))+scale_shape_manual("Offspring Produced", values = c(0, 16))+theme_light(15)+xlab("Minimum Time Scaled Distance between males and mothers")+ylab("Male Snout Vent Length (in mm)")

#a couple of other options
#lighter circles
ggplot()+geom_point(aes(x=pairs.behav$tsd[pairs.behav$Sire=="0"], y=pairs.behav$SVL[pairs.behav$Sire=="0"], shape="No"), size=2)+geom_point(aes(x=pairs.behav$tsd[pairs.behav$Sire=="1"], y=pairs.behav$SVL[pairs.behav$Sire=="1"], shape="Yes", size=pairs.behav$Count[pairs.behav$Sire=="1"]), alpha=0.2)+scale_size("Number Offspring", range=c(5,15))+scale_shape_manual("Offspring Produced", values = c(4, 19))+theme_light(15)+xlab("Minimum Time Scaled Distance between males and mothers")+ylab("Male Snout Vent Length (in mm)")




library("viridisLite")
ggplot()+theme_light(15)+geom_point(aes(x=pairs.off$SVL, y=pairs.off$prop.off, size=pairs.off$Count, color=pairs.off$mothers))+
  scale_size("Number Offspring", range=c(3,10))+scale_color_viridis_d()+
  xlab("Snout-Vent Length (in mm)")+ylab("Proportion of a Female's Offspring Sired")

#try this with tsd also, tomorrow!
        

#Figure 8. egg number vs. sire number

plot1+theme_light()+geom_point(aes(x=off.num$eggs, y=off.num$sires), color="#65BADA", size=2)+geom_smooth(aes(x=off.num$eggs, y=off.num$sires), method="lm", color="black", linetype=2)+
  xlab("Number of Eggs")+ylab("Number of Sires")

 plot1+theme_light()+geom_point(aes(x=mother.summary2$numoverlap, y=mother.summary2$numsires), color="#65BADA", size=2)+geom_smooth(aes(x=mother.summary2$numoverlap, y=mother.summary2$numsires), method="lm", color="black", linetype=2)+
   xlab("Number of Overlapping Males")+ylab("Number of Sires")  
 
 
 tiff("Fig 3.tif", height=5000, width=5000, res=600)
 