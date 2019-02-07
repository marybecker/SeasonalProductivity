setwd("/home/mkozlak/Projects/GitHub/SeasonalProductivity")

library(vegan)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)
library(lubridate)

SPP<- read.csv("data/SPP.csv",header=TRUE,row.names=1)
SPP[is.na(SPP)] <- 0
sites<- read.csv("data/sites.csv",header=TRUE)
sites$Collect_Date<-paste0(substr(sites$Collect_Date,6,9),"-0",
                          substr(sites$Collect_Date,1,1),"-",
                          substr(sites$Collect_Date,3,4))
#sites$Collect_Date<-as_date(sites$Collect_Date)
sites$month<-as.numeric(substr(sites$SID,4,4))
sites$Collect_Date<-ymd(sites$Collect_Date)
sites$Stream<- factor(sites$Stream,
                        levels=c("Salmon River","Norwalk River","Pequabuck River"))


##########Diversity index###########################
#####################################################
div<-diversity(decostand(SPP,"total"),index='shannon')
div<-as.data.frame(div)
div$SID<-row.names(div)
div<-merge(sites,div,by="SID")

ggplot(div,aes(x=Collect_Date,y=div,colour=Stream,group=Stream))+
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  labs(colour="Stream",x=NULL,y=NULL,title="Shannon Diversity Index")+
  theme(legend.position="bottom")

############Mantel Tests By Stream#######################
#########################################################
Stream<-"Salmon River"
distStream<-sites[sites$Stream==Stream,]
distSite<-distStream[c("Temp")]
rownames(distSite)<-distStream$SID
Site.dist<-vegdist(distSite,"euclidean")


SPPStream<-SPP[13:18,] #13:18 - Salmon, 7:12 - Pequabuck, 1:6 - Norwalk
SPPStream<- decostand(SPPStream,"hellinger")#Sqrt of rel abundance
SSPP.dist <- vegdist(SPPStream,"bray")

SPPTime<-mantel(SSPP.dist,Site.dist)
SPPTime


######################################################################
#######Similarity Plots###############################################
#####################################################################

# SPP5<- decostand(SPP,"pa")
# SPP5<- SPP[,colSums(SPP5)<5]
# n<- colnames(SPP5)
# SPP5<- SPP[,n]
# Other<-rowSums(SPP5)
# Oname<-names(Other)
# Other<- melt(as.data.frame(Other))
# row.names(Other)<-Oname
# Other$variable<- NULL
# colnames(Other)<-"Other"

#SPP<-sqrt(SPP)#Transformation when rel abund values
SPP<- decostand(SPP,"total")#rel abundance

SPP.dist <- vegdist(SPP,"bray")
distSite<-sites[c("TP")]
rownames(distSite)<-sites$SID
Site.dist<-vegdist(distSite,"euclidean")
mantel(SPP.dist,Site.dist)


SPP.dist<-as.matrix(SPP.dist)
Site.dist<-as.matrix(Site.dist)
levelplot(SPP.dist,at=seq(0,1,0.01),
          col.regions=topo.colors(100),scales=list(cex=0.4),
          xlab="",ylab="",main="Percent Difference Coefficient (Bray Curtis)")

SPDist<- melt(as.matrix(SPP.dist),varnames=c("SID","col"))
STDist<-melt(as.matrix(Site.dist),varnames=c("SID","col"))
DistCmb<-merge(SPDist,STDist,by=c("SID","col"))
DistCmb<-DistCmb[DistCmb$value.x!=0,]
colnames(DistCmb)[3:4]<-c("C","TP")
DistCmb$SSite<-ifelse(substr(DistCmb$SID, 1, 3)==substr(DistCmb$col, 1, 3),
                      substr(DistCmb$SID, 1, 3),"DIFF")

ggplot(DistCmb,aes(x=TP,y=((1-C)*100),colour=SSite))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="TP (mg/L) distance")

ggplot(DistCmb[DistCmb$SSite=="DIFF",],aes(x=TP,y=((1-C)*100),colour=SSite))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,50)+
  labs(y ="Species similarity (%)", x="TP (mg/L) distance")

ggplot(DistCmb[DistCmb$SSite!="DIFF",],aes(x=TP,y=((1-C)*100),colour=SSite))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="TP (mg/L) distance")

########Chlor a VS Other Parameters Across All Sample######
distChla<-sites[c("Chla")]
rownames(distChla)<-sites$SID
Chla.dist<-vegdist(distChla,"euclidean")
distSite<-sites[c("TP")]
rownames(distSite)<-sites$SID
Site.dist<-vegdist(distSite,"euclidean")
mantel(Chla.dist,Site.dist)


Chla.dist<-as.matrix(Chla.dist)
Site.dist<-as.matrix(Site.dist)
levelplot(SPP.dist,at=seq(0,1,0.01),
          col.regions=topo.colors(100),scales=list(cex=0.4),
          xlab="",ylab="",main="Percent Difference Coefficient (Bray Curtis)")

ChlaDist<- melt(as.matrix(Chla.dist),varnames=c("SID","col"))
STDist<-melt(as.matrix(Site.dist),varnames=c("SID","col"))
DistCmb<-merge(ChlaDist,STDist,by=c("SID","col"))
DistCmb<-DistCmb[DistCmb$value.x!=0,]
colnames(DistCmb)[3:4]<-c("C","Var")
DistCmb$SSite<-ifelse(substr(DistCmb$SID, 1, 3)==substr(DistCmb$col, 1, 3),
                      substr(DistCmb$SID, 1, 3),"DIFF")

ggplot(DistCmb,aes(x=Var,y=C,colour=SSite))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  labs(y ="Chla distance", x="TP (mg/L) distance")
#annotate("text", x=c(0.001,0.001), y=c(90,85), 
#label=c("mantel r = =0.1944","p = 0.7639"),hjust=0)


####Similarity between sites collected during the same month###########
#######################################################################
SPDist<- melt(as.matrix(SPP.dist),varnames=c("SID","col"))
write.csv(SPDist,"SPBCDist.csv")
SPDist<-merge(sites,SPDist,by="SID")
SPDist$CSID<-substr(SPDist$col,1,3)
SPDist$RSID<-substr(SPDist$SID,1,3)
SPDist$MS1<-as.numeric(substr(SPDist$SID,4,4))
SPDist$MS2<-as.numeric(substr(SPDist$col,4,4))
SPDist<-SPDist[which(SPDist$MS1-SPDist$MS2==0 & SPDist$CSID!=SPDist$RSID),]

site.comb<-SPDist[,10:11]
site.comb<-unique(site.comb[c("RSID","CSID")])
site.comb$mean<- 0

#Calcuated average similarity measure for each combination#
for (i in 1:dim(site.comb)[1]) {
s<- site.comb[i,]
c<-SPDist[which(SPDist$CSID==s$CSID[1] & SPDist$RSID==s$RSID[1]),]
c<-c[which(c$value!=0),]
m<-mean(c$value)
site.comb$mean[i]<-m
}

BTWSiteDist<- SPDist[which(SPDist$RSID=='NOR'&SPDist$CSID=='PEQ'),]
ggplot(BTWSiteDist,aes(x=month,y=((1-value)*100)))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="Month",
       title="Norwalk vs. Pequabuck")

####Similarity within sites collected over the POR###########
#######################################################################
SPDist<- melt(as.matrix(SPP.dist),varnames=c("SID","col"))
SPDist<-merge(sites,SPDist,by="SID")
SPDist$CSID<-substr(SPDist$col,1,3)
SPDist$RSID<-substr(SPDist$SID,1,3)
SPDist$MS1<-as.numeric(substr(SPDist$SID,4,4))
SPDist$MS2<-as.numeric(substr(SPDist$col,4,4))
SPDist<-SPDist[which(SPDist$MS1-SPDist$MS2!=0 & SPDist$CSID==SPDist$RSID),]
SPDist$MDiff<-abs(SPDist$MS1-SPDist$MS2)
SPDist<-merge(SPDist,sites,by.x=c("col","Station.ID","Stream"),
              by.y=c("SID","Station.ID","Stream"))
colnames(SPDist)[5]<-"SimDate"
colnames(SPDist)[23]<-"Collect_Date"
SPDist$SimDate<-as_date(SPDist$SimDate)
SPDist$Collect_Date<-as_date(SPDist$Collect_Date)
SPDist$DDiff<-abs(SPDist$Collect_Date-SPDist$SimDate)
SPDist$TPDiff<-abs(SPDist$TP.x-SPDist$TP.y)
SPDist$TempDiff<-abs(SPDist$Temp.x-SPDist$Temp.y)
SPDist$ChlaDiff<-abs(SPDist$Chla.x-SPDist$Chla.y)

# Spp Similarity Vs Time
Stream<-"Pequabuck River"
SPDistSite<-SPDist[SPDist$Stream==Stream,]
ggplot(SPDistSite,aes(x=DDiff,y=((1-value)*100)))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="Temporal distance (days)",
       title=Stream)
  #annotate("text", x=c(40,40), y=c(90,85), 
           #label=c("mantel r = 0.9653","p = 0.003"),hjust=0)

ggsave(paste(Stream,"SoftAlg.tiff"),width=4,height=4)

# Spp Similarity vs TP Diff
Stream<-"Salmon River"
SPDistSite<-SPDist[SPDist$Stream==Stream,]
ggplot(SPDistSite,aes(x=TPDiff,y=((1-value)*100)))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="TP (mg/L) distance",
       title=Stream)
  #annotate("text", x=c(0.001,0.001), y=c(90,85), 
           #label=c("mantel r = =0.1944","p = 0.7639"),hjust=0)

ggsave(paste(Stream,"SoftAlgSPPTPDiff.jpg"),width=4,height=4)

# Spp Similarity vs Chla Diff
Stream<-"Norwalk River"
SPDistSite<-SPDist[SPDist$Stream==Stream,]
ggplot(SPDistSite,aes(x=ChlaDiff,y=((1-value)*100)))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE,colour="black")+
  ylim(0,100)+
  labs(y ="Species similarity (%)", x="Chl a (mg/m2) distance",
       title=Stream)
#annotate("text", x=c(0.001,0.001), y=c(90,85), 
#label=c("mantel r = =0.1944","p = 0.7639"),hjust=0)

ggsave(paste(Stream,"SoftAlgSPPChlaDiff.jpg"),width=4,height=4)
  





SPDistM<-SPDist[SPDist$MS1==6,]
ggplot(SPDistM,aes(x=MS2,y=((1-value)*100),group=Stream))+
  geom_point(aes(colour=Stream))+
  ylim(0,100)+
  labs(y ="Species Similarity (%)", x="Month")

#######
SPDistM<-SPDist[SPDist$MS1==6,]
ggplot(SPDistM,aes(x=MS2,y=((1-value)*100),group=Stream))+
  geom_point(aes(colour=Stream))+
  ylim(0,100)+
  labs(y ="Species Similarity (%)", x="Month")


#Calcuated average similarity measure for each combination#
site.comb<-SPDist[,6:7]
site.comb<-unique(site.comb[c("RSID","CSID")])
site.comb$mean<- 0

for (i in 1:dim(site.comb)[1]) {
  s<- site.comb[i,]
  c<-SPDist[which(SPDist$CSID==s$CSID[1] & SPDist$RSID==s$RSID[1]),]
  c<-c[which(c$value!=0),]
  m<-mean(c$value)
  site.comb$mean[i]<-m
}

#####Similarity within sites seqential months##############
###########################################################

m<-6
sim <-as.data.frame(matrix(0.0,nrow=15,ncol=1))
mth <-as.data.frame(matrix(0,nrow=15,ncol=1))
site<-as.data.frame(matrix('',nrow=15,ncol=1),stringsAsFactors = F)
SPDist2<-cbind(site,mth,sim)
colnames(SPDist2)<-c('SID_1','SID_2','sim')

k<-1#assuming the matrix is sorted
for(i in 0:2){#number of sites 
  for(j in (i*m+1):(i*m+m-1)){ #one over diagonals
    SPDist2[k,1]<-rownames(SPP.dist)[j]
    SPDist2[k,2]<-colnames(SPP.dist)[j+1]
    SPDist2[k,3]<-SPP.dist[j,j+1]
    k<-k+1
  }
}

SPDist2$month<-substr(SPDist2$SID_1,4,4)
SPDist2$site<-substr(SPDist2$SID_1,1,3)

ggplot(SPDist2,aes(x=month,y=sim,group=site))+
  geom_line(aes(colour=site),linetype=2)+
  ylim(0,1)