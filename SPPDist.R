setwd("P:/Projects/GitHub_Prj/SeasonalProductivity")

library(vegan)
library(reshape2)
library(dplyr)
library(lattice)

SPP<- read.csv("data/SPP.csv",header=TRUE,row.names=1)
SPP[is.na(SPP)] <- 0
sites<- read.csv("data/sites.csv",header=TRUE)

SPP5<- decostand(SPP,"pa")
SPP5<- SPP[,colSums(SPP5)<5]
n<- colnames(SPP5)
SPP5<- SPP[,n]
Other<-rowSums(SPP5)
Oname<-names(Other)
Other<- melt(as.data.frame(Other))
row.names(Other)<-Oname
Other$variable<- NULL
colnames(Other)<-"Other"

#SPP<-sqrt(SPP)#Transformation when rel abund values
SPP<- decostand(SPP,"hellinger")#Sqrt of rel abundance

SPP.dist <- vegdist(SPP,"bray")
SPP.dist<-as.matrix(SPP.dist)
levelplot(SPP.dist,at=seq(0,1,0.01),
          col.regions=topo.colors(100),scales=list(cex=0.4),
          xlab="",ylab="",main="Percent Difference Coefficient (Bray Curtis)")


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

site.comb<-SPDist[,6:7]
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

####Similarity within sites collected over the POR###########
#######################################################################
SPDist<- melt(as.matrix(SPP.dist),varnames=c("SID","col"))
SPDist<-merge(sites,SPDist,by="SID")
SPDist$CSID<-substr(SPDist$col,1,3)
SPDist$RSID<-substr(SPDist$SID,1,3)
SPDist$MS1<-as.numeric(substr(SPDist$SID,4,4))
SPDist$MS2<-as.numeric(substr(SPDist$col,4,4))
SPDist<-SPDist[which(SPDist$MS1-SPDist$MS2!=0 & SPDist$CSID==SPDist$RSID),]

site.comb<-SPDist[,6:7]
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