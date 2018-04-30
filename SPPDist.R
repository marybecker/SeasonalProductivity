setwd("/Users/tbecker/Documents/Projects/GitHubProjects/SeasonalProductivity")

library(vegan)
library(reshape2)

SPP<- read.csv("data/SPP.csv",header=TRUE,row.names=1)
SPP[is.na(SPP)] <- 0
sites<- read.csv("data/sites.csv",header=TRUE)

SPP<- decostand(SPP,"total")

SPP.dist <- vegdist(SPP,"bray")

SPDist<- melt(as.matrix(SPP.dist),varnames=c("SID","col"))
write.csv(SPDist,"SPBCDist.csv")
SPDist<-merge(sites,SPDist,by="SID")
SPDist$CSID<-substr(SPDist$col,1,3)
SPDist$RSID<-substr(SPDist$SID,1,3)

unique(site.name[c("STA_SEQ","Station_Name")])
site.comb<-SPDist[,6:7]
site.comb<-unique(site.comb[c("RSID","CSID")])
site.comb$mean<- 0

##Calcuated average similarity measure for each combination##
for (i in 1:dim(site.comb)[1]) {
s<- site.comb[i,]
c<-SPDist[which(SPDist$CSID==s$CSID[1] & SPDist$RSID==s$RSID[1]),]
c<-c[which(c$value!=0),]
m<-mean(c$value)
site.comb$mean[i]<-m
}