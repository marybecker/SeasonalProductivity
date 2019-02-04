setwd("/home/mbecker/Documents/GitHub/SeasonalProductivity")

Alg<-read.csv("data/AlgBioData_121818.csv",header=TRUE)

D<-subset(Alg,Alg$Group=='N',select=c(SiteName,CollectionMonth,
                                      BioDataTaxonName,
                                      AdjLabCount))
names(D)<-c("SName","Month","Taxa","N")
D$SID<-paste0(substr(D$SName,1,3),D$Month)
D<-D[,3:5]

D<-reshape(D,idvar="SID",timevar="Taxa",
           direction="wide")
D[is.na(D)] <- 0

rownames(D) <- D[,1]
D[,1] <- NULL

colSums(D)
rowSums(D)

write.csv(D,"SoftAlgSPP.csv")