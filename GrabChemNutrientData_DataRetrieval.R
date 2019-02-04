setwd("/home/mbecker/Documents/GitHub/SeasonalProductivity")

samples<-read.csv("data/sites.csv",header=TRUE)
samples$Collect_Date<-mdy(samples$Collect_Date)
samples$Stream<- factor(samples$Stream,
                        levels=c("Salmon River","Norwalk River","Pequabuck River"))

library(dataRetrieval)
library(plyr)
library(ggplot2)
library(grid)
library(reshape2)
library(lubridate)


SID<- c("01189000","01209700","01193500")
Cd<-c("70957")##Chlor a
#Cd<-c("00605","00608","00613","00618","00631","00660","00665","00671","62855")
CdFlow<-"00060"
sdate<-"2016-04-01"
edate<-"2016-10-01"


data <- readNWISqw(siteNumbers = SID,
                     parameterCd = Cd,
                     startDate = sdate,
                     endDate = edate)

dsum<- data%>%
  group_by(site_no,parm_cd)%>%
  summarize(count=n(),avg=mean(result_va),
            min=min(result_va),max=max(result_va),
            med=median(result_va))

flowdata<-readNWISdata(siteNumbers = SID,
                       parameterCd = CdFlow,
                       startDate = sdate,
                       endDate = edate)


########2 Week Median######################

samples$endFlowDate<-samples$Collect_Date-14
samples$Flow<-NA

for (i in 1:dim(samples)[1]){
s<-samples[i,]
flowdata1<-flowdata[which(flowdata$site_no==paste0(0,s$Station.ID) & 
                            flowdata$dateTime>s$endFlowDate &
                            flowdata$dateTime<=s$Collect_Date),]
samples[i,12]<-median(flowdata1[,4])

}


#########Flow Plot###########################

flowday<-samples
flowday$site_no<-paste0(0,flowday$Station.ID)
flowday$dateTime<-flowday$Collect_Date
flowday$dateTime<-mdy(flowday$dateTime)
flowday<-merge(flowdata,flowday,by=c("site_no","dateTime"))
flowday$dateTime<-ymd(flowday$dateTime)

siteno<-'01209700'
flowdataSite<-flowdata[which(flowdata$site_no==siteno),]
flowdaysite<-flowday[which(flowday$site_no==siteno),]
flowdataSite$sampleday<-ifelse(flowdaysite$dateTime==flowdataSite$dateTime,"yes","no")

ggplot(flowdataSite,aes(dateTime,X_00060_00003))+
  geom_point(color="blue")+
  geom_point(flowdaysite,aes(dateTime,X_00060_00003,color="red"))

#########Chem Line Plot###########################
title<-c("Phosphorus mg/L","Total Nitrogen mg/L","Orthophosphate mg/L",
         "Chlorophyll a, periphyton, mg/m2")

ggplot(samples,aes(Collect_Date,samples[,6],
                   colour=factor(Stream),
                   group=factor(Stream,
                                levels=c("Salmon River", "Norwalk River","Pequabuck River"))))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_x_date(date_labels = "%m-%Y")+
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  labs(colour="Stream",x=NULL,y=NULL,title=title[1])+
  theme(legend.position="bottom")

ggsave("plots/TPpoint.tiff",width=4.5,height=3,units="in")

#########Chem Line Plot###########################
title<-c("Average Monthly Temperature (C)","Median 2 Week Flow Prior to Sample Date (cfsm)",
         "Percent Canopy Cover")
##columns 5,14,10

ggplot(samples,aes(Collect_Date,samples[,10],
                   colour=factor(Stream),
                   group=factor(Stream,
                                levels=c("Salmon River", "Norwalk River","Pequabuck River"))))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_x_date(date_labels = "%m-%Y")+
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  labs(colour="Stream",x=NULL,y=NULL,title=title[3])+
  theme(legend.position="bottom")

ggsave("plots/canopypoint.tiff",width=4.5,height=3,units="in")


######Chem Box Plot###########################################

p1<-ggplot(samples,aes(Stream,OrthoP))+
  geom_boxplot()+
  scale_y_log10()+
  labs(title="TP")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())



# #####MULTI-Plot Function##############
# #######################################
# 
# #call this with p1,p,2,... cols=4
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   plots <- c(list(...), plotlist)
#   numPlots = length(plots)
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   if (numPlots==1) {
#     print(plots[[1]])
#   } else {
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }
# #######################################
# pCd<-unique(data$parm_cd)
# pNm<- c("Organic nitrogen, wu, mg/l","Ammonia, wf, mg/l as N","Nitrite, wf, mg/l as N",
#   "Nitrate, wf, mg/l as N","NO3+NO2, wf, mg/l as N","Orthophosphate, wf, mg/l",
#   "Phosphorus, wu, mg/l as P","Orthophosphate, wf, mg/l as P","Total nitrogen, wu, mg/l")
# n<-as.data.frame(cbind(pCd,pNm))
# data$site_no<- factor(data$site_no,levels=c("01193500","01209700","01189000"))
# 
# 
# for (i in 1:dim(n)[1]){
# d<- data[data$parm_cd==n[i,1],]
# 
# assign(paste0("p",i), 
#        ggplot(d,aes(site_no,result_va))+
#           geom_boxplot()+
#           scale_y_log10()+
#           labs(title=n[i,2])+
#           theme(axis.title.x=element_blank(),axis.title.y=element_blank()))
# 
# }
# 
# pdf("NutrientRanges.pdf",width=11,height=8)
# multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,cols=3)
# dev.off()
# 
# ###########Scatterplots by site###################
# ###########**Specify site below**#################
# for (i in 1:dim(n)[1]){
#   d<- data[data$parm_cd==n[i,1],]
#   d<-d[d$site_no=='01189000',]
#   
#   assign(paste0("p",i), 
#          ggplot(d,aes(sample_dt,result_va))+
#            geom_point()+
#            geom_smooth(se=TRUE,colour="black")+
#            labs(title=n[i,2])+
#            theme(axis.title.x=element_blank(),axis.title.y=element_blank()))
#   
# }
# 
# pdf("NutrientScatterplot_01189000.pdf",width=11,height=8)
# multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,cols=3)
# dev.off()
# 
# 
# 
