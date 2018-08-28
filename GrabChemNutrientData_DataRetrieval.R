setwd("P:/Projects/GitHub_Prj/SeasonalProductivity")

samples<-read.csv("data/sites.csv",header=TRUE)

library(dataRetrieval)
library(plyr)
library(ggplot2)
library(grid)
library(reshape2)


SID<- c("01189000","01209700","01193500")
Cd<-c("00605","00608","00613","00618","00631","00660","00665","00671","62855")
CdFlow<-"00060"
sdate<-"2016-04-01"
edate<-"2016-10-01"


data <- readNWISqw(siteNumbers = SID,
                     parameterCd = Cd,
                     startDate = sdate,
                     endDate = edate)

flowdata<-readNWISdata(siteNumbers = SID,
                       parameterCd = CdFlow,
                       startDate = sdate,
                       endDate = edate)

write.csv(data,"GrabChem_NutrientData.csv")
write.csv(flowdata,"flowdata.csv")


dsum<- data%>%
          group_by(site_no,parm_cd)%>%
          summarize(count=n(),avg=mean(result_va),
                min=min(result_va),max=max(result_va),
                med=median(result_va))


#########Flow Plot###########################

ggplot(flowdata[which(flowdata$site_no=='01193500'),],aes(dateTime,X_00060_00003))+
  geom_line()







#####MULTI-Plot Function##############
#######################################

#call this with p1,p,2,... cols=4
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#######################################
pCd<-unique(data$parm_cd)
pNm<- c("Organic nitrogen, wu, mg/l","Ammonia, wf, mg/l as N","Nitrite, wf, mg/l as N",
  "Nitrate, wf, mg/l as N","NO3+NO2, wf, mg/l as N","Orthophosphate, wf, mg/l",
  "Phosphorus, wu, mg/l as P","Orthophosphate, wf, mg/l as P","Total nitrogen, wu, mg/l")
n<-as.data.frame(cbind(pCd,pNm))
data$site_no<- factor(data$site_no,levels=c("01193500","01209700","01189000"))


for (i in 1:dim(n)[1]){
d<- data[data$parm_cd==n[i,1],]

assign(paste0("p",i), 
       ggplot(d,aes(site_no,result_va))+
          geom_boxplot()+
          scale_y_log10()+
          labs(title=n[i,2])+
          theme(axis.title.x=element_blank(),axis.title.y=element_blank()))

}

pdf("NutrientRanges.pdf",width=11,height=8)
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,cols=3)
dev.off()

###########Scatterplots by site###################
###########**Specify site below**#################
for (i in 1:dim(n)[1]){
  d<- data[data$parm_cd==n[i,1],]
  d<-d[d$site_no=='01189000',]
  
  assign(paste0("p",i), 
         ggplot(d,aes(sample_dt,result_va))+
           geom_point()+
           geom_smooth(se=TRUE,colour="black")+
           labs(title=n[i,2])+
           theme(axis.title.x=element_blank(),axis.title.y=element_blank()))
  
}

pdf("NutrientScatterplot_01189000.pdf",width=11,height=8)
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,cols=3)
dev.off()
