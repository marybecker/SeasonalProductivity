library(dataRetrieval)

SID<- c("01189000","01209700","01193500")
Cd<-c("00605","00608","00613","00618","00631","00660","00665","00671","62855")
sdate<-"2016-01-01"
edate<-"2017-01-01"


data <- readNWISqw(siteNumbers = SID,
                     parameterCd = Cd,
                     startDate = sdate,
                     endDate = edate)

write.csv(data,"GrabChem_NutrientData.csv")