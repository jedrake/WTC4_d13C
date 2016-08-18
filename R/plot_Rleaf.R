#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to read and look at the data for the isotopic
#  composition of leaf respiration.
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#-- load the required libraries
source("R/loadLibraries.R")
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#- get the data
isodat <- getIso()
isoLR <- subset(isodat,type=="LR")
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#- plot chambers over time 
plotBy(d13C~Batch.DateTime|chamber,data=isoLR,legend=F)

#- calculate treatment averages over time, exclude C01 and C02
isoLR2 <- isoLR[!(isoLR$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoLR.m2 <- summaryBy(d13C~T_treatment+Batch.DateTime,data=isoLR2,FUN=c(mean,standard.error))

#- get averages of CO1 and CO2 to add
isoLR3 <- isoLR[(isoLR$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoLR.m3 <- summaryBy(d13C~Batch.DateTime,data=isoLR3,FUN=c(mean,standard.error))

#- plot treatment averages
plotBy(d13C.mean~Batch.DateTime|T_treatment,data=isoLR.m2,col=c("blue","red"),pch=16,ylim=c(-30,1300),legend="F",
       panel.first=adderrorbars(x=isoLR.m2$Batch.DateTime,y=isoLR.m2$d13C.mean,
                                SE=isoLR.m2$d13C.standard.error,direction="updown"))
#- add the non-labeled trees
points(d13C.mean~Batch.DateTime,data=isoLR.m3,pch=16,col="black",
         panel.first=adderrorbars(x=isoLR.m3$Batch.DateTime,y=isoLR.m3$d13C.mean,
                                SE=isoLR.m3$d13C.standard.error,direction="updown"))
abline(h=0)

#- add a legend
legend("topright",pch=16,col=c("blue","red","black"),legend=c("Ambient","Warmed","Unlabeled"))
#------------------------------------------------------------------------