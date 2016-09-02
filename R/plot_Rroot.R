#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to read and look at the data for the isotopic
#  composition of root respiration.
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#-- load the required libraries
source("R/loadLibraries.R")
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#- get the data
isodat <- getIso()
isoRR <- subset(isodat,type=="RR")
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#- plot chambers over time 
plotBy(d13C~Batch.DateTime|chamber,data=isoRR,legend=F)

#- calculate treatment averages over time, exclude C01 and C02
isoRR2 <- isoRR[!(isoRR$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoRR.m2 <- summaryBy(d13C~T_treatment+Batch.DateTime,data=isoRR2,FUN=c(mean,standard.error))

#- get averages of CO1 and CO2 to add
isoRR3 <- isoRR[(isoRR$chamber %in% c("C01","C02")),] # include C01 and C02
isoRR.m3 <- summaryBy(d13C~Batch.DateTime,data=isoRR3,FUN=c(mean,standard.error))

#---
#- plot treatment averages
pdf(file=paste("Output/Rroot_d13C_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(d13C.mean~Batch.DateTime|T_treatment,data=isoRR.m2,col=c("blue","red"),pch=16,ylim=c(-30,350),
       xlab="",ylab="",axes=F,legend="F",
       panel.first=adderrorbars(x=isoRR.m2$Batch.DateTime,y=isoRR.m2$d13C.mean,
                                SE=isoRR.m2$d13C.standard.error,direction="updown"))
#- add the non-labeled trees
points(d13C.mean~Batch.DateTime,data=isoRR.m3,pch=16,col="black",
       panel.first=adderrorbars(x=isoRR.m3$Batch.DateTime,y=isoRR.m3$d13C.mean,
                                SE=isoRR.m3$d13C.standard.error,direction="updown"))
abline(h=0)
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-08-05 00:00:00",tz="UTC"),
                                  to=as.POSIXct("2016-09-01 00:00:00",tz="UTC"),by="day"),las=2)
title(ylab=expression(paste(R[root]," ",delta^{13}, "C (\u2030)")),
      main="Root respiration",cex.lab=1.5)

#- add a legend
legend("topright",pch=16,col=c("blue","red","black"),legend=c("Ambient","Warmed","Unlabeled"))
dev.off()
#------------------------------------------------------------------------