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
isoLR$logd13C <- log10(isoLR$d13C+50)
isoLR2 <- isoLR[!(isoLR$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoLR.m2 <- summaryBy(d13C+logd13C~T_treatment+Batch.DateTime,data=isoLR2,FUN=c(mean,standard.error))

#- get averages of CO1 and CO2 to add
isoLR3 <- isoLR[(isoLR$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoLR.m3 <- summaryBy(d13C~Batch.DateTime,data=isoLR3,FUN=c(mean,standard.error))
isoLR.m3 <- summaryBy(d13C+logd13C~Batch.DateTime,data=isoLR3,FUN=c(mean,standard.error))

#---
#- plot treatment averages
pdf(file=paste("Output/Rleaf_d13C_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(d13C.mean~Batch.DateTime|T_treatment,data=isoLR.m2,col=c("blue","red"),pch=16,ylim=c(-30,1300),
       xlab="",ylab="",axes=F,legend="F",
       panel.first=adderrorbars(x=isoLR.m2$Batch.DateTime,y=isoLR.m2$d13C.mean,
                                SE=isoLR.m2$d13C.standard.error,direction="updown"))
#- add the non-labeled trees
points(d13C.mean~Batch.DateTime,data=isoLR.m3,pch=16,col="black",
         panel.first=adderrorbars(x=isoLR.m3$Batch.DateTime,y=isoLR.m3$d13C.mean,
                                SE=isoLR.m3$d13C.standard.error,direction="updown"))
abline(h=0)
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-08-05 00:00:00",tz="UTC"),
                                 to=as.POSIXct("2016-00-01 00:00:00",tz="UTC"),by="day"),las=2)
title(ylab=expression(paste(R[leaf]," ",delta^{13}, "C (\u2030)")),
      main="Leaf respiration",cex.lab=1.5)

#- add a legend
legend("topright",pch=16,col=c("blue","red","black"),legend=c("Ambient","Warmed","Unlabeled"))
dev.off()
#------------------------------------------------------------------------






#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- plot log-transformed data

#---
#- plot treatment averages
pdf(file=paste("Output/Rleaf_d13C_log_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(logd13C.mean~Batch.DateTime|T_treatment,data=isoLR.m2,col=c("blue","red"),pch=16,ylim=c(1,3.5),
       xlab="",ylab="",axes=F,legend="F",
       panel.first=adderrorbars(x=isoLR.m2$Batch.DateTime,y=isoLR.m2$logd13C.mean,
                                SE=isoLR.m2$logd13C.standard.error,direction="updown"))
#- add the non-labeled trees
points(logd13C.mean~Batch.DateTime,data=isoLR.m3,pch=16,col="black",
       panel.first=adderrorbars(x=isoLR.m3$Batch.DateTime,y=isoLR.m3$logd13C.mean,
                                SE=isoLR.m3$logd13C.standard.error,direction="updown"))
abline(h=0)
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-08-05 00:00:00",tz="UTC"),
                                  to=as.POSIXct("2016-09-01 00:00:00",tz="UTC"),by="day"),las=2)
title(ylab=expression(paste(log[10]~"("~R[leaf]," ",delta^{13}, "C (\u2030)",+50~")")),
      main="Leaf respiration",cex.lab=1.5)

#- add a legend
legend("topright",pch=16,col=c("blue","red","black"),legend=c("Ambient","Warmed","Unlabeled"))
dev.off()
#------------------------------------------------------------------------