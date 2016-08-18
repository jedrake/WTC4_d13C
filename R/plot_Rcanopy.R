#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to process the canopy respiration
#  data and merge the rates with the isotopic composition estimates
#  as an estimate of the total amount of label returned by canopy R.
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#-- load the required libraries
source("R/loadLibraries.R")
#------------------------------------------------------------------------




#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- Process all of the gas exchange data to estimate the canopy respiration RATE.
#  This is going to take considerable work...
#------------------------------------------------------------------------
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- Fit all the keeling plots to estimate the d13C of Rcanopy.

#- get the data
isodat <- getIso()
isoCR <- subset(isodat,type=="CR")

#- break into a list, fit the keeling plots. "Keeling_int" reflects the estimated 
#  isotopic composition of soil respiration
isoCR.l <- split(isoCR,as.factor(paste(isoCR$chamber,isoCR$Batch.DateTime,sep="-")))
CR_KP <- lapply(isoCR.l,FUN=fitKeeling)
CR_KP.df <- do.call(rbind,CR_KP)


#- plot chambers and then treatment averages

#------------------------------------------------------------------------
#- plot chambers over time 
plotBy(Keeling_int~Batch.DateTime|chamber,data=CR_KP.df,legend=F)
plotBy(r2~Batch.DateTime|chamber,data=CR_KP.df,legend=F,ylim=c(0,1))

#- calculate treatment averages over time, exclude C01 and C02. Note we haven't yet measured C01 or C02,
#  so this isn't really necessary. We did measure some pre-treatment, and we may measure them again
isoCR2 <- CR_KP.df[!(CR_KP.df$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoCR.m2 <- summaryBy(Keeling_int~T_treatment+Batch.DateTime,data=isoCR2,FUN=c(mean,standard.error))


#---
#- plot treatment averages
pdf(file=paste("Output/Rcanopy_d13C_",Sys.Date(),".pdf",sep=""))
plotBy(Keeling_int.mean~Batch.DateTime|T_treatment,data=isoCR.m2,col=c("blue","red"),pch=16,ylim=c(-30,800),legend="F",
       panel.first=adderrorbars(x=isoCR.m2$Batch.DateTime,y=isoCR.m2$Keeling_int.mean,
                                SE=isoCR.m2$Keeling_int.standard.error,direction="updown"))
abline(h=0)
title(main="Canopy R")
#- add a legend
legend("topright",pch=16,col=c("blue","red"),legend=c("Ambient","Warmed"))
dev.off()
#------------------------------------------------------------------------
#------------------------------------------------------------------------