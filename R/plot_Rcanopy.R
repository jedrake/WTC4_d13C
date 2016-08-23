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
#  This is going to take considerable work... Let's have a go.
flux.all <- getFluxes()

#- subset to just teh data after the labeling day
flux.lab <- subset(flux.all,as.Date(DateTime)>=as.Date("2016-08-05"))

#- clean up the flux data
flux <- cleanFluxData(flux.lab)


#- Calculate the median respiration rates for each chamber on each date (median is less affected by outliers)
flux.night <- summaryBy(FluxCO2+Tair_al+CO2L~Date_night+chamber+T_treatment,data=subset(flux,PAR<5),
                        FUN=median,keep.names=T,na.rm=T)

#- convert from mmol CO2 per second to mmol per hr
flux.night$R_mmol <- flux.night$FluxCO2*-1*60*60

#- average by treatment
flux.night.m <- summaryBy(R_mmol+Tair_al~Date_night+T_treatment,data=flux.night,
                          FUN=c(mean,standard.error))

pdf(file=paste("Output/Rcanopy_rate_",Sys.Date(),".pdf",sep=""))
par(mar=c(5,6,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1))
plotBy(R_mmol.mean~Date_night|T_treatment,data=flux.night.m,col=c("blue","red"),pch=16,ylim=c(0,50),legend="F",
       type="o",cex=2,xlab="Date",ylab="Rcanopy (mmol CO2 hr-1)",main="Rcanopy",
       panel.first=adderrorbars(x=flux.night.m$Date_night,y=flux.night.m$R_mmol.mean,
                                SE=flux.night.m$R_mmol.standard.error,direction="updown"))

plotBy(R_mmol.mean~Tair_al.mean|T_treatment,data=flux.night.m,col=c("blue","red"),pch=16,ylim=c(0,50),legend="F",
       type="p",cex=2,xlab="Tair (deg C)",ylab="R (mmol CO2 hr-1)",
       panel.first=adderrorbars(x=flux.night.m$Tair_al.mean,y=flux.night.m$R_mmol.mean,
                                SE=flux.night.m$R_mmol.standard.error,direction="updown"))
dev.off()

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


#------------------------------------------------------------------------
#- plot chambers and then treatment averages

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
par(mar=c(6,7,1,4))
plotBy(Keeling_int.mean~Batch.DateTime|T_treatment,data=isoCR.m2,col=c("blue","red"),pch=16,ylim=c(-30,800),
       legend="F",axes=F,xlab="",ylab="",
       panel.first=adderrorbars(x=isoCR.m2$Batch.DateTime,y=isoCR.m2$Keeling_int.mean,
                                SE=isoCR.m2$Keeling_int.standard.error,direction="updown"))
abline(h=0)
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-08-05 00:00:00",tz="UTC"),
                                  to=as.POSIXct("2016-08-17 00:00:00",tz="UTC"),by="day"),las=2)
title(ylab=expression(paste(R[canopy]," ",delta^{13}, "C (\u2030)")),
      main="Canopy respiration",cex.lab=1.5)

#- add a legend
legend("topright",pch=16,col=c("blue","red"),legend=c("Ambient","Warmed"))
dev.off()
#------------------------------------------------------------------------
#------------------------------------------------------------------------







#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- merge the data regarding daily mean canopy respiration rates and isotopic 
#  composition. Estimate the total amount of label respired aboveground.


#- subtract one day from the date of night-time measurements after midnight but before 5am
#This makes it easier to calculate nightly-mean respiration rates
CR_KP.df$Date <- as.Date(CR_KP.df$Batch.DateTime)
CR_KP.df$Date_night <- CR_KP.df$Date
toadd <- which(hour(CR_KP.df$Batch.DateTime)<=5)
CR_KP.df$Date_night[toadd] <- CR_KP.df$Date[toadd]-1

#- daily average canopy respiration isotopic composition
CR_KP_m <- summaryBy(Keeling_int~chamber+T_treatment+Date_night,data=CR_KP.df,FUN=mean,keep.names=T)


#- daily average canopy respiration rate
head(flux.night)
flux.night <- subset(flux.night,chamber %in% c("C04","C05","C06","C07","C08","C09"))

#- mergedatasets
Rcanopy <- merge(CR_KP_m,flux.night,by=c("Date_night","chamber","T_treatment"),all.y=T)
Rcanopy$chamber <- factor(Rcanopy$chamber)

#- calculate the amount of label respired, in units of mg 13C.
# Assumes the night is 13 hours long (sunrise at 6:30am, sunset at 5:30pm)



#------------------------------------------------------------------------
#------------------------------------------------------------------------



#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- loop over each chamber, gapfill missing data, calculate cumulative sum
#- get and plot the cumulative sum of 13C respiration
Rcanopy.l <- split(Rcanopy,Rcanopy$chamber)
for (i in 1:length(Rcanopy.l)){
  #- subset to just dates I want
  Rcanopy.l[[i]] <- subset(Rcanopy.l[[i]],Date_night %in% seq.Date(from=as.Date("2016-08-05"),
                                                                   to=as.Date("2016-08-16"),by=1))
  
  #- gapfill
  Rcanopy.l[[i]]$Keeling_int <- na.approx(Rcanopy.l[[i]]$Keeling_int)
  
  #- calculate 13C respired
  Rcanopy.l[[i]]$AP <- getAP(Rcanopy.l[[i]]$Keeling_int) # get the atom percent 13C from per mill data
  Rcanopy.l[[i]]$AP_natural <- getAP(-30)         # get the atom percent 13C of "background" respiration
  Rcanopy.l[[i]]$Rcanopy_13C <- with(Rcanopy.l[[i]],R_mmol*(AP-AP_natural)/100*13/1000*13*1000) # convert to units of mg 13C of EXCESS 13C
  
  #- cumulative sum
  Rcanopy.l[[i]]$Rcanopy_13C_cumsum <- cumsum(Rcanopy.l[[i]]$Rcanopy_13C)
}
Rcanopy2 <- do.call(rbind,Rcanopy.l)

#------------------------------------------------------------------------
#------------------------------------------------------------------------



#- calculate treatment averages over time, exclude C01 and C02. Note we haven't yet measured C01 or C02,
#  so this isn't really necessary. We did measure some pre-treatment, and we may measure them again
isoCS <- summaryBy(Rcanopy_13C+Rcanopy_13C_cumsum~T_treatment+Date_night,data=Rcanopy2,FUN=c(mean,standard.error))


#---
#- plot treatment averages
pdf(file=paste("Output/Rcanopy_d13C_sums_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(Rcanopy_13C_cumsum.mean~Date_night|T_treatment,data=isoCS,col=c("blue","red"),pch=16,
       legend="F",axes=F,xlab="",ylab="",ylim=c(0,100),
       panel.first=adderrorbars(x=isoCS$Date_night,y=isoCS$Rcanopy_13C_cumsum.mean,
                                SE=isoCS$Rcanopy_13C_cumsum.standard.error,direction="updown"))
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-08-05"),
                                  to=as.Date("2016-08-17"),by="day"),las=2)
title(ylab=expression(R[canopy]~(mg~13*C~excess)),main="Cumulative sum, Rcanopy",cex.lab=1.5)

#- add a legend
legend("topleft",pch=16,col=c("blue","red"),legend=c("Ambient","Warmed"))
dev.off()

#- sum, then average across treatments
iso.s <- summaryBy(Rcanopy_13C~chamber+T_treatment,data=Rcanopy2,FUN=c(sum),keep.names=T)
iso.m <- summaryBy(Rcanopy_13C~T_treatment,data=iso.s,FUN=mean)

#------------------------------------------------------------------------
#------------------------------------------------------------------------
