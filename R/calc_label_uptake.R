#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to calculate the amount of 13C label taken
#  up by each chamber during the labeling event on 2016-08-05 (5 Aug)
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#- source required libraries
source("R/loadLibraries.R")




#------------------------------------------------------------------------
#----- Flux data
#- read in the flux data (quicker than read.csv for large files).
fluxes <- getFluxes()

#- extract just the labeling period from the flux data. Note that this varies a bit for each chamber, so could be improved.
starttime <- as.POSIXct("2016-08-05 12:50:00",tz="UTC")
endtime <-  as.POSIXct("2016-08-05 17:15:00",tz="UTC")

#-- extract only the labeled chambers on the labeling day
labelday <- subset(fluxes,
                   DateTime>starttime & DateTime<endtime & chamber %in% c("C04","C05","C06","C07","C08","C09"),
)
labelday$chamber <- factor(labelday$chamber) #- drop the levels associated with un-labeled chambers
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#--- Isotope data

#- get all of the isotope data (note, this function needs edition when all data have been entered)
isodat <- getIso()

#- pull out just the labeling data, calculate atom percent C
labeling <-subset(isodat,type=="LAB")
labeling$AP <- getAP(labeling$d13C) # AP has units of %

#------------------------------------------------------------------------




#------------------------------------------------------------------------
#- plot par, co2 uptake, and canopy isotopic composition. Note that ambient chambers have lower Co2 uptake.
windows(100,60);par(mfrow=c(4,1),mar=c(1,6,1,1),oma=c(5,2,0,0),cex.lab=2)
#palette(c("black",brewer.pal(5,"Spectral")))
palette(c("red","blue","orange","darkgrey","brown","lightblue"))

plotBy(PAR~DateTime|chamber,data=labelday,type="o",ylab="PAR",legendwhere="topright",pch=16,cex=1.5)   ;abline(h=0)
plotBy(FluxCO2~DateTime|chamber,data=labelday,ylim=c(-0.1,0.3),,pch=16,cex=1.5,type="o",ylab="CO2 Flux (mmol s-1)",legend=F)   ;abline(h=0)
plotBy(d13C~Collection.DateTime|chamber,data=labeling,pch=16,cex=1.5,type="b")
plotBy(AP~Collection.DateTime|chamber,data=labeling,pch=16,cex=1.5,type="b")

#------------------------------------------------------------------------



#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- Break into a list, interpolate into a minutely-resolved and interpolated
#  time series, and calculate the amount of label assimilated by each chamber
library(zoo)

#- vector of minutes between measurements, made into a dataframe for merging 
times <- expand.grid(Collection.DateTime=seq(min(starttime), max(endtime), "min"),chamber=c("C04","C05","C06","C07","C08","C09"))

#- merge the isotope data with the minutes
iso.dat <- merge(labeling,times,by=c("Collection.DateTime","chamber"),all=T)[,c("Collection.DateTime","chamber","d13C","AP")]
#- merge the flux data
iso.flux <- merge(labelday,iso.dat,by.x=c("DateTime","chamber"),
                  by.y=c("Collection.DateTime","chamber"),all=T)[,c("DateTime","chamber","d13C","AP","FluxCO2")]



#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- Interpolate data into minutely-resolved numbers
#- now we have a dataframe with minutely-resolved data, ready for interpolation
#- list, with one element per chamber
iso.l <- split(iso.flux,iso.flux$chamber)
iso.interp.l <- list()
for(i in 1:length(iso.l)){
  
  print(i)
  
  #- extract the data to work with, for convenience
  isodat <- iso.l[[i]]
  
  #- replace <NA> values with NA
  isodat[which(is.na(isodat$d13C)),"d13C"] <- NA
  
  
  
  
  isodat2 <- isodat[,c("d13C","AP","FluxCO2")]
  
  #- time series of isotope data
  zoo.iso <- zoo(isodat2)
 
  #- NA approximation by linear interpolation
  iso.interp <- zoo::na.approx(zoo.iso,na.rm=F) 
  
  #- put it back to a normal dataframe (get out of zoo format), add informative variables back
  iso.interp2 <- data.frame(coredata(iso.interp))
  iso.interp2$DateTime <- isodat$DateTime
  iso.interp2$chamber <- isodat$chamber
 
  
  
  
  #- much of the flux data are missing for chamber 8. I'm not sure what happened here, but I think we lost the local IRGA
  #  and had to reset the PSOIX module (and thus lost diffP and CO2L) during the labeling. This means that the system
  #  can't calculate the fluxes as it normally does.
  #- As a placeholder solution, I'm estimating the C08 fluxes using the C06 fluxes, which is defensible.
  if (i == 5){
    c6data <- iso.interp.l[[3]]
    startfill <- as.POSIXct("2016-08-05 14:17:00",tz="UTC")
    stopfill <- as.POSIXct("2016-08-05 16:00:00",tz="UTC")
    
    tofill <- which(iso.interp2$DateTime >= startfill & iso.interp2$DateTime <= stopfill)
    iso.interp2$FluxCO2[tofill] <- c6data$FluxCO2[tofill]
  }
  
  #- assign the dataframe to the list that was "preallocated" to catch the output
  iso.interp.l[[i]] <- iso.interp2
  
}
#- merge the lists back into a big dataframe
isodat.int <- do.call(rbind,iso.interp.l)
isodat.int <- isodat.int[complete.cases(isodat.int),] # get rid of NA's
#------------------------------------------------------------------------
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- calculate and sum 13C uptake, convert from mmol 13CO2 s-1 to mmol 13CO2 min-1
isodat.int$uptake13C <- with(isodat.int,FluxCO2*AP/100*60) 

#- zero fill any negative numbers
negs <- which(isodat.int$uptake13C < 0)
isodat.int$uptake13C[negs] <- 0

isodat.sums <- summaryBy(uptake13C~chamber,data=isodat.int,FUN=sum,keep.names=T)
isodat.sums$T_treatment <- factor(rep(c("elevated","ambient"),3),levels=c("ambient","elevated"))

#- convert mmol 13CO2 to g 13C
isodat.sums$uptake13C_g <- isodat.sums$uptake13C/1000*13
names(isodat.sums)[2] <- c("uptake13CO2_mmol")

#- rearrange to make it pretty
isodat.sums <- isodat.sums[c(3,1,2,4)]
#------------------------------------------------------------------------
#------------------------------------------------------------------------