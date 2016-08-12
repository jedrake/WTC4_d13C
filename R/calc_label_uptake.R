#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to calculate the amount of 13C label taken
#  up by each chamber during the labeling even on 2016-08-05 (5 Aug)
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
labelday$chamber <- factor(labelday$chamber)
#------------------------------------------------------------------------





#------------------------------------------------------------------------
#--- Isotope data

#- get all of the isotope data (note, this function needs edition when all data have been entered)
isodat <- getIso()

#- pull out just the labeling data, calculate atom percent C
labeling <-subset(isodat,type=="LAB")
ARC <- 0.0111803 # absolute ratio of VPDB. See http://www.isotope.uottawa.ca/guides/guides-atom-percent-en.html
labeling$AP <- with(labeling,(100*ARC*(d13C/1000+1)) / (1+ARC*(d13C/1000+1)) ) 

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

#- seconds 
times <- data.frame(seq(min(starttime), max(endtime), "min"))
names(times)[1] <- "Collection.DateTime" 

#- list, with one element per chamber
iso.l <- split(labeling,labeling$chamber)
for(i in 1:length(iso.l)){
  
  #- merge isotope data with dummy vector of times
  isodat <- iso.l[[i]][,c("Collection.DateTime","d13C")]
  isodat2 <- merge(isodat,times,by="Collection.DateTime",all.y=T)
  
  #- time series of isotope data
  zoo.iso <- zoo(isodat2,order.by=isodat2$Collection.DateTime)
  zoo.iso2 <- merge(zoo.iso,zoo.times)
  
}
#------------------------------------------------------------------------