#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to process the viasala data for soil respiration
#  and return estimates of the soil respiration rate
#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#-- load the required libraries
source("R/loadLibraries.R")
#------------------------------------------------------------------------






#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- Data processing

#- find all of the soil respiration files, measured post-labeling
files <- list.files(path="Data/Viasala/Post-labeling/",pattern="SR")

fits.all <- list() #- preallocate a list to "catch" the output of the fits
#- look over each file, read in the data
for (i in 1:length(files)){
  
  #- read in the data, do some housekeeping
  dat <- read.csv(file=paste("Data/Viasala/Post-labeling/",files[i],sep=""))[,1:2]
  names(dat) <- c("DateTime","CO2")
  
  #- does dat$DateTime have an AM or PM in it?
  #AMflag <- any(grep("AM",dat$DateTime))# | grep("PM",dat$DateTime))
  #PMflag <- any(grep("PM",dat$DateTime))# | grep("PM",dat$DateTime))
  
  dat$DateTime <- as.POSIXct(dat$DateTime,format="%d/%m/%Y %H:%M:%S",tz="UTC")
  mintime <- dat$DateTime[1]
  dat <- subset(dat,DateTime>mintime+60) # remove the first minute, as many files have bad data before leuvers close fully
  dat$dTime <- as.numeric(difftime(dat$DateTime,mintime,units="mins"))
  
  #- extract the chamber number and DateTime from the file name.
  cc <- base::strsplit(as.character(files[i]),split="-")
  hour <- cc[[1]][1]
  chamber <- cc[[1]][3]
  year <- cc[[1]][4]
  month <- cc[[1]][5]
  day <- substr(cc[[1]][6],start=1,stop=2)
  DateTime <- as.POSIXct(paste(year,month,day,hour,sep="-"),format="%Y-%m-%d-%H",tz="UTC")
  
  #- actually do the fitting (call to function in processingFunction.R), assign output as dataframe to preallocated list
  fit.temp <- fitRsoil(dat=dat,chamber=chamber,plotson=F)
  outDF <- data.frame("chamber"=chamber,"DateTime"=DateTime,"Rsoil"=fit.temp[1],"pseudoR2"=fit.temp[2])
  fits.all[[i]] <- outDF
  
}

#- compile the output. Note that "Rsoil" has units of umol CO2 m-2 s-1.
outdat <- do.call(rbind,fits.all)
outdat$chamber <- factor(outdat$chamber,levels=c("C01","C02","C04","C05","C06","C07","C08","C09"))
link <- data.frame(chamber=levels(outdat$chamber),
                   T_treatment=c("ambient","elevated","elevated","ambient","elevated","ambient","elevated","ambient"))
outdat <- merge(outdat,link,by="chamber")
#------------------------------------------------------------------------
#------------------------------------------------------------------------




#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- plot output

#- plot Rsoil for each chamber
pdf(file=paste("Output/Rsoil_rates_",Sys.Date(),".pdf",sep=""))
palette(c("green","darkred","red","blue","orange","darkgrey","brown","lightblue"))
plotBy(Rsoil~DateTime|chamber,data=outdat,type="o",pch=16)

#- plot Rsoil for each treatment
outdat.m <- summaryBy(Rsoil~T_treatment+DateTime,data=subset(outdat,
                                                             chamber %in% c("C04","C05","C06","C07","C08","C09")),
                      FUN=c(mean,standard.error))

plotBy(Rsoil.mean~DateTime|T_treatment,data=outdat.m,col=c("blue","red"),pch=16,ylim=c(0,5),type="o",
       panel.first=adderrorbars(x=outdat.m$DateTime,y=outdat.m$Rsoil.mean,
                                SE=outdat.m$Rsoil.standard.error,direction="updown"),
       ylab="Rsoil (umol CO2 m-2 s-1)")
dev.off()
#------------------------------------------------------------------------
#------------------------------------------------------------------------











#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- okay, now fit all the keeling plots to estimate the d13C of Rsoil.

#- get the data
isodat <- getIso()
isoRsoil <- subset(isodat,type=="SR")

#- break into a list, fit the keeling plots. "Keeling_int" reflects the estimated 
#  isotopic composition of soil respiration
isoRsoil.l <- split(isoRsoil,as.factor(paste(isoRsoil$chamber,isoRsoil$Batch.DateTime,sep="-")))
RsoilKP <- lapply(isoRsoil.l,FUN=fitKeeling)
RsoilKP.df <- do.call(rbind,RsoilKP)


#-- pull out two examples to fit
ids <- which(isoRsoil$chamber=="C08" & isoRsoil$Batch.DateTime %in% c(as.POSIXct("2016-08-06 00:00:00",tz="UTC"),c(as.POSIXct("2016-08-06 16:00:00",tz="UTC"))))
Rsoil.examples <- isoRsoil[ids,]
Rsoil.examples$CO2i <- 1/Rsoil.examples$CO2

windows();par(mar=c(6,6,1,1),cex.lab=2)
plotBy(d13C~CO2i|as.factor(Batch.DateTime),data=Rsoil.examples,legend=F,col=c("black","red"),pch=16,
       xlab="1/[CO2]",ylab=expression(paste(delta^{13}, "C (\u2030)")))
lm1 <- lm(d13C~CO2i,data=Rsoil.examples[c(1,4,6),])
lm2 <- lm(d13C~CO2i,data=Rsoil.examples[c(2,3,5),])
abline(lm1,col="black")
abline(lm2,col="red")


#- plot chambers and then treatment averages

#------------------------------------------------------------------------
#- plot chambers over time 
plotBy(Keeling_int~Batch.DateTime|chamber,data=RsoilKP.df,legend=F)
plotBy(r2~Batch.DateTime|chamber,data=RsoilKP.df,legend=F,ylim=c(0,1))

#- calculate treatment averages over time, exclude C01 and C02
isoSR2 <- RsoilKP.df[!(RsoilKP.df$chamber %in% c("C01","C02")),] # exclude C01 and C02
isoSR.m2 <- summaryBy(Keeling_int~T_treatment+Batch.DateTime,data=isoSR2,FUN=c(mean,standard.error))

#- get averages of CO1 and CO2 to add
isoSR3 <- RsoilKP.df[(RsoilKP.df$chamber %in% c("C01","C02")),] # only include C01 and C02
isoSR.m3 <- summaryBy(Keeling_int~Batch.DateTime,data=isoSR3,FUN=c(mean,standard.error))


#---
#- plot treatment averages
pdf(file=paste("Output/Rsoil_d13C_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(Keeling_int.mean~Batch.DateTime|T_treatment,data=isoSR.m2,col=c("blue","red"),pch=16,ylim=c(-30,200),
       legend="F",axes=F,xlab="",ylab="",
       panel.first=adderrorbars(x=isoSR.m2$Batch.DateTime,y=isoSR.m2$Keeling_int.mean,
                                SE=isoSR.m2$Keeling_int.standard.error,direction="updown"))
#- add the non-labeled trees
points(Keeling_int.mean~Batch.DateTime,data=isoSR.m3,pch=16,col="black",
       panel.first=adderrorbars(x=isoSR.m3$Batch.DateTime,y=isoSR.m3$Keeling_int.mean,
                                SE=isoSR.m3$Keeling_int.standard.error,direction="updown"))
abline(h=0)
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-08-05 00:00:00",tz="UTC"),
                                  to=as.POSIXct("2016-08-19 00:00:00",tz="UTC"),by="day"),las=2)
title(ylab=expression(paste(R[soil]," ",delta^{13}, "C (\u2030)")),
      main="Soil respiration",cex.lab=1.5)

#- add a legend
legend("topright",pch=16,col=c("blue","red","black"),legend=c("Ambient","Warmed","Unlabeled"))
dev.off()
#------------------------------------------------------------------------
#------------------------------------------------------------------------




#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- okay then, let's merge the estimates of teh Rsoil flux
#  and isotopic composition to estimate the total amount of 13C
#  label respired belowground.
head(isoSR2)
head(outdat)

#- daily average Rsoil rate
outdat$Date <- as.Date(outdat$DateTime)
Rsoil.d <- summaryBy(Rsoil~Date+chamber+T_treatment,data=subset(outdat,chamber %in% c("C04","C05","C06","C07","C08","C09")),
                                                                FUN=mean,keep.names=T)
Rsoil.d$chamber <- factor(Rsoil.d$chamber)

#- daily average Rsoil d13C
isoSR2$Date <- as.Date(isoSR2$Batch.DateTime)
Rsoil.d13C.d <- summaryBy(Keeling_int~Date+chamber+T_treatment,data=subset(isoSR2,chamber %in% c("C04","C05","C06","C07","C08","C09")),,FUN=mean,keep.names=T)

#- merge, add a vector of empty dates
Rsoil <- merge(Rsoil.d,Rsoil.d13C.d,by=c("Date","chamber","T_treatment"))
key <- expand.grid(chamber=levels(Rsoil$chamber),Date=seq.Date(from=as.Date("2016-08-05"),
                                                               to=as.Date("2016-08-16"),by=1))
Rsoil <- merge(key,Rsoil,by=c("chamber","Date"),all.x=T)




#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- loop over each chamber, gapfill missing data, calculate cumulative sum
#- get and plot the cumulative sum of rsoil 13C respiration
Rsoil.l <- split(Rsoil,Rsoil$chamber)
for (i in 1:length(Rsoil.l)){
  #- subset to just dates I want
  Rsoil.l[[i]] <- subset(Rsoil.l[[i]],Date %in% seq.Date(from=as.Date("2016-08-06"),
                                                                   to=as.Date("2016-08-15"),by=1))
  
  Rsoil.l[[i]]$T_treatment <- Rsoil.l[[i]]$T_treatment[1]
  #- gapfill
  Rsoil.l[[i]]$Keeling_int <- na.approx(Rsoil.l[[i]]$Keeling_int)
  Rsoil.l[[i]]$Rsoil <- na.approx(Rsoil.l[[i]]$Rsoil)
  
  
  #- calculate 13C flux rate. Note that the chambers have a diameter of 3.25m
  Rsoil.l[[i]]$AP <- getAP(Rsoil.l[[i]]$Keeling_int) # get the atom percent 13C from per mill data
  Rsoil.l[[i]]$AP_natural <- getAP(-30)         # get the atom percent 13C of "background" respiration
  area <- pi*(3.25/2)^2
  Rsoil.l[[i]]$Rsoil_13C <- with(Rsoil.l[[i]],Rsoil*(AP-AP_natural)/100*area*60*60*24*1e-6*13*1000) #convert to units of mg 13C day-1 EXCESS 13C
  
  
  #- cumulative sum
  Rsoil.l[[i]]$Rsoil_13C_cumsum <- cumsum(Rsoil.l[[i]]$Rsoil_13C)
  
  #------------------------------------------------------------------------
  #------------------------------------------------------------------------
  
  
  
  

}
Rsoil2 <- do.call(rbind,Rsoil.l)

#------------------------------------------------------------------------
#------------------------------------------------------------------------


isoSR <- summaryBy(Rsoil_13C+Rsoil_13C_cumsum~T_treatment+Date,data=Rsoil2,FUN=c(mean,standard.error))


#---
#- plot treatment averages
pdf(file=paste("Output/Rsoil_d13C_sums_",Sys.Date(),".pdf",sep=""))
par(mar=c(6,7,1,4))
plotBy(Rsoil_13C_cumsum.mean~Date|T_treatment,data=isoSR,col=c("blue","red"),pch=16,
       legend="F",axes=F,xlab="",ylab="",ylim=c(0,250),
       panel.first=adderrorbars(x=isoSR$Date,y=isoSR$Rsoil_13C_cumsum.mean,
                                SE=isoSR$Rsoil_13C_cumsum.standard.error,direction="updown"))
magaxis(side=c(2,4),las=1,frame.plot=T)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-08-05"),
                             to=as.Date("2016-08-17"),by="day"),las=2)
title(ylab=expression(R[soil]~(mg~13*C~excess)),main="Cumulative sum, Rsoil",cex.lab=1.5)

#- add a legend
legend("topleft",pch=16,col=c("blue","red"),legend=c("Ambient","Warmed"))
dev.off()

#- sum, then average across treatments
iso.s <- summaryBy(Rsoil_13C~chamber+T_treatment,data=Rsoil2,FUN=c(sum),keep.names=T)
iso.m <- summaryBy(Rsoil_13C~T_treatment,data=iso.s,FUN=mean)

#------------------------------------------------------------------------
#------------------------------------------------------------------------