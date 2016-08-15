#------------------------------------------------------------------------
#------------------------------------------------------------------------
#- This script is meant to process the viasala data for soil respiration
#  and return estimates of the soil respiration rate
#------------------------------------------------------------------------
#------------------------------------------------------------------------



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

#- compile the output
outdat <- do.call(rbind,fits.all)
outdat$chamber <- factor(outdat$chamber,levels=c("C01","C02","C04","C05","C06","C07","C08","C09"))
link <- data.frame(chamber=levels(outdat$chamber),
                   T_treatment=c("ambient","elevated","elevated","ambient","elevated","ambient","elevated","ambient"))
outdat <- merge(outdat,link,by="chamber")


#- plot chambers
windows(60,40)
palette(c("green","darkred","red","blue","orange","darkgrey","brown","lightblue"))
plotBy(Rsoil~DateTime|chamber,data=outdat,type="o",pch=16)

#- plot treatments
outdat.m <- summaryBy(Rsoil~T_treatment+DateTime,data=outdat,FUN=c(mean,standard.error))

plotBy(Rsoil.mean~DateTime|T_treatment,data=outdat.m,col=c("blue","red"),pch=16,ylim=c(0,5),
       panel.first=adderrorbars(x=outdat.m$DateTime,y=outdat.m$Rsoil.mean,
                                SE=outdat.m$Rsoil.standard.error,direction="updown"))
