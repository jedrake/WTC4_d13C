#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- A collection of functions that do useful bits of analysis, such as fit Keeling plots
#  or fit Rsoil data.
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------
#- funtion to fit Rsoil data
fitRsoil <- function(dat,chamber,plotson=F){
  

  
  #- fit non-linear model
  nls1 <- nls(CO2 ~ SSasymp( dTime, Asym, R0, lrc), data = dat)
  preds.nls <- predict(nls1,newdata=data.frame(dTime=seq(0,30,length=101)))
  summary(nls1)
  
  #- plot data and predictions
  if(plotson==T){
    plot(CO2~dTime,data=dat,ylim=c(400,600),xlim=c(0,40))
    lines(preds.nls~seq(1,30,length=101),col="red",lwd=2)
    legend("topright",chamber,bty="n",cex=1.2)
    legend("bottom",legend=as.Date(dat$DateTime[1]),bty="n")
    title(xlab="Time (minutes)",outer=T,cex.lab=2)
    title(ylab="[CO2] (ppm)",outer=T,cex.lab=2)
    
  }
  
  #- calculate derivative based on the first two predictions. units of umol mol-1 min-1
  flux <- preds.nls[6]-preds.nls[5]
  
  vol <- pi*(3.25/2)^2*0.45 # subfloor volume (45 cm tall, 3.25m diameter). Units of m3
  flux2 <- flux/22.4*1000*vol/60 # units of umol CO2 m-2 s-1 (just like chamber Rsoil measures)
  
  #- extract a pseudo r2 value
  preds.forr2 <- predict(nls1,newdata=data.frame(dTime=dat$dTime))
  lm.r2 <- lm(dat$CO2~preds.forr2)
  r2 <- summary(lm.r2)$r.squared
  
  #- return a simple vector with the flux estimate and the pseudo-r2
  return(c(flux2,r2))
}
#-------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------
#- function to fit keeling plots. Accepts and returns a dataframe.
fitKeeling <- function(dat,plotson=F){
  
  #- exctract some bits from the dataframe
  CO2 <- dat$CO2
  d13C <- dat$d13C
  CO2i <- 1/CO2
  
  #- fit the keeling plot, extract the intercept
  lm1 <- lm(d13C~CO2i)
  intercept <- coef(lm1)[1]
  names(intercept) <- "Keeling_int"
  r2 <- summary(lm1)$r.squared
  names(r2) <- "r2"
  
  #- plot, if requested
  if(plotson==T){
    plot(d13C~CO2i)
    abline(lm1)
  }
  

  #- aggregate output
  outputdf <- data.frame(chamber=dat$chamber[1],T_treatment=dat$T_treatment[1],Batch.DateTime=dat$Batch.DateTime[1],
                         Keeling_int=intercept,r2=r2)
  return(outputdf)
}
#-------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
#- function to remove some known-bad flux data
cleanFluxData <- function(flux){

  #- door open
  torm <- which(flux$DoorCnt > 0)
  flux[torm,"FluxCO2"] <- NA
  
  # The central IRGA pump died at ~7pm on 16 Aug. NA-fill those data
  torm1 <- which(flux$DateTime > as.POSIXct("2016-08-16 18:30:00",tz="UTC") & flux$DateTime < as.POSIXct("2016-08-16 21:30:00",tz="UTC"))
  flux[torm1,"FluxCO2"] <- NA
  
  # We also have crazy data after 20:30 on Tuesday 8 Aug continuing to 7am on 9 Aug.
  torm2 <- which(flux$DateTime > as.POSIXct("2016-08-8 20:30:00",tz="UTC") & flux$DateTime < as.POSIXct("2016-08-9 07:15:00",tz="UTC"))
  flux[torm2,"FluxCO2"] <- NA
  
  #- people we in the chambers doing the growth measurements on 17 Aug 2016. The DoorCnt didn't catch all the bad data,
  #  so remove some based on crazy CO2 concentrations
  torm3 <- which(flux$FluxCO2 < -0.03)
  torm4 <- which(flux$FluxCO2 > 0.3)
  flux[c(torm3,torm4),"FluxCO2"] <- NA
  
  #- subtract one day from the date of night-time measurements after midnight but before 5am
  #  This makes it easier to calculate nightly-mean respiration rates
  flux$Date_night <- flux$Date
  toadd <- which(hour(flux$DateTime)<=5)
  flux$Date_night[toadd] <- flux$Date[toadd]-1
  
  return(flux)
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
# Adds error bars to a plot
adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  #if(direction == "updown")
  #  direction <- c("up","down")
  if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("updown" %in% direction) 
    arrows(x0=x, x1=x, y0=y+SE, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#----------------------------------------------------------------------------------------------------------------

