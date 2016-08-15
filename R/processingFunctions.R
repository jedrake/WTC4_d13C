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
#- function to fit keeling plots
fitKeeling <- function(CO2,d13C){
  CO2i <- 1/CO2
  lm1 <- lm(d13C~CO2i)
  
  plot(d13C~CO2i)
  abline(lm1)
  
  intercept <- coef(lm1)[1]
  names(intercept) <- "Keeling_int"
  return(intercept)
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

