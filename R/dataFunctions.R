#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- A collection of functions that return useful bits of data. These functions can 
#  then be called by any script, to keep that script nice and tidy.
#  Naming convention is just "getXXX" where XXX is a short description of the data type.
#  For example, getIso() or getFluxes()
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#---- read in the isotope data, do some manipulation, and return a dataframe
getIso <- function(){
  #-- read in the isotope data
  isodat <- read.csv("Data/13C pulse chase data entry_18Aug.csv")
  names(isodat) <- c("sample.no","Date","Sample.name","Time.collected","Picarro.start.time","Picarro.end.time",
                     "CO2","CO2.sd","d13C","d13C.sd","CH4","CH4.sd","Notes")
  isodat$Date <- as.Date(isodat$Date)
  isodat$Collection.DateTime <- as.POSIXct(paste(isodat$Date,isodat$Time.collected,sep=" "),format="%Y-%m-%d %H:%M",tz="UTC")
  
  
  #- dump variables that are no longer useful
  isodat$sample.no <- isodat$CH4 <- isodat$CH4.sd <- NULL
  
  #- isolate useful information from the sample name
  cc <- base::strsplit(as.character(isodat$Sample.name),split="-")
  isodat$batchtime <- isodat$type <- isodat$chamber <- isodat$timefactor <- NA
  for (i in 1:nrow(isodat)){
    isodat$batchtime[i] <- unlist(cc[i])[1]
    isodat$type[i] <- unlist(cc[i])[2]
    isodat$chamber[i] <- unlist(cc[i])[3]
    isodat$timefactor[i] <- unlist(cc[i])[4]
  }
  
  #- add a time variable for the batch
  isodat$Batch.DateTime <- as.POSIXct(paste(isodat$Date,isodat$batchtime,sep=" "),format="%Y-%m-%d %H",tz="UTC")
  
  #- add the treatment
  linkdf <- data.frame(chamber=c("C01","C02","C04","C05","C06","C07","C08","C09"),
                       T_treatment=c("ambient","elevated",rep(c("elevated","ambient"),3)))
  isodat <- merge(isodat,linkdf,by="chamber")
  
  
  
  
  #--------------------------------------
  #- remove some known bad data
  
  #- leaf respiration on Saturday night (Aug 7), measured at hour 17, were bad 
  #    (bag's weren't closed, as Mike was new to the valves). Remove these data.
  toremove <- which(isodat$Batch.DateTime==as.POSIXct("2016-08-07 17",format="%Y-%m-%d %H",tz="UTC") & isodat$type=="LR")
  isodat2 <- isodat[-toremove,]
  
  #- remove a bad point of soil respiraiton measurement. C05-C at 6am on 6 Aug 2016
  toremove2 <- which(isodat2$Batch.DateTime==as.POSIXct("2016-08-06 06",format="%Y-%m-%d %H",tz="UTC") 
                          & isodat2$type=="SR" & isodat2$chamber=="C05" & isodat2$timefactor=="C")
  isodat3 <- isodat2[-toremove2,]
  
  #--------------------------------------
  
  
  return(isodat3)
}
#-------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------
#---- read in the flux data, do some manipulation, and return a dataframe
getFluxes <- function(){
  file <- list.files(path="Data",pattern="WTCFLUX",full.names=T)
  if (length(file)>1){
    print("error- too many flux datafiles in data folder. Remove all WTCFLUX files except the target file.")
  }
  flux <- as.data.frame(data.table::fread(file))
  
  #- make factors
  flux$chamber <- as.factor(flux$chamber)
  flux$T_treatment <- as.factor(flux$T_treatment)
  
  #- deal with date and time formatting
  flux$DateTime <- as.POSIXct(flux$DateTime,format="%Y-%m-%d %T",tz="UTC")
  flux$Date <- as.Date(flux$DateTime)
  return(flux)
}
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#- function to return the atom percent 13C from a vector of per-mill data
getAP <- function(d13C){
  ARC <- 0.0111803 # absolute ratio of VPDB. See http://www.isotope.uottawa.ca/guides/guides-atom-percent-en.html
  #AP <- (100*ARC*(d13C/1000+1)) / (1+ARC*(d13C/1000+1))
  
  AP <- (d13C/1000+1)*ARC*100 # alternative formulation. Is this right?
  return(AP)
}
#-------------------------------------------------------------------------------------