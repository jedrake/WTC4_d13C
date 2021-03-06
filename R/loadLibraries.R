#----------------------------------------------------------------------------------------------------------------
#- This script loads all of the required R libraries
#----------------------------------------------------------------------------------------------------------------


# do this once
#devtools::install_bitbucket("remkoduursma/HIEv")
#devtools::install_bitbucket("remkoduursma/plotBy")

library(plotBy)
library(HIEv)
#- little function to check for required libraries, and install them if required
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}

#- load (and install if necessary) packages
Library(dplyr)
Library(stringr)
Library(doBy)
Library(plyr)
Library(RColorBrewer)
Library(reshape2)
Library(data.table)
Library(zoo)
Library(magicaxis)

#- load the data manipulation functions
source("R/dataFunctions.R")
source("R/processingFunctions.R")

#- create directories, if needed
if(!dir.exists("output"))dir.create("Output")
#if(!dir.exists("cache"))dir.create("cache")

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
