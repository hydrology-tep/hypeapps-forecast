#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-additional-utils.R
#
# Copyright 2016-2017 SMHI
#
# This file is part of H-TEP Hydrological Modelling Application, which is open source 
# and distributed under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version. The Hydrology TEP Hydrological Modelling Application is distributed 
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU 
# General Public License for more details. You should have received a copy of the Lesser 
# GNU General Public License along with the Hydrology TEP Hydrological Modelling Application. 
# If not, see <http://www.gnu.org/licenses/>.
#
# hypeapps-additional-utils.R: Additional supporting utility functions for offline preparation of model data, etc.
# Author:                      David Gustafsson, SMHI
# Version:                     2017-07-13

# some libraries
library(sp)
library(rgdal)
library(rgeos)

# function to make a hype2csv data frame from standard WHIST generated HYPE model shapefiles
shape2csvkey<-function(shpdf){
  csvdf=data.frame("SUBID"=shpdf@data$SUBID)
  csvdf$CENTERX=shpdf@data$CENTERX
  csvdf$CENTERY=shpdf@data$CENTERY
  csvdf$POURX=shpdf@data$POURX
  csvdf$POURY=shpdf@data$POURY
  csvdf$LON=round(coordinates(shpdf)[,1],digits=6)
  csvdf$LAT=round(coordinates(shpdf)[,2],digits=6)
  csvdf$ELEV=shpdf@data$elevation
  return(csvdf)
}
# function to write the hype2csv data frame to table
writeCsvKey<-function(csvdf,filename="hype2csv.txt"){
  write.table(x=csvdf,file = filename,quote = F,sep = "\t",row.names = F,col.names = T)
  return(0)
}

# function transform model shapefile to hype2csv
hypeShapefile2hype2csv<-function(appName="historical",modelName="my-hype",usrName="myuser"){
  dsnStr   = paste("/home/",usrName,"/hypeapps-",appName,"/src/main/app-resources/model/",modelName,"/shapefiles",sep="")
  layerStr = modelName
  # read shapefile
  model.shp = readOGR(dsn=dsnStr,layer=layerStr)

  # make csv key
  csvdf = shape2csvkey(model.shp)
  
  # make [model-name]/hype2csv directory if missing
  
  # write hype2csv.txt
  hype2csvFile = paste("/home/",usrName,"/hypeapps-",appName,"/src/main/app-resources/model/",modelName,"/hype2csv/",modelName,"2csv.txt",sep="")
  writeCsvKey(csvdf = csvdf,filename = hype2csvFile)
  return(0)
}


# function to transform model shapefile to Rdata
hypeShapefile2Rdata<-function(appName="historical",modelName="my-hype",usrName="myuser"){
  dsnStr   = paste("/home/",usrName,"/hypeapps-",appName,"/src/main/app-resources/model/",modelName,"/shapefiles",sep="")
  layerStr = modelName
  # read shapefile
  model.shp = readOGR(dsn=dsnStr,layer=layerStr)
  # save as Rdata
  rdataStr = paste("/home/",usrName,"/hypeapps-",appName,"/src/main/app-resources/model/",modelName,"/shapefiles/",modelName,".shp.Rdata",sep="")
  save(model.shp,file = rdataStr)
  return(0)
}

##
## ReadPTQobs
##
ReadPTQobs <- function (filename, dt.format = "%Y-%m-%d", nrows = -1) {
  
  ## import ptqobs file header, extract attribute
  # import
  xattr <- readLines(filename,n = 1)
  # extract obsids
  sbd <- as.integer(strsplit(xattr, split = "\t")[[1]][-1])
  
  # read the data
  x <- fread(filename,  na.strings = "-9999", sep = "\t", header = T, data.table = F, nrows = nrows)
  #colClasses = c("NA", rep("numeric", length(sbd))))
  
  attr(x, which = "obsid") <- sbd
  
  # date conversion 
  xd <- as.POSIXct(strptime(x[, 1], format = dt.format), tz = "GMT")
  x[, 1] <- xd
  #x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
  #  print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])
  # })
  
  return(x)
}

##
## makeGFDArchive - function to build GFD forcing data archive
##
makeGFDArchive<-function(workDir="",resDir="",useLocal=T){
  
  # download archive from catalogue
  if(!useLocal){
    sysCmd=paste("curl ", " -o ", workdir,"/archive.zip https://store.terradue.com//smhi/gfd/niger-hype/hindcast/files/v1/archive.zip",sep="")
    a=system(sysCmd,intern=T)
  }
  
  # unzip archive to workDir
  archiveFile=paste(workDir,"archive.zip",sep="/")
  if(file.exists(archiveFile)){
    unzip(zipfile = archiveFile,overwrite = T,exdir = workDir)
    archiveFound = T
  }else{
    return(-1)
    archiveFound = F
  }
  
  # load obs-files one by one, and save to separate Rdata files
  if(archiveFound){
    obsFiles = c("Pobs.txt","Tobs.txt","TMAXobs.txt","TMINobs.txt")
    for(i in 1:length(obsFiles)){
      obsFile = paste(workDir,obsFiles[i],sep="/")
      if(file.exists(obsFile)){
        obsData   = ReadPTQobs(obsFile)
        rdataFile = paste(resDir,paste(obsFiles[i],".Rdata",sep=""),sep="/")
        save(obsData,file = rdataFile)
        file.copy(from = obsFile,to = paste(resDir,obsFiles[i],sep="/"),overwrite = T)
      }
    }
    return(0)
  }else{
    return(-1)
  }
}
