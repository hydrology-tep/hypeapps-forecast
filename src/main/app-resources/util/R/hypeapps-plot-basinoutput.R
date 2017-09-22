#!/opt/anaconda/envs/cairo-env/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-plot-basinoutput.R
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
# hypeapps-plot-basinoutput.R: Script to plot HYPE basinoutput data, used for TEP Hydrology.
# Author:                      David Gustafsson, SMHI
# Version:                     2017-08-25

# -------------------------------------------------------------------
# 1 - Argument inputs
# -------------------------------------------------------------------
args         = commandArgs(trailingOnly=TRUE)
fileDir      = args[1]   # folder with basinoutput files
plotDir      = args[2]   # folder where the plot files should be written
basinFile    = args[3]   # basinoutput filename
modelName    = args[4]   # modelName (used for plot main title)
hype2csvFile = args[5]   # path to hype2csv file (with geolocations)
prefix.jpg   = args[6]   # csv filename prefix
assimOn      = args[7]   # assim on/off
nObsSim      = args[8]   # number of pairs observation/simulation, to be ploted in same figure
nObsSim      = as.integer(nObsSim)
if(nObsSim>0){
  for(i in 1:nObsSim){
    if(i==1){
      obsSimVars     = args[8+i]   # observation/simulation pairs
    }else{
      obsSimVars     = c(obsSimVars,args[8+i])   # observation/simulation pairs
    }
  }
}
# -------------------------------------------------------------------
# 2 - Working directory
# -------------------------------------------------------------------
setwd(plotDir)

# -------------------------------------------------------------------
# 3 - Dependancies
# -------------------------------------------------------------------
library(rciop)
library(data.table) # to read the basinoutput file
library(Cairo)      # graphics device for output files

pngORjpg = 0        # 0=> png, 1=>jpg

# -------------------------------------------------------------------
# 4- Local function to read HYPE basinoutput file
# -------------------------------------------------------------------
ReadBasinOutput <- function(filename, dt.format = "%Y-%m-%d", type = "df", subid = NULL) {
  
  # handling output type user choice
  if (type == "df") {
    d.t <- F
  } else if (type %in% c("dt", "hmv")) {
    d.t <- T
  } else {
    stop(paste0("Unknown type ", type, "."))
  }
  nm <- strsplit(readLines(filename, n = 1),split = "\t")[[1]]
  x <- fread(filename,  na.strings = c("-9999", "****************"), skip = 2, sep = "\t", header = F, data.table = d.t, 
             colClasses = c("NA", rep("numeric", length(nm) - 1)))      
  names(x) <- c("DATE", nm[-1])
  
  
  ## Date string handling, conditional on import format (HYPE allows for matlab or posix type, without or with hyphens),
  ## handles errors which might occur if the date string differs from the specified format, on error, strings are returned.
  
  
  # if user-requested, hop over date-time conversion
  if (!is.null(dt.format)) {
    
    # convert to posix string if possible, catch failed attempts with error condition and return string unchanged
    # conditional on class of imported data (different syntax for data.table)
    if (is.data.table(x)) {
      
      if (dt.format == "%Y-%m") {
        xd <- as.POSIXct(strptime(paste(x[, DATE], "-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")
        x[, DATE := tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, DATE])})]
      } else if (dt.format == "%Y%m") {
        xd <- as.POSIXct(strptime(paste(x[, DATE], "-01", sep = ""), format = "%Y%m-%d"), tz = "GMT")
        x[, DATE := tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, DATE])})]
      } else if (dt.format == "%Y") {
        xd <- as.POSIXct(strptime(paste(x[, DATE], "-01-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")
        x[, DATE := tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, DATE])})]
      } else {
        xd <- as.POSIXct(strptime(x[, DATE], format = dt.format), tz = "GMT")
        x[, DATE := tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, DATE])})]
      }
      
    } else {
      
      if (dt.format == "%Y-%m") {
        xd <- as.POSIXct(strptime(paste(x[, 1], "-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")
        x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])})
      } else if (dt.format == "%Y%m") {
        xd <- as.POSIXct(strptime(paste(x[, 1], "-01", sep = ""), format = "%Y%m-%d"), tz = "GMT")
        x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])})
      } else if (dt.format == "%Y") {
        xd <- as.POSIXct(strptime(paste(x[, 1], "-01-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")
        x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])})
      } else {
        xd <- as.POSIXct(strptime(x[, 1], format = dt.format), tz = "GMT")
        x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
          print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])})
      }
    }
  } else {
    # dummy date vector as there is always one needed in timestep attribute derivation below
    xd <- NA
  }
  
  ## extract attributes to hold measurement units and SUBID
  munit <- readLines(filename, n = 2)
  munit <- strsplit(munit[2], split = "\t")[[1]][-1]
  # subid conditional on user argument
  if (is.null(subid)) {
    sbd <- strsplit(filename, "/")[[1]]
    sbd <- as.integer(substr(sbd[length(sbd)], start = 1, stop = 7))
    #as.integer(gsub("[[:alpha:][:punct:]]", "", sbd[length(sbd)]))
  } else {
    sbd  <- subid
  }
  # conditional: timestep attribute identified by difference between first two entries
  tdff <- as.numeric(difftime(xd[2], xd[1], units = "hours"))
  if (!is.na(tdff)) {
    if (tdff == 24) {
      tstep <- "day"
    } else if (tdff == 168) {
      tstep <- "week"
    } else if (tdff %in% c(744, 720, 696, 672)) {
      tstep <- "month"
    } else if (tdff %in% c(8760, 8784)) {
      tstep <- "year"
    } else {
      tstep <- paste(tdff, "hour", sep = "")
    }
  } else {
    # add timestep attribute with placeholder value
    tstep <- "none"
  }
  
  # conditional on user choice: output formatting
  if (type %in% c("dt", "df")) {
    
    # update with new attributes
    attr(x, which = "unit") <- munit
    attr(x, which = "subid") <- sbd
    attr(x, which = "timestep") <- tstep
    
    
  } else {
    ## HypeMultiVar formatting
    hvar <- toupper(names(x)[-1])
    # remove dates
    x <- x[, !"DATE", with = F]
    # convert to array (straigtht conversion to array gives error, therefore intermediate matrix)
    x <- as.array(as.matrix(x))
    # adding 'iteration' dimension
    dim(x) <- c(dim(x), 1)
    x <- HypeMultiVar(x = x, date = xd, hype.var = hvar, subid = sbd, tstep = tstep)
  }
  
  return(x)
}

# --------------------------------------------------------------------------
# 5- function to write a world file to geotag an image file
# --------------------------------------------------------------------------
writeWorldFile<-function(fileName, pxWidth, pxHeight, degHeight, lonBasin, latBasin, plotPos="below"){
  # write world file according to definition 
  # http://www.gdal.org/frmt_various.html#WLD
  
  # pixel X size
  pixelSize = degHeight/pxHeight
  writeLinesData=as.character(round(pixelSize,digits = 5)) 
  # rotation about the Y axis (usually 0.0)
  writeLinesData=c(writeLinesData,"0.0")
  # rotation about the X axis (usually 0.0)
  writeLinesData=c(writeLinesData,"0.0")
  # negative pixel Y size
  writeLinesData=c(writeLinesData,as.character(-round(pixelSize,digits = 5)))
  # X coordinate of upper left pixel center
  # Y coordinate of upper left pixel center
  if(plotPos=="below"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin,digits = 5)))
  }else if(plotPos=="upper"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin+degHeight,digits = 5)))
  }else if(plotPos=="center"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin-pixelSize*0.5*pxWidth,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin+degHeight*0.5,digits = 5)))
  }else{
    # default, plotPos=="below"
    writeLinesData=c(writeLinesData,as.character(round(lonBasin,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin,digits = 5)))
  }
  # open file
  fileConn<-file(fileName)
  # write the lines
  writeLines(writeLinesData,con=fileConn,sep="\n")
  # close file
  close(fileConn)
  return(0)  
}

# --------------------------------------------------------------------------
# 5- plot content of the basinoutput as png-file to output dir
# --------------------------------------------------------------------------

# echo arguments to the TEP log file'
rciop.log ("DEBUG", paste(" plot-basinoutput, fileDir: ",fileDir,sep=""), "/util/R/hypeapps-plot-basinoutput.R")
rciop.log ("DEBUG", paste(" plot-basinoutput, plotDir: ",plotDir,sep=""), "/util/R/hypeapps-plot-basinoutput.R")
rciop.log ("DEBUG", paste(" plot-basinoutput, basinFile: ",basinFile,sep=""), "/util/R/hypeapps-plot-basinoutput.R")

# some parameters
graphScale=1.6
lineWidth=2
maxmarg=0.1
lgndScale=0.98

# get subid from basinoutput filename
subid=as.integer(strsplit(basinFile,split = ".txt"))

# read the hype2csv file for geolocations, etc
if(substr(hype2csvFile,1,4)!="9999"){
  hype2csv=read.table(file=hype2csvFile,header=T)
  centerx = hype2csv$CENTERX[which(hype2csv$SUBID==subid)]
  centery = hype2csv$CENTERY[which(hype2csv$SUBID==subid)]
}else{
  centerx = 0
  centery = 0
}
# create plotDir if missing
if(!file.exists(plotDir)) {dir.create(plotDir)}

# Read basin output file
basinData = ReadBasinOutput(paste(fileDir,basinFile,sep="/"))

# Read min/max files if assimOn="on"
if(assimOn=="on"){
  basinDataMin = ReadBasinOutput(paste(fileDir,paste(substr(basinFile,1,nchar(basinFile)-4),"_002.txt",sep=""),sep="/"))
  basinDataMax = ReadBasinOutput(paste(fileDir,paste(substr(basinFile,1,nchar(basinFile)-4),"_003.txt",sep=""),sep="/"))
}

# List variables and their units
variables = colnames(basinData)[2:ncol(basinData)]
units     = attr(basinData,"unit")

rciop.log ("DEBUG", paste(" plot-basinoutput, variables: ",variables,sep=""), "/util/R/hypeapps-plot-basinoutput.R")
rciop.log ("DEBUG", paste(" plot-basinoutput, units: ",units,sep=""), "/util/R/hypeapps-plot-basinoutput.R")

# image filename base
jpegFileBase = paste(plotDir,paste(prefix.jpg,"_",substr(basinFile,1,nchar(basinFile)-4),sep=""),sep="/")

rciop.log ("DEBUG", paste(" plot-basinoutput, FileBase: ",jpegFileBase,sep=""), "/util/R/hypeapps-plot-basinoutput.R")

# main title base
mainBase = paste(substr(basinFile,1,nchar(basinFile)-4)," (",modelName,")",sep="")

# loop over variables
for(i in 1:length(variables)){
  
  # check that there is data in the variable
  if(length(which(!is.na(basinData[,variables[i]])))>0){
  
    # initiate jpeg file for plotting using Cairo graphics device
    if(pngORjpg==1){
      jpegFileName = paste(jpegFileBase,"_",variables[i],".jpg",sep="")
      CairoJPEG(filename = jpegFileName, width = 640, height = 480, units = "px",bg = "white")
    }else{
      jpegFileName = paste(jpegFileBase,"_",variables[i],".png",sep="")
      CairoPNG(filename = jpegFileName, width = 640, height = 480, units = "px",bg = "white")
    }
    # simple plot
    plot(basinData[,1],basinData[,variables[i]],type="l",col="black",
        main = paste(variables[i], "(",units[i],")"," - ",mainBase,sep=""),
        xlab="",ylab=paste(variables[i], "(",units[i],")",sep=""),
        cex.lab=graphScale,cex.axis=graphScale,cex.main=graphScale,
        lwd=lineWidth)
    
    # add min/max
    if(assimOn=="on"){
      lines(basinDataMin[,1],basinDataMin[,variables[i]],lty=2,col="black")
      lines(basinDataMax[,1],basinDataMax[,variables[i]],lty=2,col="black")
    }
    
    # add observation
    if(nObsSim>0){
      for(j in 1:nObsSim){
        obsSimVar1 = tolower(substr(obsSimVars[j],5,8))
        if(obsSimVar1==variables[i]){
          obsSimVar2 = tolower(substr(obsSimVars[j],1,4))
          jObs = which(variables==obsSimVar2)
          if(length(jObs)>0){
            points(basinData[,1],basinData[,jObs+1],pch=18,col="black")
          }
        }
      }
    }

    # close plot device
    dev.off()
  
    # prepare and write worldfile to support the display on the TEP portal map
    if(pngORjpg==1){
      jpegwFileName = paste(jpegFileBase,"_",variables[i],".jpgw",sep="")
    }else{
      jpegwFileName = paste(jpegFileBase,"_",variables[i],".pngw",sep="")
    }
    wfres = writeWorldFile(jpegwFileName, pxWidth=640, pxHeight=480,
                  degHeight=4, lonBasin=centerx, latBasin=centery, plotPos="upper")
  }
}

