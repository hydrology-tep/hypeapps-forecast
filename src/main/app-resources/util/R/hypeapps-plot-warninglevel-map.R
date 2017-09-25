#!/opt/anaconda/envs/cairo-env/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-plot-forecast-map.R
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
# hypeapps-plot-warninglevel-map.R: Script to plot map with forecast warning levels, used for TEP Hydrology.
# Author:                           David Gustafsson, Jafet Andersson, SMHI
# Version:                          2017-09-25
# -------------------------------------------------------------------
# 1 - Argument inputs
# -------------------------------------------------------------------
args         = commandArgs(trailingOnly=TRUE)
fileDir      = args[1]   # folder with timeOutput files (inputs to this script)
plotDir      = args[2]   # folder where the output files should be written (warning levels by text and by plot)
name.hypeout = args[3]   # timeOutput filename (input) 
name.retlev  = args[4]   # return period level file (input, full path) 
name.wl.txt  = args[5]   # warning level text output filename
name.wl.png  = args[6]   # warning level png  output filename
rdataFile    = args[7]   # path to rdata file with subbasin spatial polygons data frame
modelNameIN  = args[8]

# arguments for development
# fileDir      = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test/wl"   # folder with timeOutput files (inputs to this script)
# plotDir      = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test/wl"    # folder where the output files should be written (warning levels by text and by plot)
# name.hypeout = "005_forecast_timeCOUT.txt"   # timeOutput filename (input) 
# name.retlev  = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test/wl/niger-hype-rp-cout.txt"   # return period level file (input, full path) 
# name.wl.txt  = "004_mapWarningLevel.txt"   # warning level text output filename
# name.wl.png  = "004_mapWarningLevel.png"   # warning level text output filename
# modelName    = "niger-hype"   # modelName (used for plot main title)
# rdataFile    = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test/shp.Rdata"   # path to rdata file with subbasin spatial polygons data frame
# modelNameIN  = "niger-hype"

# -------------------------------------------------------------------
# 2 - Dependancies
# -------------------------------------------------------------------
library(rciop)
library(data.table) # to read the basinoutput file
library(Cairo)      # graphics device for output files
library(sp)         # sp for reading Rdata file

# ----------------------------------
# 2.1 ReadTimeOutput from HYPEtools
#--- -------------------------------
ReadTimeOutput <- function(filename, dt.format = "%Y-%m-%d", hype.var = NULL, type = "df", select = NULL, nrows = -1L) {
  
  # argument checks
  if (!is.null(select) && !(1 %in% select)) {
    stop("Argument 'select' must include column 1.")
  }
  
  # handling output type user choice
  if (type == "df") {
    d.t <- F
  } else if (type %in% c("dt", "hsv")) {
    d.t <- T
  } else {
    stop(paste("Unknown type", type, "."))
  }
  
  # import subids, prepare subid attribute vector
  xattr <- readLines(filename, n = 2)
  sbd <- as.numeric(strsplit(xattr[2], split = "\t")[[1]][-1])
  if (!is.null(select)) {
    sbd <- sbd[select[-1] - 1]
  }
  
  # create select vector for fread, workaround for suspected bug in data.table (reported at https://github.com/Rdatatable/data.table/issues/2007)
  if (is.null(select)) {
    select <- 1:(length(sbd) + 1)
  }
  
  #read.table(filename, header = T, na.strings = "-9999", skip = 1)      
  x <- fread(filename,  na.strings = c("-9999", "****************"), skip = 2, sep = "\t", header = F, data.table = d.t, 
             select = select, nrows = nrows)
  #x <- fread(filename,  na.strings = c("-9999", "****************"), skip = 2, sep = "auto", header = F, data.table = d.t, 
  #           select = select, nrows = nrows)
  
  
  # read hype.var from filename, if not provided by user
  if (is.null(hype.var)) {
    hype.var <- substr(strsplit(filename, "time")[[1]][2], start = 1, stop = 4)
  }
  
  # create column names
  names(x) <- c("DATE", paste0("X", sbd))
  
  ## Date string handling, conditional on import format (HYPE allows for matlab or posix type, without or with hyphens),
  ## handles errors which might occur if the date string differs from the specified format. On error, strings are returned.
  
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
  
  # conditional on user choice: output formatting
  if (type %in% c("dt", "df")) {
    
    attr(x, which = "subid") <- sbd
    attr(x, "variable") <- toupper(hype.var)
    
    # conditional: timestep attribute identified by difference between first two entries
    tdff <- as.numeric(difftime(xd[2], xd[1], units = "hours"))
    if (!is.na(tdff)) {
      if (tdff == 24) {
        attr(x, which = "timestep") <- "day"
      } else if (tdff == 168) {
        attr(x, which = "timestep") <- "week"
      } else if (tdff %in% c(744, 720, 696, 672)) {
        attr(x, which = "timestep") <- "month"
      } else if (tdff %in% c(8760, 8784)) {
        attr(x, which = "timestep") <- "year"
      } else {
        attr(x, which = "timestep") <- paste(tdff, "hour", sep = "")
      }
    } else {
      # add timestep attribute with placeholder value
      attr(x, which = "timestep") <- "none"
    }
    
  } else {
    ## HypeSingleVar formatting
    # remove dates
    x <- x[, !"DATE", with = F]
    # convert to array (straigtht conversion to array gives error, therefore intermediate matrix)
    x <- as.array(as.matrix(x))
    # adding 'iteration' dimension
    dim(x) <- c(dim(x), 1)
    x <- HypeSingleVar(x = x, date = xd, subid = sbd, hype.var = toupper(hype.var))
  }
  
  return(x)
}

# --------------------------------
# 2.2 ReadMapOutput from HYPEtools
# --------------------------------
ReadMapOutput <- function(filename, dt.format = NULL, hype.var = NULL, type = "df", nrows = -1L) {
  
  # handling output type user choice
  if (type == "df") {
    d.t <- F
  } else if (type %in% c("dt", "hsv")) {
    d.t <- T
  } else {
    stop(paste("Unknown type", type, "."))
  }
  
  #x <- read.table(filename, header = T, sep = ",", na.strings = "-9999", skip = 1)      
  x <- fread(filename,  na.strings = c("-9999", "****************"), skip = 2, sep = ",", header = F, data.table = d.t, 
             nrows = nrows)
  
  
  # read hype.var from filename, if not provided by user
  if (is.null(hype.var)) {
    hype.var <- substr(strsplit(filename, "map")[[1]][2], start = 1, stop = 4)
  }
  
  # import comment and dates, prepare date attribute vector
  xattr <- readLines(filename, n = 2)
  xd <- strsplit(xattr[2], split = ",")[[1]][-1]
  
  # create column names
  names(x) <- c("SUBID", paste0("X", gsub(pattern = "-", replacement = ".", x = xd)))
  
  ## update with new attributes to hold POSIX dates and timestep keyword, create from column names
  
  
  # if user-requested, hop over date-time conversion
  if (!is.null(dt.format)) {
    # temporary copy to fall back to
    te <- xd
    # convert to posix string if possible, catch failed attempts with error condition and return string unchanged
    if (dt.format == "%Y-%m") {
      xd <- tryCatch(na.fail(as.POSIXct(strptime(paste(xd, "-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")), error = function(e) {
        print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(te)})
    } else if (dt.format == "%Y%m") {
      xd <- tryCatch(na.fail(as.POSIXct(strptime(paste(xd, "-01", sep = ""), format = "%Y%m-%d"), tz = "GMT")), error = function(e) {
        print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(te)})
    } else if (dt.format == "%Y") {
      xd <- tryCatch(na.fail(as.POSIXct(strptime(paste(xd, "-01-01", sep = ""), format = "%Y-%m-%d"), tz = "GMT")), error = function(e) {
        print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(te)})
    } else {
      xd <- tryCatch(na.fail(as.POSIXct(strptime(xd, format = dt.format), tz = "GMT")), error = function(e) {
        print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(te)})
    }
  }
  
  # conditional on user choice: output formatting
  if (type %in% c("dt", "df")) {
    
    attr(x, which = "date") <- xd
    attr(x, "variable") <- toupper(hype.var)
    attr(x, "comment") <- xattr[1]
    
    # conditional: timestep attribute identified by difference between first two entries
    tdff <- tryCatch(as.numeric(difftime(xd[2], xd[1], units = "hours")), error = function(e) {NA})
    if (!is.na(tdff)) {
      if (tdff == 24) {
        attr(x, which = "timestep") <- "day"
      } else if (tdff == 168) {
        attr(x, which = "timestep") <- "week"
      } else if (tdff %in% c(744, 720, 696, 672)) {
        attr(x, which = "timestep") <- "month"
      } else if (tdff %in% c(8760, 8784)) {
        attr(x, which = "timestep") <- "year"
      } else {
        attr(x, which = "timestep") <- paste(tdff, "hour", sep = "")
      }
    } else {
      # add timestep attribute with placeholder value
      attr(x, which = "timestep") <- "none"
    }
    
  } else {
    ## HypeSingleVar formatting
    # copy and remove subids
    sbd <- x[, SUBID]
    x <- x[, !"SUBID", with = F]
    # transpose and convert to array (straigtht conversion to array gives error, therefore intermediate matrix)
    x <- transpose(x)
    x <- as.array(as.matrix(x))
    # adding 'iteration' dimension
    dim(x) <- c(dim(x), 1)
    x <- HypeSingleVar(x = x, date = xd, subid = sbd, hype.var = toupper(hype.var))
  }
  
  return(x)
}

# --------------------------------------------------------------------------
# 2.3 function to write a world file to geotag an image file
# --------------------------------------------------------------------------
writeWorldFile<-function(fileName, pxWidth, pxHeight, degWidth, degHeight, lonBasin, latBasin, plotPos="below"){
  # write world file according to definition 
  # http://www.gdal.org/frmt_various.html#WLD
  
  # pixel X size
  pixelXSize = degWidth/pxWidth
  writeLinesData=as.character(round(pixelXSize,digits = 5)) 
  # rotation about the Y axis (usually 0.0)
  writeLinesData=c(writeLinesData,"0.0")
  # rotation about the X axis (usually 0.0)
  writeLinesData=c(writeLinesData,"0.0")
  # negative pixel Y size
  pixelYSize = degHeight/pxHeight
  writeLinesData=c(writeLinesData,as.character(-round(pixelYSize,digits = 5)))
  # X coordinate of upper left pixel center
  # Y coordinate of upper left pixel center
  if(plotPos=="below"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin,digits = 5)))
  }else if(plotPos=="upper"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin,digits = 5)))
    writeLinesData=c(writeLinesData,as.character(round(latBasin+degHeight,digits = 5)))
  }else if(plotPos=="center"){
    writeLinesData=c(writeLinesData,as.character(round(lonBasin-degWidth*0.5,digits = 5)))
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

# ------------------------------------------------------
# 3 - Derive warning classes and save to a file
# ------------------------------------------------------
{
  # ger return-period magnitudes (return levels)
  rciop.log("DEBUG","Reading return-period magnitudes...", "/util/R/hypeapps-plot-warninglevel-map.R")
  retlev<-read.table(name.retlev,header=T)
  retlev2<-t(retlev[,-1])
  colnames(retlev2)<-retlev[,"SUBID"];rm(retlev)
  # retlev2[,1:5]
  wl.rp<-as.numeric(sub("RP","",rownames(retlev2)))
  if(is.null(warnings())){
    rciop.log("DEBUG","Done.", "/util/R/hypeapps-plot-warninglevel-map.R")
  }else{
    rciop.log("DEBUG","Some warning occurred.", "/util/R/hypeapps-plot-warninglevel-map.R")
  }
  
  # Define function to check wich warning level we should assign
  # Note the principle here is to check for the highest warning level during the entire forecast period (all days)
  # If no return level value is present there is never a warning issued
  wldef<-function(subid) {  
    #subid<-"816"
    myf<-thisq[,subid]
    mywl<-0
    if(!any(is.na(retlev2[,subid]))) {  # ignoring if the return levels could not be estimated, then never warn
      for(k in 1:length(wl.rp)) { #k<-1
        if(any(myf>retlev2[k,subid])) {mywl<-k}
      }}
    return(mywl)
  }
  
  # Derive warning levels for the current forecast and save output
  rciop.log("DEBUG","Deriving warning levels for current forecast...", "/util/R/hypeapps-plot-warninglevel-map.R")

  thisq<-ReadTimeOutput(paste(fileDir,name.hypeout,sep="/"))
  
  # determine first and last date in file
  cdate=as.character(thisq[1,1],format="%Y-%m-%d")
#  edate=as.character(thisq[nrow(thisq),1],format="%Y%m%d")
  
  # continue
  colnames(thisq)<-sub("X","",colnames(thisq))
  thisq<-thisq[,-1]  # remove date column as we don't use it here
  
  if(all(colnames(thisq)==colnames(retlev2))){ # check that the columns match
    thiswl<-sapply(colnames(thisq),FUN = wldef)
    thiswl.df<-data.frame(SUBID=names(thiswl),WarningLevel=thiswl)
    writeLines(text=paste(" Warning levels based on magnitudes with return-period:", paste(wl.rp,collapse=", "), "years"),con=paste(plotDir,name.wl.txt,sep="/"))
    suppressWarnings(write.table(thiswl.df,file=paste(plotDir,name.wl.txt,sep="/"),append=T,row.names=F,quote=F,sep=","))
    rm(thiswl,thiswl.df)
  }else{
    rciop.log("DEBUG","Error: column names in thisq and retlev2 don't match", "/util/R/hypeapps-plot-warninglevel-map.R")
  } 
  
  rm(thisq)  
  
  if(is.null(warnings())){
    rciop.log("DEBUG","Done.", "/util/R/hypeapps-plot-warninglevel-map.R")
  }else{
    rciop.log("DEBUG","Some warning occurred.", "/util/R/hypeapps-plot-warninglevel-map.R")
  }
}

# -----------------------------------------------------
# 4 - Map warning levels
# -----------------------------------------------------
# load subbasin spatial points data frame (shapefileData)
rciop.log("DEBUG","Reading input shapefile...", "/util/R/hypeapps-plot-warninglevel-map.R")
load(rdataFile)

# # for debugging
# shapefileData=shp

if(is.null(warnings())){
  rciop.log("DEBUG","Done.", "/util/R/hypeapps-plot-warninglevel-map.R")
}else{
  rciop.log("DEBUG","Some warning occurred.", "/util/R/hypeapps-plot-warninglevel-map.R")
}

# define colours / graphis
graphScale<-1.6  
wl.alpha=200  # out of 255
wl.col <- rev(heat.colors(length(wl.rp)))
# wl.col=c("yellow","orange","red")  
# plot(1:length(wl.rp),col=wl.col,cex=10,pch=20)

.makeTransparent <- function(someColor, alpha=60) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, maxColorValue = 255)})
}

# get warning levels and define colours to plot
#oldpar<-par()
thiswl<-ReadMapOutput(paste(plotDir,name.wl.txt,sep="/"), hype.var = "WarningLevel")
colnames(thiswl)<-sub("X","",colnames(thiswl))
thiswl[,"col"]<-NA
thiswl[which(thiswl[,"WarningLevel"]==0),"col"]<-NA

for(k in 1:length(wl.col)) { #k<-3
  thiswl[which(thiswl[,"WarningLevel"]==k),"col"]<- .makeTransparent(wl.col[k],wl.alpha)
}
sm<-match(slot(shapefileData,"data")[,"SUBID"],thiswl[,"SUBID"])  #all(thiswl[sm,"SUBID"]==slot(subs,"data")$SUBID)  # match subids

# ---------------------------------
# straightforward PNG plotting
# ---------------------------------
rciop.log("DEBUG","Plotting PNG map...", "/util/R/hypeapps-plot-warninglevel-map.R")

plotFileName = paste(plotDir,name.wl.png,sep="/")

# Plot dimensions and World file
ydiff = bbox(shapefileData)[2,2]-bbox(shapefileData)[2,1]
xdiff = bbox(shapefileData)[1,2]-bbox(shapefileData)[1,1]

width = round(xdiff/(1/60),digits=0)
height = round(width * ydiff/xdiff*1.03,digits=0)


cx = bbox(shapefileData)[1,2] - 0.5 * (bbox(shapefileData)[1,2]-bbox(shapefileData)[1,1])
cy = bbox(shapefileData)[2,2] - 0.5 * (bbox(shapefileData)[2,2]-bbox(shapefileData)[2,1])
wfres = writeWorldFile(paste(plotFileName,"w",sep=""), pxWidth=round(width,digits=0), pxHeight=height,
                       degWidth=xdiff,degHeight=ydiff, lonBasin=cx, latBasin=cy, plotPos="center")

# legend text
modelName=paste(toupper(substr(modelNameIN,1,1)),tolower(substr(modelNameIN,2,nchar(modelNameIN))),sep="")
# hype to HYPE or Hype to HYPE
testHYPE = regexpr(pattern ="hype",modelName)
if(testHYPE[1]>0){
  substr(modelName,testHYPE[1],testHYPE[1]+3)<-"HYPE"
}else{
  testHYPE = regexpr(pattern ="ype",mainTitle)
  if(testHYPE[1]>0){
    substr(modelName,testHYPE[1],testHYPE[1]+2)<-"YPE"
  }
}

# initiate jpeg or png file for plotting using Cairo graphics device
CairoPNG(filename = plotFileName, width = width, height = height, units = "px",bg = "white")

par(xaxs = "i", yaxs = "i", lend = 1,mar=c(0,0,0,0),cex=graphScale)
plot(shapefileData,col=NA,border="grey")  # subbasin boundaries
plot(shapefileData,col=thiswl[sm,"col"],border="NA",add=T)  # warninglevels

legend("topleft",inset=c(0,0.07),title=paste(modelName," 10 day forecast    \n","Issue date ",cdate,"             ",
                                             "\n","River discharge warning levels",sep=""),
       pt.cex=graphScale*4,cex=2,pch=15,bty="n",
       col=.makeTransparent(wl.col,wl.alpha),
       legend=c(paste(" Warning 1   (",as.character(wl.rp[1])," yr RP)",sep=""),
                paste(" Warning 2   (",as.character(wl.rp[2])," yr RP)",sep=""),
                paste(" Warning 3   (",as.character(wl.rp[3])," yr RP)",sep="")))
dev.off()


if(is.null(warnings())){
  rciop.log("DEBUG","Done.", "/util/R/hypeapps-plot-warninglevel-map.R")
}else{
  rciop.log("DEBUG","Some warning occurred.", "/util/R/hypeapps-plot-warninglevel-map.R")
}

