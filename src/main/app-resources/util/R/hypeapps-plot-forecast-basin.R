#!/opt/anaconda/envs/cairo-env/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-forecast/src/main/app-resources/util/R/hypeapps-plot-forecast-basin.R
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
# hypeapps-plot-forecast-basin.R: Script to plot HYPE forecasts basinoutputs, used for TEP Hydrology.
# Author:                         David Gustafsson, Jafet Andersson SMHI
# Version:                        2017-05-31

# -------------------------------------------------------------------
# 1 - Argument inputs
# -------------------------------------------------------------------
args         = commandArgs(trailingOnly=TRUE)
# at least the first 4 arguments are needed
if(length(args)<4){
  return(-1)
}
for(i in 1:length(args)){
  if(i==1){hindcastDir  = args[i]}             # folder with hindcast output files
  else if(i==2){forecastDir  = args[i]}        # folder with forecast output files
  else if(i==3){plotDir      = args[i]}        # folder where the plot files should be written
  else if(i==4){basinFile    = args[i]}        # basinoutput filename
  else if(i==5){modelName    = args[i]}        # modelName (used for plot main title)
  else if(i==6){hype2csvFile = args[i]}        # path to hype2csv file (with geolocations)
  else if(i==7){name.retlev  = args[i]}        # path to warning levels file
  else if(i==8){prefix.fn    = args[i]}        # filename prefix
#  else if(i==7){longHindFile = args[i]}        # path to longterm hindcast file for regime background
}

# check and flag/replace missing arguments
if(length(args)<5){modelName=""}
if(length(args)<6){
  doWorldFile=F
}else{
  if(file.exists(hype2csvFile)){
    doWorldFile=T
  }else{
    doWorldFile=F
  }
}
# if(length(args)<7){
#   doRegime=F
# }else{
#   if(file.exists(longHindFile)){
#     doRegime=T
#   }else{
#     doRegime=F
#   }
# }
if(length(args)<7){
  doWarningLevels=F
}else{
  if(file.exists(name.retlev)){
    if(length(args)<8){
      doWarningLevels=F
    }else{
      doWarningLevels=T
    }
  }else{
    doWarningLevels=F
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
# 6 - various supporting functions
# --------------------------------------------------------------------------
## makeTransparent plotting color
.makeTransparent <- function(someColor, alpha=60) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, maxColorValue = 255)})
}
## function to make character vector with days from a posix vector
dates2days<-function(dates){
  days=as.character(dates)
  days=substr(days,6,10)
  return(days)
}
# --------------------------------------------------------------------------
# 7 - function to plot streamflow hindcast, forecast, observation, with 
#     warning levels and simulated regime in the background
# --------------------------------------------------------------------------
plotHindForWL<-function(hindcast=NULL,
                        forecast=NULL,
#                        obsdata=NULL,
#                        annsim=NULL,
#                        annobs=NULL,
                        wl.sim=NULL,
#                        wl.obs=NULL,
                        wl.col=c("yellow","orange","red"),
#                        maxsim=NULL,
#                        maxobs=NULL,
                        maintitle=NULL,
                        hinddays=NULL,
                        maxmarg=0.10,
                        wl.alpha=120,
#                        wl.alpha.obs=120,
                        hind.col="blue",
                        fore.col="red",
#                        obs.col="black",
                        lineWidth=3,
                        graphScale=2,
                        lgndScale=1,
                        lgndDY=-0.2){
  
  ## initiate plot
  if(!is.null(hinddays)){
    nh=nrow(hindcast)
    dateAxis=c(hindcast$DATE[(nh-hinddays+1):nh],forecast$DATE)
    iniData=c(hindcast[(nh-hinddays+1):nh,2],forecast[,2])
  }else{
    hinddays=nrow(hindcast)
    nh=nrow(hindcast)
    dateAxis=c(hindcast$DATE,forecast$DATE)
    iniData=c(hindcast[,2],forecast[,2])
  }
#  maxy=max(c(maxobs,maxsim,max(hindcast[,2]),max(forecast[,2]),max(obsdata[,2],na.rm = T)),na.rm = T)*(1+maxmarg)
#  miny=0
  maxy=max(c(max(hindcast[,2]),max(forecast[,2]),max(wl.sim,na.rm = T)),na.rm = T)*(1+maxmarg)
  miny=0
  
  plot(dateAxis,iniData,type="n",ylim=c(miny,maxy),main=maintitle,xlab="",ylab="Discharge (m3/s)",xaxs="i",cex.lab=graphScale,cex.axis=graphScale,cex.main=graphScale,xpd=T)
  axis(side = 1,at=dateAxis-60*60*2,labels=F,tcl=-0.2)  # add daily tickmarks for simplicity
  
  ## warning levels as polygons with transparent colors
  if(!is.null(wl.sim)){
    polcol <- .makeTransparent(wl.col[1], wl.alpha)
    if(length(wl.col)>1){
      for(j in 2:length(wl.col)){
        polcol <- c(polcol,.makeTransparent(wl.col[j], wl.alpha))
      }
    }
    for(i in 1:length(wl.sim)){
      x=c(dateAxis[1],dateAxis[length(dateAxis)],dateAxis[length(dateAxis)],dateAxis[1])
      y=c(wl.sim[i],wl.sim[i])
      if(i==length(wl.sim)){
        y=c(y,maxy*2,maxy*2)
      }else{
        y=c(y,wl.sim[i+1],wl.sim[i+1])
      }
      polygon(x, y, col = polcol[i], border = NA)
    }
  }

#  ## warning levels from observations
#  if(!is.null(wl.obs)){
#    polcol <- .makeTransparent(wl.col[1], wl.alpha.obs)
#    if(length(wl.col)>1){
#      for(j in 2:length(wl.col)){
#        polcol <- c(polcol,.makeTransparent(wl.col[j], wl.alpha))
#      }
#    }
#    for(i in 1:length(wl.obs)){
#      lines(c(min(dateAxis),max(dateAxis)),c(1,1)*wl.obs[i],col=polcol[i],lty=1,lwd=lineWidth)
#    }
#  }    
  
#  ## maximum ever observed and simulated
#  if(!is.null(maxsim)){
#    lines(c(min(dateAxis),max(dateAxis)),c(1,1)*maxsim,col=hind.col,lty=2,lwd=lineWidth)
#  }
#  if(!is.null(maxobs)){
#    lines(c(min(dateAxis),max(dateAxis)),c(1,1)*maxobs,col=obs.col,lty=2,lwd=lineWidth)
#  }
  
#  ## regime ribbons
#  days2plot = dates2days(dateAxis)
#  ind2plot  = match(days2plot,annsim$mean$day)
#  
#  # simulation regime
#  polcol <- .makeTransparent(hind.col, 30)
#  mindata=cbind(dateAxis,annsim$minimum[ind2plot, 3])
#  maxdata=cbind(dateAxis,annsim$maximum[ind2plot,3])
#  polcoor.minmax <- rbind(mindata, maxdata[nrow(maxdata):1,])
#  mindata=cbind(dateAxis,annsim$p25[ind2plot, 3])
#  maxdata=cbind(dateAxis,annsim$p75[ind2plot,3])
#  polcoor.p25p75 <- rbind(mindata, maxdata[nrow(maxdata):1,])
#  polygon(polcoor.minmax[, 1], polcoor.minmax[, 2], col = polcol, border = NA)

  # observation regime
  #  polcol <- .makeTransparent("black", 50)
  #  mindata=cbind(dateAxis,annobs$minimum[ind2plot, 3])
  #  maxdata=cbind(dateAxis,annobs$maximum[ind2plot,3])
  #  polcoor.minmax <- rbind(mindata, maxdata[nrow(maxdata):1,])
  #  mindata=cbind(dateAxis,annobs$p25[ind2plot, 3])
  #  maxdata=cbind(dateAxis,annobs$p75[ind2plot,3])
  #  polcoor.p25p75 <- rbind(mindata, maxdata[nrow(maxdata):1,])
  #  polygon(polcoor.minmax[, 1], polcoor.minmax[, 2], col = polcol, border = NA)
  #  polygon(polcoor.p25p75[, 1], polcoor.p25p75[, 2], col = polcol, border = NA)
  
  # hindcast
  lines(rbind(hindcast[(nh-hinddays+1):nh,],forecast[1,]),col=.makeTransparent(hind.col, 200),lwd=lineWidth,lty=1)
  # forecast
  lines(forecast,col=.makeTransparent(fore.col, 200),lwd=lineWidth,lty=1,)
#  # observation
#  if(!is.null(obsdata)){
#    lines(obsdata,col=.makeTransparent(obs.col, 200),lwd=lineWidth,lty=1)
#    #  points(obsdata,col=.makeTransparent("black", 200),pch=20)
#  }

#  # Legend 1: forecast, hindcast, maxima and observations (lines)
#  legend(x=par("usr")[2],y=par("usr")[3],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth, yjust=0,pt.cex=graphScale*2,
#         legend = c("10 day forecast","Current hindcast","Hindcast max (1979-2015)","Hindcast daily range (1979-2015)",if(!is.null(obsdata)) {"Observed"},if(!is.null(maxobs)) {"Observed max (1979-2015)"}), 
#         col = c(.makeTransparent(fore.col, alpha = 200),.makeTransparent(hind.col, alpha = 200),hind.col,.makeTransparent(hind.col, 30),.makeTransparent(obs.col, 200),obs.col),
#         lty = c(1,1,2,NA,1,2),
#         pch = c(NA,NA,NA,15,NA,NA))
  
  # Legend 1: forecast, hindcast (lines)
  legend(x=par("usr")[2],y=par("usr")[3],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth, yjust=0,pt.cex=graphScale*2,
         legend = c("10 day forecast","Current hindcast"), 
         col = c(.makeTransparent(fore.col, alpha = 200),.makeTransparent(hind.col, alpha = 200)),
         lty = c(1,1),
         pch = c(NA,NA))
  
#  # Legend 2: warning levels and hindcast range for this period (polygons)
#  legend(x=par("usr")[2],y=par("usr")[4],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth,lty=NA,pch=15, pt.cex=graphScale*2,
#         legend = c(paste("Warning 3 (",wl.rp[3]," yr RP)",sep=""),paste("Warning 2 (",wl.rp[2]," yr RP)",sep=""),paste("Warning 1 (",wl.rp[1]," yr RP)",sep="")),
#         col=c(.makeTransparent(wl.col[3], wl.alpha),.makeTransparent(wl.col[2], wl.alpha),.makeTransparent(wl.col[1], wl.alpha)))

  # Legend 2: warning levels (polygons)
  legend(x=par("usr")[2],y=par("usr")[4],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth,lty=NA,pch=15, pt.cex=graphScale*2,
         legend = c(paste("Warning 3 (",wl.rp[3]," yr RP)",sep=""),paste("Warning 2 (",wl.rp[2]," yr RP)",sep=""),paste("Warning 1 (",wl.rp[1]," yr RP)",sep="")),
         col=c(.makeTransparent(wl.col[3], wl.alpha),.makeTransparent(wl.col[2], wl.alpha),.makeTransparent(wl.col[1], wl.alpha)))
  
}

# --------------------------------------------------------------------------
# SCRIPT: plot content of the basinoutput as jpg-file to output dir
# --------------------------------------------------------------------------

# echo arguments to the TEP log file'
rciop.log ("DEBUG", paste(" plot-forecast-basin, 1 hindcastDir: ",hindcastDir,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 2 forecastDir: ",forecastDir,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 3 plotDir: ",plotDir,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 4 basinFile: ",basinFile,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 5 modelName: ",modelName,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 6 hype2csvFile: ",hype2csvFile,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 7 name.retlev: ",name.retlev,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")
rciop.log ("DEBUG", paste(" plot-forecast-basin, 8 prefix.fn: ",prefix.fn,sep=""), "/util/R/hypeapps-plot-forecast-basin.R")

# get subid from basinoutput filename
subid=as.integer(strsplit(basinFile,split = ".txt"))

# read the hype2csv file for geolocations, etc
if(doWorldFile){
  hype2csv=read.table(file=hype2csvFile,header=T)
  centerx = hype2csv$CENTERX[which(hype2csv$SUBID==subid)]
  centery = hype2csv$CENTERY[which(hype2csv$SUBID==subid)]
}

# Read the longterm hindcast data for background regime plot (time file)
#if(doRegime){
#  longHindData=ReadTimeOutput(longHindFile)
#}

# Read the return level file
retlev<-read.table(name.retlev,header=T)

# return periods in file
rp.cols = colnames(retlev)[2:ncol(retlev)]
rp.file = NULL
for(i in 1:length(rp.cols)){
  rp.file = c(rp.file,as.integer(substr(rp.cols[i],3,nchar(rp.cols[i]))))
}

# (first three) return levels for this subid
iretlev = which(retlev$SUBID==subid)
wl.rp=c(0,0,0)*NA
wl.sim=c(0,0,0)
for(i in 1:min(3,length(rp.file))){
  wl.rp[i]=rp.file[i]
  wl.sim[i]=retlev[iretlev,i+1]
}
#for(i in 1:nwl){
#  jrp=which(rp.file==wl.rp[i])
#  wl.sim[i]=retlev[iretlev,jrp+1]
#}

# Read basin output files (hindcast and forecast)
hindcastData = ReadBasinOutput(paste(hindcastDir,basinFile,sep="/"))
forecastData = ReadBasinOutput(paste(forecastDir,basinFile,sep="/"))

# Get issueDate from forecast time series
issueDate = forecastData$DATE[1]

# create plotDir if missing
if(!file.exists(plotDir)) {dir.create(plotDir)}

pngORjpg = 0
if(pngORjpg==1){
  fileExt = "jpg"
}else{
  fileExt = "png"
}

# image filename
plotFileName = paste(plotDir,paste(prefix.fn,substr(basinFile,1,nchar(basinFile)-4),sep="_"),sep="/")
plotFileName = paste(plotFileName,"_","discharge-forecast",".",fileExt,sep="")

if(doWorldFile){
  # prepare and write worldfile to support the display on the TEP portal map
  plotwFileName = paste(plotFileName,"w",sep="")
  wfres = writeWorldFile(plotwFileName, pxWidth=640, pxHeight=480,
                         degHeight=4, lonBasin=centerx, latBasin=centery, plotPos="upper")
}

# open file for plotting with Cairo device
if(pngORjpg==1){
  CairoJPEG(filename = plotFileName, width = 640, height = 480, units = "px",bg = "white")
}else{
  CairoPNG(filename = plotFileName, width = 640, height = 480, units = "px",bg = "white")
}
# Modify plot parameters
oldpar<-par()
par(mgp=c(2.5,0.8,0),mar=c(3.1,4.1,3.1,23))

# Make the plot
plotHindForWL(hindcast=hindcastData[,c("DATE","cout")],
              forecast=forecastData[,c("DATE","cout")],
#              obsdata= NULL,
#              annsim=NULL,
#              annobs = NULL,
              wl.sim=wl.sim,
#              wl.obs=NULL,
              wl.col=c("yellow","orange","red"),
#              maxsim=maxsim,
#              maxobs=NULL,
              maintitle=paste(paste("Sub-basin ",as.character(subid),sep=""),paste("issue date ",as.character(issueDate),sep=""),sep=" - "),
              hinddays=50,
              maxmarg=0.1,
              wl.alpha=90,
#              wl.alpha.obs=240,
              hind.col="blue",
              fore.col="red",
 #             obs.col="black",
              lineWidth=4,
              graphScale=1.6,
              lgndScale = 0.98)

# Reset plot parameters
suppressWarnings(par(oldpar)) # warns for some parameters that cannot be modified... 

# close plotting device
dev.off()









