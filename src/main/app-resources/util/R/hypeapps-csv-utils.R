#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-csv-utils.R
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
# hypeapps-csv-utils.R: Utilities to read/write csv time series data in the TEP Hydrology format.
# Author:               Anna Kuenz, David Gustafsson, SMHI
# Version:              2017-08-25

# -------------------------------------------------------------------
# dependancies
# -------------------------------------------------------------------
# source the hypeapps-hype-utils.R

if(app.sys=="tep"){
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-hype-utils.R",sep="/"))
}else if(app.sys=="win"){
  source(paste("application","util/R/hypeapps-hype-utils.R",sep="/"))
}else{
}

## --------------------------------------------------------------------
## Transform HYPE basin output files to TEP csv time series format
## --------------------------------------------------------------------
basinfiles2csv<-function(hypeFile=NULL,csvFile=NULL,hype2csv=NULL,hype2csvFile=NULL,assimOn="off"){
  
  #hype2csvFile="D:/TEP/ModellingService/sandbox/home/dgustafsson/hypeapps-historical/src/main/app-resources/model/niger-hype/shapefiles/niger-hype2csv.txt"
  #hype2csv=NULL
  #hypeFile="D:/TEP/ModellingService/sandbox/home/dgustafsson/hypeapps-historical/src/main/app-resources/modeldata/0000253.txt"
  #csvFile=paste(hypeFile,".csv",sep="")
  
  # List of HYPE variables and units, variables linked to outlet
  # ------------------------------------------------------
  # list of variables and units
  units = c("cout"="m3/s","crun"="m3/s","prec"="mm","temp"="cel","soim"="mm","wcom"="m")
  # list of variables that should be linked to outlet coordinates instead of center of catchment
  var.outlet = c("cout","rout")
  
  # check critical inputs
  # ------------------------------------------------------
  if(is.null(hypeFile)|is.null(csvFile)){
    if(app.sys=="tep"){rciop.log ("DEBUG", " ERROR: NULL in some critical input to writeHype2CSV", "/util/R/hypeapps-csv-utils.R")}
    return(-1)
  }
  # check that hypeFile exists
  if(!file.exists(hypeFile)){
    if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ERROR: writeHype2CSV input hypeFile does not exist =", hypeFile,sep=""), "/util/R/hypeapps-csv-utils.R")}
    return(-1)
  }
  # check that either hype2csv or hype2csvFile exists 
  # read hype2csvFile if hype2csv is missing on input
  # ------------------------------------------------------
  if(is.null(hype2csv)){
    if(is.null(hype2csvFile)){
      if(app.sys=="tep"){rciop.log ("DEBUG", " ERROR: writeHype2CSV input hype2csv and hype2csvFile are both missing, at least one should be provided!", "/util/R/hypeapps-csv-utils.R")}
      return(-1)
    }else{
      if(!file.exists(hype2csvFile)){
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ERROR: writeHype2CSV input file hype2csvFile does not exist = ",hype2csvFile,sep=""), "/util/R/hypeapps-csv-utils.R")}
        return(-1)
      }else{
        # read hype2csvFile
        hype2csv=read.table(file=hype2csvFile,header=T)
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" hype2csv read from file = ",hype2csvFile,sep=""), "/util/R/hypeapps-csv-utils.R")}
      }
    }
  }else{
    if(app.sys=="tep"){rciop.log ("DEBUG", " hype2csv give as function inputfile ", "/util/R/hypeapps-csv-utils.R")}
  }
  
  # read hypeFile
  # ------------------------------------------------------
  hypeData=ReadBasinOutput(filename = hypeFile)
  # also read min and max if assimOn=T
  if(assimOn=="on"){
    hypeDataMin=ReadBasinOutput(filename = paste(substr(hypeFile,1,nchar(hypeFile)-4),"_002.txt",sep=""))
    hypeDataMax=ReadBasinOutput(filename = paste(substr(hypeFile,1,nchar(hypeFile)-4),"_003.txt",sep=""))
  }
  # read subid from filename
  splitFile=strsplit(hypeFile,split="/")[[1]]
  basinFile=splitFile[length(splitFile)]
  basinStr=strsplit(basinFile,split=".txt")[[1]][1]
  subid=as.integer(basinStr)

  if(app.sys=="tep"){rciop.log ("DEBUG", paste(" basinoutput file read: ",hypeFile,sep=""), "/util/R/hypeapps-csv-utils.R")}
  
  # variable and unit names from basinFile attributes
  varStr=colnames(hypeData)
  varStr=varStr[2:length(varStr)]
  unitStr=attr(hypeData,"unit")
  
  # format values to a vector
  # ------------------------------------------------------
  values = unlist(hypeData[,2:ncol(hypeData)])
  if(assimOn=="on"){
    valuesMin = unlist(hypeDataMin[,2:ncol(hypeData)])
    valuesMax = unlist(hypeDataMax[,2:ncol(hypeData)])
  }
  # get coordinates from hype2csv
  # ------------------------------------------------------
  mat = match(subid,hype2csv$SUBID)
  lat=NULL
  long=NULL
  for(i in 1:length(varStr)){
    if(varStr[i] %in% var.outlet){long = c(long,hype2csv$POURX[mat])}else{long = c(long,hype2csv$CENTERX[mat])}
    if(varStr[i] %in% var.outlet){lat = c(lat,hype2csv$POURY[mat])}else{lat = c(lat,hype2csv$CENTERY[mat])}
  }  
  # prepare output table
  # ------------------------------------------------------
  numVars=length(varStr)
  numVals=length(values)
  numDates=nrow(hypeData)
  if(assimOn=="on"){
    csvdata=data.frame(id = seq(1,numVals*3,1),
                       timestamp = c(rep(hypeData$DATE,times=numVars),rep(hypeData$DATE,times=numVars),rep(hypeData$DATE,times=numVars)),
                       longitude = c(rep(long,each=numDates),rep(long,each=numDates),rep(long,each=numDates)),
                       latitude  = c(rep(lat,each=numDates),rep(lat,each=numDates),rep(lat,each=numDates)),
                       uom       = c(rep(unitStr,each=numDates),rep(unitStr,each=numDates),rep(unitStr,each=numDates)),
                       value     = c(values,valuesMin,valuesMax),
                       class     = c(rep(varStr,each=numDates),rep(paste(varStr,"MIN",sep=""),each=numDates),rep(paste(varStr,"MAX",sep=""),each=numDates)),
                       subbasin  = c(rep(subid,times=numVals),rep(subid,times=numVals),rep(subid,times=numVals)))
  }else{
    csvdata=data.frame(id = seq(1,numVals,1),
                   timestamp = rep(hypeData$DATE,times=numVars),
                   longitude = rep(long,each=numDates),
                   latitude  = rep(lat,each=numDates),
                   uom       = rep(unitStr,each=numDates),
                   value     = values,
                   class     = rep(varStr,each=numDates),
                   subbasin  = rep(subid,times=numVals))
  }
  # write output file
  # ------------------------------------------------------
  write.table(csvdata,file=csvFile,sep=",",row.names=FALSE,quote=FALSE)
  if(app.sys=="tep"){rciop.log ("DEBUG", paste(" csv file written: ",csvFile,sep=""), "/util/R/hypeapps-csv-utils.R")}
  return(0)
}

## --------------------------------------------------------------------
## Transform HYPE time output data files to tep-hydro csv format
## --------------------------------------------------------------------
timefiles2csv<-function(hypeFile=NULL,csvFile=NULL,subid=NULL,hype2csv=NULL,hype2csvFile=NULL){
  
  #hype2csvFile="D:/TEP/ModellingService/sandbox/home/dgustafsson/hypeapps-historical/src/main/app-resources/model/niger-hype/shapefiles/niger-hype2csv.txt"
  #hype2csv=NULL
  #hypeFile="D:/TEP/ModellingService/sandbox/home/dgustafsson/hypeapps-historical/src/main/app-resources/modeldata/timeCOUT.txt"
  #csvFile=paste(hypeFile,".csv",sep="")
  #subid=c(37,253)
  
  # List of variables and units, variables linked to outlet
  # ------------------------------------------------------
  # list of variables and units
  units = c("cout"="m3/s","crun"="m3/s","prec"="mm","temp"="cel","soim"="mm","wcom"="m")
  # list of variables that should be linked to outlet coordinates instead of center of catchment
  var.outlet = c("cout","rout")
  
  # check critical inputs
  # ------------------------------------------------------
  if(is.null(hypeFile)|is.null(csvFile)){
    if(app.sys=="tep"){rciop.log ("DEBUG", " ERROR: NULL in some critical input to writeHype2CSV", "/util/R/hypeapps-hype-csv-utils.R")}
    return(-1)
  }
  # check that hypeFile exists
  if(!file.exists(hypeFile)){
    if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ERROR: writeHype2CSV input hypeFile does not exist =", hypeFile,sep=""), "/util/R/hypeapps-hype-csv-utils.R")}
    return(-1)
  }
  # check that either hype2csv or hype2csvFile exists 
  # read hype2csvFile if hype2csv is missing on input
  # ------------------------------------------------------
  if(is.null(hype2csv)){
    if(is.null(hype2csvFile)){
      if(app.sys=="tep"){rciop.log ("DEBUG", " ERROR: writeHype2CSV input hype2csv and hype2csvFile are both missing, at least one should be provided!", "/util/R/hypeapps-hype-csv-utils.R")}
      return(-1)
    }else{
      if(!file.exists(hype2csvFile)){
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ERROR: writeHype2CSV input file hype2csvFile does not exist = ",hype2csvFile,sep=""), "/util/R/hypeapps-hype-csv-utils.R")}
        return(-1)
      }else{
        # read hype2csvFile
        hype2csv=read.table(file=hype2csvFile,header=T)
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" hype2csv read from file = ",hype2csvFile,sep=""), "/util/R/hypeapps-hype-csv-utils.R")}
      }
    }
  }else{
    if(app.sys=="tep"){rciop.log ("DEBUG", " hype2csv give as function inputfile ", "/util/R/hypeapps-csv-utils.R")}
  }
  
  # check requested subid input (if hypeType==timeFile)
  # ------------------------------------------------------
  if(is.null(subid)){writeAllSubid=T}else{writeAllSubid=F}

  # read hypeFile
  # ------------------------------------------------------
  hypeData=ReadTimeOutput(filename = hypeFile)
  # extract subid from data frame attribute
  subidFile = attr(hypeData,"subid")
  # which data columns to include in the csv file?
  if(writeAllSubid){
    timeColOut = 1:length(subidFile)
    subid=subidFile
  }else{
    timeColOut = match(subid,subidFile)
    # check for missing subid in file
    iOk=which(!is.na(timeColOut))
    subid=subid[iOk]
    timeColOut=timeColOut[iOk]
  }

  # variable and unit names
  varStr=tolower(attr(hypeData,"variable"))
  # unit from pre-defined list of units, if variable not found, unit is set to "unknown"
  unitStr=ifelse(!is.na(units[varStr]),units[varStr],"unknown")

  # format values to a vector
  # ------------------------------------------------------
  values = unlist(hypeData[,timeColOut+1])
  
  # get coordinates from hype2csv
  # ------------------------------------------------------
  mat = match(subid,hype2csv$SUBID)
  if(varStr %in% var.outlet){long = hype2csv$POURX[mat]}else{long = hype2csv$CENTERX[mat]}
  if(varStr %in% var.outlet){lat  = hype2csv$POURY[mat]}else{lat  = hype2csv$CENTERY[mat]}

  # prepare output table
  # ------------------------------------------------------
  numBasins=length(subid)
  numDates=nrow(hypeData)
  numVals=length(values)
  csvdata=data.frame(id = seq(1,numVals,1),
                     timestamp = rep(hypeData$DATE,times=numBasins),
                     longitude = rep(long,each=numDates),
                     latitude  = rep(lat,each=numDates),
                     uom       = rep(unitStr,each=numVals),
                     value     = values,
                     class     = rep(varStr,times=numVals),
                     subbasin  = rep(subid,each=numDates))
  
  # write output file
  # ------------------------------------------------------
  write.table(csvdata,file=csvFile,sep=",",row.names=FALSE,quote=FALSE)
  if(app.sys=="tep"){rciop.log ("DEBUG", paste(" csv file written: ",csvFile,sep=""), "/util/R/hypeapps-csv-utils.R")}

  return(0)
}

