#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-eo-utils.R
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
# hypeapps-eo-utils.R: EO data pre-processing utilities for the HTEP hydrological modelling applications.
# Author:              David Gustafsson, SMHI
# Version:             2018-01-18

# --------------------------------------------------------------------
# Read a csv file from Water Level application and transform to HYPE
# xobs format
# --------------------------------------------------------------------
readWaterLevelCSV<-function(fname=NULL,subid=253,vars="wstr",wl0=NULL,t1=NULL,t2=NULL,modelSHP=NULL,limData=F,lowLim=-0.9999,highLim=9999){
  
  # fname    = path and filename to the water level csv file
  # subid    = subid for which the lake water level data is valid (default 253, which corresponds to lake Kainji)
  # vars     = HYPE model variable name (used as header in XOBS file) - wstr (olake) for instance
  # wl0      = Offset correction (m)
  # t1       = requested start time
  # t2       = requested end time
  # modelSHP = in the future, we can use subbasin polygon data to assign subid
  
  # set offset=0 if missing
  if(is.null(wl0)){
    wl0=0
  }
  
  # read csv file
  if(file.exists(fname)){
    wldata <- read.csv(fname)
  }else{
    return(-1)
  }

  # remove missing data if requested
  if(limData & !is.null(lowLim) & !is.null(highLim)){
    if(highLim>lowLim){
      iData=which(wldata$value>=lowLim & wldata$value<=highLim)
      if(length(iData)>0){
        wldata=wldata[iData,]
      }else{
        return(-1)
      }
    }
  }

  # create time axis
  #
  # Check time format - ouch, the Water level service provide timestep in two formats [YYYY-mm-dd] or [dd-mm-yy]
  thirdChar = substr(as.character(wldata$timestamp[1]),3,3)
  if(thirdChar=="-"){
    timeFormat="dd-mm-yyyy"
  }else{
    timeFormat="yyyy-mm-dd"
  }
  # extract year, month, day
  if(timeFormat=="dd-mm-yyyy"){
    years=substr(as.character(wldata$timestamp),7,10)
    months=substr(as.character(wldata$timestamp),4,5)
    days=substr(as.character(wldata$timestamp),1,2)
  }else{
    years=substr(as.character(wldata$timestamp),1,4)
    months=substr(as.character(wldata$timestamp),6,7)
    days=substr(as.character(wldata$timestamp),9,10)
  }
  data.dateVector=as.POSIXct(paste(years,months,days,sep="-"),tz="GMT")
  startDate=min(data.dateVector)
  endDate=max(data.dateVector)
  if(!is.null(t1)){
    startDate=max(t1,startDate)
  }
  if(!is.null(t2)){
    endDate=max(t2,endDate)
  }
  xobs.dateVector=seq(startDate,endDate,by="days")

  # initiate xobs data frame
  xobs.df=data.frame("DATE"=xobs.dateVector)
  varName=NULL
  subidVec=NULL
  varVec=NULL
  # For now, only read the first variable and first subid
  for(j in 1:1){ #length(vars)){
    for(i in 1:1){ #length(subid)){
      varName=paste(vars[j],as.character(subid[i]),sep="")
      xobs.df[,varName]=NA
      varVec=c(varVec,vars[j])
      subidVec=c(subidVec,subid[i])
    }
  }
  attr(xobs.df,"subid")<-subidVec
  attr(xobs.df,"variable")<-varVec

  # loop over days and assign data to xobs.df
  for(i in 1:nrow(xobs.df)){
    #find data
    jData=which(data.dateVector==xobs.df$DATE[i])
    if(length(jData)>0){
      xobs.df[i,2]=mean(wldata$value[jData])+wl0
    }
  }
  
  # final adjustment, add offset level as attribute as well
  attr(xobs.df,"wl0")<-wl0

  return(xobs.df)
}
