#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-utils.R
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
# hypeapps-utils.R: R tools for the HTEP hydrological modelling application 
# Author:           David Gustafsson, SMHI
# Version:          2018-02-07 
#

## --------------------------------------------------------------------------------
## initial settings
## --------------------------------------------------------------------------------
# set system flag if not set
if(!exists("app.sys")){
  app.sys ="tep"
}

# load rciop library
if(app.sys=="tep"){
  library("rciop")
}

# source hypeapps environment file, if needed
if(!exists("app.envset")){
  if(app.sys=="tep"){
    source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-environment.R",sep="/"))
  }else if(app.sys=="win"){
    source("application/util/R/hypeapps-environment.R")  
  }
  if(app.sys=="tep"){rciop.log ("DEBUG", paste("hypeapps-environment.R sourced"), "/util/R/hypeapps-utils.R")}
}

# source csv and hype-utils.R
if(app.sys=="tep"){
  
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-hype-utils.R",sep="/"))
  rciop.log ("DEBUG", paste("hypeapps-hype-utils.R sourced"), "/util/R/hypeapps-utils.R")
  
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-csv-utils.R",sep="/"))
  rciop.log ("DEBUG", paste("hypeapps-csv-utils.R sourced"), "/util/R/hypeapps-utils.R")
  
  if(app.name=="eodata"){
    source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-eo-utils.R",sep="/"))
    rciop.log ("DEBUG", paste("hypeapps-eo-utils.R sourced"), "/util/R/hypeapps-utils.R")
  }
  if(app.name=="returnperiod"){
    source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-returnperiod-utils.R",sep="/"))
    rciop.log ("DEBUG", paste("hypeapps-returnperiod-utils.R sourced"), "/util/R/hypeapps-utils.R")
  }
}else if(app.sys=="win"){
  source(paste("application","util/R/hypeapps-hype-utils.R",sep="/"))
  source(paste("application","util/R/hypeapps-csv-utils.R",sep="/"))
  if(app.name=="eodata"){source(paste("application","util/R/hypeapps-eo-utils.R",sep="/"))}
  if(app.name=="returnperiod"){source(paste("application","util/R/hypeapps-returnperiod-utils.R",sep="/"))}
}else{
}

## --------------------------------------------------------------------------------
## functions
## -------------------------------------------------------------------------------
## getHypeAppInput - function to load user parameter inputs, depending on application
getHypeAppInput<-function(appName){
  
  ## HISTORICAL ##  
  if(appName=="historical"){
    if(app.sys=="tep"){
      # get parameters with the rciop function when running on the TEP system
      cdate       <- rciop.getparam("cdate")       # start of results output
      edate       <- rciop.getparam("edate")       # end of simulation
      outvarsIN   <- rciop.getparam("variables")   # output variables
      basinselect <- rciop.getparam("basinselect") # output subbasins
      basinset    <- rciop.getparam("basinset")    # output subbasins
      xobs        <- rciop.getparam("xobs")        # EO/Insitu data Xobs file(s)
      assimOn     <- rciop.getparam("assimOn")     # Assimilation on/off
      assimVarIN  <- rciop.getparam("assimVars")   # Assimilation variables
      
      # parse the outvar inputs
      outvarSplit = trimws(strsplit(outvarsIN,split = ",")[[1]])
      nOut=length(outvarSplit)
      for(i in 1:nOut){
        splitFirst = trimws(strsplit(outvarSplit[i],split = "(",fixed=T)[[1]])
        oName=splitFirst[1]
        splitSecond = trimws(strsplit(splitFirst[2],split = ")",fixed=T)[[1]])
        oUnit=splitSecond[1]
        oID=substr(splitSecond[2],nchar(splitSecond[2])-3,nchar(splitSecond[2]))
        oName2=substr(splitSecond[2],1,nchar(splitSecond[2])-5)
        if(i==1){
          outvars.name=oName
          outvars.name2=oName2
          outvars.unit=oUnit
          outvars.id=oID
        }else{
          outvars.name=c(outvars.name,oName)
          outvars.name2=c(outvars.name2,oName2)
          outvars.unit=c(outvars.unit,oUnit)
          outvars.id=c(outvars.id,oID)
        }
      }
      outvars.num=nOut
      outvars=outvars.id[1]
      if(nOut>1){
        for(i in 2:nOut){
          outvars=paste(outvars,outvars.id[i],sep=",")
        }
      }
      # parse the basinselect and basinset intputs
      basinSplit = trimws(strsplit(basinselect,split = ",")[[1]])
      nBasin=length(basinSplit)
      for(i in 1:nBasin){
        basinSplit2=trimws(strsplit(basinSplit[i],split = " ")[[1]])
        bID=basinSplit2[length(basinSplit2)]
        bName=substr(basinSplit[i],1,nchar(basinSplit[i])-nchar(bID)-1)
        if(i==1){
          basins.name=bName
          basins.id=bID
        }else{
          basins.name=c(basins.name,bName)
          basins.id=c(basins.id,bID)
        }
      }
      basins.num=nBasin
      basins=basins.id[1]
      if(nBasin>1){
        for(i in 1:(nBasin-1)){
          basins=paste(basins,basins.id[i+1],sep=",")   
        }
      }
      if(basinset!="-9999" & nchar(basinset)>0){
        basinSplit=trimws(strsplit(basinset,split = ",")[[1]])
        # check if any of the basins in basinset was already selected by the basinselect
        iBasin=match(basinSplit,basins.id)
        iAdd=which(is.na(iBasin))
        if(length(iAdd)>0){
          basinAdd=basinSplit[iAdd]
          for(i in 1:length(iAdd)){
            basins=paste(basins,basinAdd[i],sep=",")
            basins.name=c(basins.name,"NN")
            basins.num=basins.num+1
            basins.id=c(basins.id,basinAdd[i])
          }      
        }
      }
      # If xobs !=-9999, parse the input to URLs
      if(length(xobs)>0){
        xobsNum=0
        xobsURL=NULL
        for(i in 1:length(xobs)){
          if(nchar(xobs[i])>1 & xobs[i]!="-9999"){
            # update the number of xobs inputs
            xobsNum=xobsNum+1
            
            # extract the URL from xobs input string:
            if(xobsNum==1){
              xobsURL=strsplit(xobs[i],split = "&")[[1]][1]
            }else{
              xobsURL=c(xobsURL,strsplit(xobs[i],split = "&")[[1]][1])
            }
          }
        }
      }else{
        xobsNum=0
        xobsURL=NULL
      }
      
      # parse the AssimVarIn input
      if(assimOn=="on"){
        # <option>Lake Water Level - altimetry AOWL WCOM</option>
        assimVarSplit = trimws(strsplit(assimVarIN,split = ",")[[1]])
        nOut=length(assimVarSplit)
        for(i in 1:nOut){
          aName=substr(assimVarSplit[i],1,nchar(assimVarSplit[i])-10)
          asID=substr(assimVarSplit[i],nchar(assimVarSplit[i])-3,nchar(assimVarSplit[i]))
          aoID=substr(assimVarSplit[i],nchar(assimVarSplit[i])-8,nchar(assimVarSplit[i])-5)
          if(i==1){
            assimVar.name=aName
            assimVar.simid=asID
            assimVar.obsid=aoID
          }else{
            assimVar.name=c(assimVar.name,aName)
            assimVar.simid=c(assimVar.simid,asID)
            assimVar.obsid=c(assimVar.obsid,aoID)
          }
        }
        assimVar.num=nOut
        assimVar = paste(assimVar.obsid[1],assimVar.simid[1],sep=",")
        if(nOut>1){
          for(i in 1:(nOut-1)){
            assimVar=c(assimVar,paste(assimVar.obsid[i+1],assimVar.simid[i+1],sep=","))
          }
        }
      }else{
        assimVar="9999,9999"
      }
      
    }else if(app.sys=="win"){
      # set to default values, if not set
      # set some test values for development on windows:
      if(#!exists("win.bdate")|
        !exists("win.cdate")|
        !exists("win.edate")|
        !exists("win.outvars")|
        !exists("win.outbasins")|
        !exists("win.xobs")){
        #  bdate     <- "2006-01-01"
        cdate     <- "2016-01-01"
        edate     <- "2016-12-01"
        outvars   <- "cout,cprc,ctmp,evap,epot"
        outbasins <- "37"
        xobs <- "-9999"
        xobsNum=0
        xobsURL=NULL
        assimOn="off"
        assimVar="9999,9999"
      }else{
        #bdate     <- win.bdate
        cdate     <- win.cdate
        edate     <- win.edate
        outvars   <- win.outvars
        outbasins <- win.outbasins
        outbasins <- win.xobs
        xobs <- "-9999"
        xobsNum=0
        xobsURL=NULL
        assimOn="off"
        assimVar="9999,9999"
      }
    }else{
      #bdate     <- NULL
      cdate     <- NULL
      edate     <- NULL
      outvars   <- NULL
      outbasins <- NULL
      xobs <- NULL
      xobsNum=0
      xobsURL=NULL
      assimOn="off"
      assimVar="9999,9999"
      print("WARNING: hypeapps.sys not set, allowed values are 'tep' or 'win' ")
    }
    
    # Return Input List
    appInput=list(#"bdate"     = bdate,
      "cdate"     = cdate,
      "edate"     = edate,
      "outvars"   = outvars,
      "outvars.name" = outvars.name,
      "outvars.name2" = outvars.name2,
      "outvars.unit" = outvars.unit,
      "outvars.id"   = outvars.id,
      "outvars.num"  = outvars.num,
      "outbasins" = basins,
      "outbasins.name"=basins.name,
      "outbasins.id"=basins.id,
      "outbasins.num"=basins.num,
      "xobs"      = xobs,
      "xobsNum"   = xobsNum,
      "xobsURL"   = xobsURL,
      "assimOn"   = assimOn,
      "assimVar"  = assimVar)
    
    ## FORECAST ##
  }else if(appName=="forecast"){
    
    if(app.sys=="tep"){
      # get parameters with rciop function when running on the TEP system
      idate       <- rciop.getparam("idate")       # Forecast issue date
      outvarsIN   <- rciop.getparam("variables")   # output variables
      basinselect <- rciop.getparam("basinselect") # output subbasins
      basinset    <- rciop.getparam("basinset")    # output subbasins
      rpcout      <- rciop.getparam("rpcout")      # Return periods levels file
      xobs        <- rciop.getparam("xobs")        # EO/Insitu data Xobs file
      
      # temporarily commenting out the assimilation in the forecast application (David 20170827)
      assimOn     <- rciop.getparam("assimOn")     # Assimilation on/off
      assimVarIN  <- rciop.getparam("assimVars")   # Assimilation variables
      # parse the outvar inputs
      outvarSplit = trimws(strsplit(outvarsIN,split = ",")[[1]])
      nOut=length(outvarSplit)
      for(i in 1:nOut){
        splitFirst = trimws(strsplit(outvarSplit[i],split = "(",fixed=T)[[1]])
        oName=splitFirst[1]
        splitSecond = trimws(strsplit(splitFirst[2],split = ")",fixed=T)[[1]])
        oUnit=splitSecond[1]
        oID=substr(splitSecond[2],nchar(splitSecond[2])-3,nchar(splitSecond[2]))
        oName2=substr(splitSecond[2],1,nchar(splitSecond[2])-5)
        if(i==1){
          outvars.name=oName
          outvars.name2=oName2
          outvars.unit=oUnit
          outvars.id=oID
        }else{
          outvars.name=c(outvars.name,oName)
          outvars.name2=c(outvars.name2,oName2)
          outvars.unit=c(outvars.unit,oUnit)
          outvars.id=c(outvars.id,oID)
        }
      }
      outvars.num=nOut
      outvars=outvars.id[1]
      if(nOut>1){
        for(i in 2:nOut){
          outvars=paste(outvars,outvars.id[i],sep=",")
        }
      }
      # parse the basinselect and basinset intputs
      basinSplit = trimws(strsplit(basinselect,split = ",")[[1]])
      nBasin=length(basinSplit)
      for(i in 1:nBasin){
        basinSplit2=trimws(strsplit(basinSplit[i],split = " ")[[1]])
        bID=basinSplit2[length(basinSplit2)]
        bName=substr(basinSplit[i],1,nchar(basinSplit[i])-nchar(bID)-1)
        if(i==1){
          basins.name=bName
          basins.id=bID
        }else{
          basins.name=c(basins.name,bName)
          basins.id=c(basins.id,bID)
        }
      }
      basins.num=nBasin
      basins=basins.id[1]
      if(nBasin>1){
        for(i in 1:(nBasin-1)){
          basins=paste(basins,basins.id[i+1],sep=",")   
        }
      }
      if(basinset!="-9999" & nchar(basinset)>0){
        basinSplit=trimws(strsplit(basinset,split = ",")[[1]])
        # check if any of the basins in basinset was already selected by the basinselect
        iBasin=match(basinSplit,basins.id)
        iAdd=which(is.na(iBasin))
        if(length(iAdd)>0){
          basinAdd=basinSplit[iAdd]
          for(i in 1:length(iAdd)){
            basins=paste(basins,basinAdd[i],sep=",")
            basins.name=c(basins.name,"NN")
            basins.num=basins.num+1
            basins.id=c(basins.id,basinAdd[i])
          }      
        }
      }
      # If xobs !=-9999, parse the input to URLs
      if(length(xobs)>0){
        xobsNum=0
        xobsURL=NULL
        for(i in 1:length(xobs)){
          if(nchar(xobs[i])>1 & xobs[i]!="-9999"){
            # update the number of xobs inputs
            xobsNum=xobsNum+1
            
            # extract the URL from xobs input string:
            if(xobsNum==1){
              xobsURL=strsplit(xobs[i],split = "&")[[1]][1]
            }else{
              xobsURL=c(xobsURL,strsplit(xobs[i],split = "&")[[1]][1])
            }
          }
        }
      }else{
        xobsNum=0
        xobsURL=NULL
      }
      
      # parse the AssimVarIn input
      if(assimOn=="on"){
        # <option>Lake Water Level - altimetry AOWL WCOM</option>
        assimVarSplit = trimws(strsplit(assimVarIN,split = ",")[[1]])
        nOut=length(assimVarSplit)
        for(i in 1:nOut){
          aName=substr(assimVarSplit[i],1,nchar(assimVarSplit[i])-10)
          asID=substr(assimVarSplit[i],nchar(assimVarSplit[i])-3,nchar(assimVarSplit[i]))
          aoID=substr(assimVarSplit[i],nchar(assimVarSplit[i])-8,nchar(assimVarSplit[i])-5)
          if(i==1){
            assimVar.name=aName
            assimVar.simid=asID
            assimVar.obsid=aoID
          }else{
            assimVar.name=c(assimVar.name,aName)
            assimVar.simid=c(assimVar.simid,asID)
            assimVar.obsid=c(assimVar.obsid,aoID)
          }
        }
        assimVar.num=nOut
        assimVar = paste(assimVar.obsid[1],assimVar.simid[1],sep=",")
        if(nOut>1){
          for(i in 1:(nOut-1)){
            assimVar=c(assimVar,paste(assimVar.obsid[i+1],assimVar.simid[i+1],sep=","))
          }
        }
      }else{
        assimVar="9999,9999"
      }
      
    }else if(app.sys=="win"){
      # set to default values, if not set
      # set some test values for development on windows:
      if(!exists("win.idate")|
         !exists("win.outvars")|
         !exists("win.outbasins")|
         !exists("win.wlperiod")|
         !exists("win.statfile")|
         !exists("win.xobs")|
         !exists("win.assimOn")|
         !exists("win.assimVar")){
        idate     <- "2017-01-01"  # Forecast issue date
        outvars   <- "cout"        # Output variables
        outbasins <- "37"          # Output basins
        rpfile    <- "default"     # Return periods levels file
        xobs      <- "-9999"       # EO/Insitu data Xobs file
        xobsNum=0
        xobsURL=NULL
        assimOn   <- "off"         # Assimilation on/off
        assimVar  <- "9999,9999"   # Assimilation variables (pairs, obs/sim)
      }else{
        idate     <- win.idate  # Forecast issue date
        outvars   <- win.outvars  # Output variables
        outbasins <- win.outbasins  # Output basins
        rpfile    <- win.statfile  # Return periods levels file
        xobs      <- win.xobs  # EO/Insitu data Xobs file
        xobsNum=0
        xobsURL=NULL
        assimOn   <- win.assimOn  # Assimilation on/off
        assimVar  <- win.assimVar  # Assimilation variables (pairs, obs/sim)
      }
    }else{
      idate     <- NULL  # Forecast issue date
      outvars   <- NULL  # Output variables
      outbasins <- NULL  # Output basins
      rpfile    <- NULL  # Return periods levels file
      xobs      <- NULL  # EO/Insitu data Xobs file
      xobsNum   <- 0
      xobsURL   <- NULL
      assimOn   <- NULL  # Assimilation on/off
      assimVar  <- NULL  # Assimilation variables (pairs, obs/sim)
      
      print("WARNING: hypeapps.sys not set, allowed values are 'tep' or 'win' ")
    }
    appInput=list("idate"     = idate,     # Forecast issue date
                  "outvars"   = outvars,   # Output variables
                  "outvars.name" = outvars.name,
                  "outvars.name2" = outvars.name2,
                  "outvars.unit" = outvars.unit,
                  "outvars.id"   = outvars.id,
                  "outvars.num"  = outvars.num,
                  "outbasins" = basins, # Output basins
                  "outbasins.name"=basins.name,
                  "outbasins.id"=basins.id,
                  "outbasins.num"=basins.num,
                  "rpfile"    = rpcout,    # Return periods levels file
                  "xobs"      = xobs,      # EO/Insitu data Xobs file
                  "xobsNum"   = xobsNum,
                  "xobsURL"   = xobsURL,
                  "assimOn"   = assimOn,   # Assimilation on/off
                  "assimVar"  = assimVar)  # Assimilation variables (pairs, obs/sim)
    
  }else if(appName=="returnperiod"){
    if(app.sys=="tep"){
      
      # get parameters and data files
      timeFileIN <- rciop.getparam("timeOutputData") # time output input files
      returnPeriodIN   <- rciop.getparam("returnPeriod")   # return period inputs
      
      # parse the list of input timeXXXX files
      # If xobs !=-9999, parse the input to URLs
      if(nchar(timeFileIN)>0 & timeFileIN[1]!="-9999"){
        timeFileNum=0
        timeFileURL=NULL
        for(i in 1:length(timeFileIN)){
          if(nchar(timeFileIN[i])>1){
            
            # update the number of timeFileNum inputs
            timeFileNum=timeFileNum+1
            
            # extract the URL from input string:
            if(timeFileNum==1){
              timeFileURL=strsplit(timeFileIN[i],split = "&")[[1]][1]
            }else{
              timeFileURL=c(timeFileURL,strsplit(timeFileIN[i],split = "&")[[1]][1])
            }
          }
        }
      }else{
        timeFileNum=0
        timeFileURL=NULL
      }
      
      # parse the list of return perids
      returnPeriod = as.integer(trimws(strsplit(returnPeriodIN,split = ",")[[1]]))
      
      appInput=list("timeFileNum"=timeFileNum,"timeFileURL"=timeFileURL,
                    "returnPeriod"=returnPeriod)
    }else{
      appInput=list("timeFileNum"=NULL,"timeFileURL"=NULL,
                    "returnPeriod"=NULL)
    }
    
  }else if(appName=="eodata"){
    
    if(app.sys=="tep"){
      # get parameters and data files
      wlDataIn   <- as.character(rciop.getparam("wlData"))     # Water level data
      wlSubid    <- as.character(rciop.getparam("wlSubid"))    # water level data subbasin identifiers
      wlVariable <- as.character(rciop.getparam("wlVariable")) # water level Xobs variable name
      wlOffset   <- as.character(rciop.getparam("wlOffset"))   # water level offset correction (m)
      
      wlDataInput=F
      
      # Flood Monitoring data commented out for now, DG 20170615      
      #      flData     <- as.character(rciop.getparam("flData"))      # flood map data
      #      flVariable <- as.character(rciop.getparam("flVariable"))  # flood map Xobs variable name
      
      # parse the water level inputs, first the data URLs:
      if(wlDataIn!="-9999" & nchar(wlDataIn)>0){
        wlData = strsplit(x = wlDataIn,split = ",")[[1]]
        wlDataNum=0
        for(i in 1:length(wlData)){
          if(nchar(wlData[i])>1){
            # update the number of wldata inputs
            wlDataNum=wlDataNum+1
            
            # extract the URL from wlData input string:
            if(wlDataNum==1){
              wlDataURL=strsplit(wlData[i],split = "&")[[1]][1]
            }else{
              wlDataURL=c(wlDataURL,strsplit(wlData[i],split = "&")[[1]][1])
            }
          }
        }
        # secondly the SUBIDs and Variables 
        if(wlDataNum>0){
          # SUBIDs:
          if(nchar(wlSubid)>0){
            wlDataSubid=as.integer(strsplit(wlSubid,",")[[1]])
          }else{
            wlDataSubid=-9999
          }
          if(length(wlDataSubid)==wlDataNum & length(wlDataURL)==wlDataNum & wlDataSubid[1]!=-9999){
            # flag that input is ok
            wlDataInput=T
            
            # variables:
            wlDataVariableTemp=trimws(strsplit(wlVariable,",")[[1]])
            nVar=length(wlDataVariableTemp)
            for(i in 1:nVar){
              splitFirst = trimws(strsplit(wlDataVariableTemp[i],split = "(",fixed=T)[[1]])
              oName=splitFirst[1]
              splitSecond = trimws(strsplit(splitFirst[2],split = ")",fixed=T)[[1]])
              oUnit=splitSecond[1]
              oID=splitSecond[2]
              if(i==1){
                wlDataVariable.name=oName
                wlDataVariable.unit=oUnit
                wlDataVariable=oID
              }else{
                wlDataVariable.name=c(wlDataVariable.name,oName)
                wlDataVariable.unit=c(wlDataVariable.unit,oUnit)
                wlDataVariable=c(wlDataVariable,oID)
              }
            }

            # Offsets:
            wlDataOffset=as.integer(strsplit(wlOffset,",")[[1]])
            if(length(wlDataVariable)<wlDataNum){
              wlDataVariable=c(wlDataVariable,rep(wlDataVariable[1],wlDataNum-length(wlDataVariable)))
            }else if(length(wlDataVariable)>wlDataNum){
              wlDataVariable=wlDataVariable[1:wlDataNum]
            }
            if(length(wlDataOffset)<wlDataNum){
              wlDataOffset=c(wlDataOffset,rep(wlDataOffset[1],wlDataNum-length(wlDataOffset)))
            }else if(length(wlDataOffset)>wlDataNum){
              wlDataOffset=wlDataOffset[1:wlDataNum]
            }
          }else{
            wlDataNum=0
            wlDataInput=F
            wlDataURL=NULL
            wlDataSubid=NULL
            wlDataVariable=NULL
            wlDataVariable.name <- NULL
            wlDataVariable.unit <- NULL
            wlDataOffset=NULL
          }
        
        }else{
          wlDataNum=0
          wlDataInput=F
          wlDataURL=NULL
          wlDataSubid=NULL
          wlDataVariable=NULL
          wlDataVariable.name <- NULL
          wlDataVariable.unit <- NULL
          wlDataOffset=NULL
        }
      
      #      # parse the flood mapping inputs
      #      if(nchar(flData)>1){
      #        flDataInput=T
      #        
      #        # extract the URL from wlData input string:
      #        flDataURL=strsplit(flData[1],split = "&")[[1]][1]
      #        
      #        # subid list
      #        flDataSubid = flSubid
      #        
      #      }else{
      #      flDataNum=0
      #      flDataInput=F
      #      flDataURL=NULL
      #      flDataSubid=NULL
      #      }
      }else{
        wlDataNum    <- 0
        wlDataInput  <- F
        wlDataURL    <- NULL
        wlDataSubid  <- NULL
        wlDataVariable <- NULL
        wlDataVariable.name <- NULL
        wlDataVariable.unit <- NULL
        wlDataOffset <- NULL
      #     flDataNum    <- 0
      #     flDataInput  <- F
      #     flDataURL    <- NULL
      #     flDataSubid  <- NULL
      }
    }else{
      wlDataNum=0
      wlDataInput=F
      wlDataURL=NULL
      wlDataSubid=NULL
      wlDataVariable=NULL
      wlDataVariable.name <- NULL
      wlDataVariable.unit <- NULL
      wlDataOffset=NULL
    }
    appInput=list("wlDataNum"=wlDataNum,
                  "wlDataInput"=wlDataInput,
                  "wlDataURL"=wlDataURL,
                  "wlDataSubid"=wlDataSubid,
                  "wlDataVariable"=wlDataVariable,
                  "wlDataVariable.name"=wlDataVariable.name,
                  "wlDataVariable.unit"=wlDataVariable.unit,
                  "wlDataOffset"=wlDataOffset) #,
    #  "flDataNum"=flDataNum,
    #  "flDataInput"=flDataInput,
    #  "flDataURL"=flDataURL,
    #  "flDataSubid"=flDataSubid)
    
  }else{
    appInput=list("appName"=NULL)
  }
  return(appInput)
}

## -------------------------------------------------------------------------------
## prepare work directories and copy basic model files
getHypeAppSetup<-function(modelName,modelBin,tmpDir,appDir,appName,appInput,modelFilesURL,forcingArchiveURL=NULL,stateFilesURL=NULL,stateFilesIN=NULL){
  
  ## model files run directory (for all applications, except returnperiod)
  if(appName=="historical"|appName=="forecast"|appName=="eodata"|appName=="returnperiod"){
    modelFilesRunDir=paste(tmpDir,'model',modelName,sep="/")
    dir.create(modelFilesRunDir,recursive = T,showWarnings = F)
  }else{
    modelFilesRunDir=NULL
  }
  
  # model files results directory (not necessary for "eodata" and "returnperiod")
  if(appName=="historical"){
    modelResDir=paste(modelFilesRunDir,'results',sep="/")
    dir.create(modelResDir,recursive = T,showWarnings = F)
  }else if(appName=="forecast"){
    modelResDir=paste(modelFilesRunDir,'hindcast',sep="/")
    modelResDir=c(modelResDir,paste(modelFilesRunDir,'forecast',sep="/"))
    for(i in 1:2){
      dir.create(modelResDir[i],recursive = T,showWarnings = F)
    }
  }else{
    modelResDir=NULL
  }
  
  # eodata results directory
  if(appName=="eodata"){
    eodataResDir=paste(modelFilesRunDir,"eodata",sep="/")
    dir.create(eodataResDir,recursive = T,showWarnings = F)
  }else{
    eodataResDir=NULL
  }
  
  # return period run directory
  if(appName=="returnperiod"){
    returnperiodResDir=paste(modelFilesRunDir,"returnperiod",sep="/")
    dir.create(returnperiodResDir,recursive = T,showWarnings = F)
  }else{
    returnperiodResDir=NULL
  }
  
  # copy model files to working directory
  if(appName=="historical"|appName=="forecast"|appName=="eodata"){
    fileNames=c("par.txt",
                "GeoData.txt",
                "GeoClass.txt",
                "BranchData.txt",
                "FloodData.txt",
                "LakeData.txt")
    
    if(appName=="historical"|appName=="eodata"){
      fileNames=c(fileNames,"info-historical.txt")
    }else if(appName=="forecast"){
      fileNames=c(fileNames,"info-hindcast.txt","info-forecast.txt")
    }
    
    if(appName=="historical"|appName=="forecast"){
      if(appInput$assimOn=="on"){
        fileNames=c(fileNames,"info-hindcast-assimilation.txt","info-historical-assimilation.txt","AssimInfo-AOWL.txt","AssimInfo-Openloop.txt","AssimInfo-Openloop-inibin.txt")
      }
    }
    
    for(i in 1:length(fileNames)){
      if(app.sys=="tep"){
        res <- rciop.copy(paste(modelFilesURL,fileNames[i],sep="/"), modelFilesRunDir, uncompress=TRUE)
      }else{
        file.copy(from=paste(modelFilesURL,fileNames[i],sep="/"),
                  to =paste(modelFilesRunDir,fileNames[i],sep="/"),
                  overwrite = T)
      }
    }
  }
  
  ## model binary file (stays in application folder)
  if(appName=="historical"|appName=="forecast"){
    # model binary source file
    modelBinaryFile=paste(appDir,'util/bin',modelBin,sep="/")
    
    # command line to run the model "hype rundir"
    sysCommand = paste(paste(modelBinaryFile,modelFilesRunDir,sep=" "),"/",sep="")
    
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("HYPE RUN system command=",sysCommand,sep=""), "/util/R/hypeapps-utils.R")}
  }else{
    modelBinaryFile=NULL
    sysCommand=NULL
  }
  
  ## State file filenames and dates, set bdate limits
  if(appName=="historical"|appName=="forecast"){
    if(app.sys=="tep"){
      stateFiles = stateFilesIN
      stateDates = as.POSIXct(paste(substr(stateFiles,11,14),substr(stateFiles,15,16),substr(stateFiles,17,18),sep="-"),tz="GMT")
      bdateMax   = max(stateDates)
      bdateMin   = min(stateDates)
    }else{
      stateFiles = NULL
      stateDates = NULL
      bdateMax   = as.POSIXct("2017-01-01", tz="GMT")
      bdateMin   = as.POSIXct("1979-01-01", tz="GMT")
    }
  }else{
    stateFiles = NULL
    stateDates = NULL
    bdateMax   = NULL
    bdateMin   = NULL
  }
  
  ## check existance of forcing archive, and if existing, it's first and last date.
  if(!is.null(forcingArchiveURL)){
    forcingArchiveExist=T
  }else{
    forcingArchiveExist=F
  }
  
  ## Sub-basin shapefiles (potentially for all applications for map plots)
  if(!is.null(shapefile.url)){
    # libraries needed to read the shapefile
    library(sp)
    library(rgdal)
    library(rgeos)
    
    # subfolder for shapefiles
    shapefileDir = paste(modelFilesRunDir,"shapefile",sep="/")
    dir.create(shapefileDir,recursive = T,showWarnings = F)
    
    #download shapefile from storage
    for(i in 1:length(shapefile.ext)){
      rciop.copy(paste(paste(shapefile.url,shapefile.layer,sep="/"),shapefile.ext[i],sep=""), shapefileDir)
    }
    
    # open and save shapefile as Rdata
    shapefileData = readOGR(dsn = shapefileDir, layer = shapefile.layer)
    shapefileRdata = paste(shapefileDir,"/",shapefile.layer,".Rdata",sep="")
    save(list = "shapefileData",file = shapefileRdata)
    
  }else{
    shapefileRdata=NULL
    shapefileData=NULL
    shapefileDir=NULL
  }
  
  ## return period magnitudes default files OR file from input
  if(appName=="forecast"){
    rpFileCOUT=NULL

    rciop.log ("DEBUG", paste(" appInput$rpfile= ", appInput$rpfile, sep=""), "getHypeSetup")
    
    if(appInput$rpfile=="default"){
      # download default file from data storage
      rpFileURL = paste(modelFilesURL,"returnlevels",paste(modelName,"-rp-cout.txt",sep=""),sep="/")
      # download rpfile to forecast output folder - using rciop.copy since we already have the URL to the file
      rciop.copy(rpFileURL, modelResDir[2])
      # path to downloaded rpfile
      rpFileCOUT = paste(modelResDir[2],paste(paste(modelName,"-rp-cout.txt",sep="")),sep="/")
    }else{
      # download the file specified by user input
      #
      # in this case we have to use the opensearch-client command, since the user 
      # input is a opensearch URL and not the URL to the file
      
      # make subfolder for download:
      targetFolder = paste(tmpDir,"/rpFile_1",sep="")
      dir.create(targetFolder,recursive = T,showWarnings = F)

      rciop.log ("DEBUG", paste(" targetFolder = ", targetFolder, sep=""), "getHypeSetup")
      
      # get file using opensearch-client
      sysCmd=paste("opensearch-client '",appInput$rpfile,"' enclosure | ciop-copy -s -U -O ",targetFolder,"/ -", sep="")
      
      rciop.log ("DEBUG", paste(" sysCmd = ", sysCmd, sep=""), "getHypeSetup")
      
      rpFileCOUT=system(command = sysCmd,intern = T)
      
      rciop.log ("DEBUG", paste(" rpFileCOUT = ", rpFileCOUT, sep=""), "getHypeSetup")
      
    }
    # check existance of rpFileCOUT and check available return periods in the file
    if(!file.exists(rpFileCOUT)){
      rpFileCOUT=NULL
    }
    
  }else{
    rpFileCOUT=NULL
  }
  
  ## return list with application setup
  appSetup = list("runDir"=modelFilesRunDir,
                  "resDir"=modelResDir,
                  "runCommand"=sysCommand,
                  "eodataResDir"=eodataResDir,
                  "returnperiodResDir"=returnperiodResDir,
                  "tmpDir"=tmpDir,
                  "appDir"=appDir,
                  "appName"=appName,
                  "modelName"=modelName,
                  "modelBin"=modelBin,
                  "stateFilesURL"=stateFilesURL,
                  "stateFiles"=stateFiles,
                  "stateDates"=stateDates,
                  "forcingArchiveURL"=forcingArchiveURL,
                  "forcingArchiveExist"=forcingArchiveExist,
                  "hype2csvURL"=hype2csv.url,
                  "hype2csvFile"=hype2csv.file,
                  "shapefileDir"=shapefileDir,
                  "shapefileLayer"=shapefile.layer,
                  "shapefileExt"=shapefile.ext,
                  "shapefileRdata"=shapefileRdata,
                  "shapefileData"=shapefileData,
                  "rpFileCOUT"=rpFileCOUT)
  return(appSetup)
}

## -------------------------------------------------------------------------------
## get eo data from URL
getEoData<-function(appInput,appSetup){
  
  #  appInput=list("wlDataNum"=wlDataNum
  #                "wlDataInput"=wlDataInput,
  #                "wlDataURL"=wlDataURL,
  #                "wlDataSubid"=wlDataSubid,
  #                "wlDataVariable"=wlDataVariable,
  #                "wlDataOffset"=wlDataOffset,
  #                "flDataNum"=flDataNum  
  #                "flDataInput"=flDataInput,
  #                "flDataURL"=flDataURL,
  #                "flDataSubid"=flDataSubid)
  
  #  appSetup = list("runDir"=modelFilesRunDir,
  #                  "resDir"=modelResDir,
  #                  "runCommand"=sysCommand,
  #                  "eodataResDir"=eodataResDir,
  #                  "tmpdir"=runDir)
  
  # loop over wlDataInput
  nDownLoad=0
  if(appInput$wlDataNum>0){
    for(i in 1:appInput$wlDataNum){
      rciop.log("DEBUG", paste("i=",as.character(i),sep=""), "getEoData")
      
      # make subfolder for download:
      targetFolder = paste(appSetup$tmpDir,"/wlData_",as.character(i),sep="")
      dir.create(targetFolder,recursive = T,showWarnings = F)

      rciop.log("DEBUG", paste("targetFolder=",targetFolder,sep=""), "getEoData")
      
      # get file using opensearch-client
      sysCmd=paste("opensearch-client '",appInput$wlDataURL[i],"' enclosure | ciop-copy -s -U -O ",appSetup$tmpDir,"/wlData_",as.character(i),"/ -", sep="")
      
      rciop.log("DEBUG", paste("sysCmd=",sysCmd,sep=""), "getEoData")
      
      wlFile=system(command = sysCmd,intern = T)
      
      if(file.exists(wlFile)){
        nDownLoad=nDownLoad+1
        if(nDownLoad==1){
          wlFiles=wlFile
          wlSubid=appInput$wlDataSubid[i]
          wlVariable=appInput$wlDataVariable[i]
          wlVariable.name=appInput$wlDataVariable.name[i]
          wlVariable.unit=appInput$wlDataVariable.unit[i]
          wlOffset=appInput$wlDataOffset[i]
        }else{
          wlFiles=c(wlFiles,wlFile)
          wlSubid=c(wlSubid,appInput$wlDataSubid[i])
          wlVariable=c(wlVariable,appInput$wlDataVariable[i])
          wlVariable.name=c(wlVariable.name,appInput$wlDataVariable.name[i])
          wlVariable.unit=c(wlVariable.unit,appInput$wlDataVariable.unit[i])
          wlOffset=c(wlOffset,appInput$wlDataOffset[i])
        }
      }
    }
    eoData=list("wlDataNum"=nDownLoad,
                "wlDataFile"=wlFiles,
                "wlDataSubid"=wlSubid,
                "wlDataVariable"=wlVariable,
                "wlDataVariable.name"=wlVariable.name,
                "wlDataVariable.unit"=wlVariable.unit,
                "wlDataOffset"=wlOffset)
  }else{
    eoData=list("wlDataNum"=0,
                "wlDataFile"=NULL,
                "wlDataSubid"=NULL,
                "wlDataVariable"=NULL,
                "wlDataVariable.name"=NULL,
                "wlDataVariable.unit"=NULL,
                "wlDataOffset"=NULL)
  }
  #      #wlData="http://sb-10-15-36-31.hydro.terradue.int/sbws/production/run/water-levels/0000081-170411124531167-oozie-oozi-W/products/search?uid=0000081-170411124531167-oozie-oozi-W/outputs/lakes_summary_multi_ATK_Kainji_L2.csv&format=atom"
  #      if(nchar(wlData)>1){
  #        # download data
  #        sysCmd=paste("opensearch-client",wlDataURL,"enclosure | ciop-copy -s -U -O /tmp/ -", sep=" ")
  #        wlFile=system(command = sysCmd,intern = T)
  
  return(eoData)
}

## -------------------------------------------------------------------------------
## Read eodata and prepare data on the Xobs format
readEoData<-function(appSetup,eoData){
  for(i in 1:eoData$wlDataNum){
    xobs = readWaterLevelCSV(fname=eoData$wlDataFile[i],
                             subid=eoData$wlDataSubid[i],
                             vars=eoData$wlDataVariable[i],
                             wl0=eoData$wlDataOffset[i],
                             t1=NULL,
                             t2=NULL,
                             modelSHP=NULL,
                             limData=T,
                             lowLim=-0.9999,
                             highLim=9999)
    if(i==1){
      xobsData=xobs
    }else{
      #      xobsData= mergeXobs
    }
  }
  return(xobsData)
}

## -------------------------------------------------------------------------------
## write eo data in Xobs format
writeEoData<-function(appSetup,xobsData,appDate,prefix="001"){
  
  # Create folder for data to be published
  outDir = paste(appSetup$tmpDir,'output',sep="/")
  dir.create(outDir,recursive = T,showWarnings = F)
  
  # variable and subid to include in output filename (we use only the first variable and subid)
  varName=attr(xobsData,"variable")[1]
  subID=as.character(attr(xobsData,"subid")[1])
  
  # output filename
  xobsFile=paste(outDir,"/",prefix,"_",appDate,"_Xobs_",varName,"_",subID,".txt",sep="")
  
  # write to file
  outres = WriteXobs(xobsData, filename = xobsFile)
  
  # return output list
  appOutput = list("outDir"=outDir,"files"=xobsFile)
  return(appOutput)
}

## -------------------------------------------------------------------------------
## get Xobs data from open catalogue
getXobsData<-function(appInput,appSetup){
  
  #  appInput=list("xobsNum"=xobsNum
  #                "xobs"=xobs,
  #                "xobsURL"=xobsURL)
  
  #  appSetup = list("runDir"=modelFilesRunDir,
  #                  "resDir"=modelResDir,
  #                  "runCommand"=sysCommand,
  #                  "tmpdir"=runDir)
  
  # loop over xobs
  nDownLoad=0
  if(appInput$xobsNum>0){
    for(i in 1:appInput$xobsNum){
      
      # make subfolder for download:
      targetFolder = paste(appSetup$tmpDir,"/xobsData_",as.character(i),sep="")
      dir.create(targetFolder,recursive = T,showWarnings = F)
      
      # get file using opensearch-client
      sysCmd=paste("opensearch-client '",appInput$xobsURL[i],"' enclosure | ciop-copy -s -U -O ",appSetup$tmpDir,"/xobsData_",as.character(i),"/ -", sep="")
      xobsFile=system(command = sysCmd,intern = T)
     
      if(file.exists(xobsFile)){
        nDownLoad=nDownLoad+1
        # add xobs file id number
        xobsFileNew=paste(substr(xobsFile,1,nchar(xobsFile)-4),"_",as.character(i),".txt",sep="")
        file.copy(from = xobsFile, to = xobsFileNew,overwrite = T)
        if(nDownLoad==1){
          xobsFiles=xobsFileNew
        }else{
          xobsFiles=c(xobsFiles,xobsFileNew)
        }
      }
    }
    xobsData=list("xobsNum"=nDownLoad,
                  "xobsFile"=xobsFiles)
  }else{
    xobsData=list("xobsNum"=0,
                  "xobsFile"=NULL)
  }
  return(xobsData)
}

## -------------------------------------------------------------------------------
## read downloaded Xobs input file(s) - merge into one Xobs.txt in the model run folder
readXobsData<-function(appSetup,xobsData){
  
  # copy the first downloaded file into model run directory
  if(xobsData$xobsNum>0){
    if(file.exists(xobsData$xobsFile[1])){
      # copy file
      file.copy(from = xobsData$xobsFile[1],to = paste(appSetup$runDir,"Xobs.txt",sep="/") , overwrite = T)
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" Xobs.txt copied to runDir from  >> ",xobsData$xobsFile[1],sep=""), "/util/R/hypeapps-utils.R")}
      # read xobs file and extract the variable names and subbasin ids
      xobs.data = ReadXobs(filename = xobsData$xobsFile[1])
      xobs.var  = unique(attr(xobs.data,"variable"))
      xobs.subid= unique(attr(xobs.data,"subid"))
      if(is.na(xobs.var[1])){
        if(app.sys=="tep"){rciop.log("DEBUG", paste(" somethings wrong with Xobs.txt content, xobs.var =  ",xobs.var,sep=""), "/util/R/hypeapps-utils.R")}
        if(app.sys=="tep"){rciop.log("DEBUG", paste(" somethings wrong with Xobs.txt content, xobs.subid =  ",xobs.subid,sep=""), "/util/R/hypeapps-utils.R")}
      }else{
        if(app.sys=="tep"){rciop.log("DEBUG", paste(" Xobs.txt content, xobs.var =  ",xobs.var,sep=""), "/util/R/hypeapps-utils.R")}
        if(app.sys=="tep"){rciop.log("DEBUG", paste(" Xobs.txt content, xobs.subid =  ",xobs.subid,sep=""), "/util/R/hypeapps-utils.R")}
        # output list
        xobsInput = list("xobsFile"=paste(appSetup$runDir,"Xobs.txt",sep="/"),
                       "xobsVar"=xobs.var,
                       "xobsSubid"=xobs.subid)
      }
    }else{
      if(app.sys=="tep"){rciop.log("DEBUG", paste(" Xobs.txt did not exist: ",xobsData$xobsFile[1],sep=""), "/util/R/hypeapps-utils.R")}

      xobsInput = list("xobsFile"=NA,
                       "xobsVar"=NA,
                       "xobsSubid"=NA)
    }
  }else{
    if(app.sys=="tep"){rciop.log("DEBUG", paste(" xobsNum = 0, no xobs file given on input: ",as.character(xobsData$xobsNum),sep=""), "/util/R/hypeapps-utils.R")}
    
    xobsInput = list("xobsFile"=NA,
                     "xobsVar"=NA,
                     "xobsSubid"=NA)
  }
  return(xobsInput) 
}

## -------------------------------------------------------------------------------
## getTimeOutputData - get time ouput data from data catalogue
getTimeOutputData<-function(appInput,appSetup){
  
  # input variables:
  # appInput$timeFileNum
  # appInput$timeFileURL
  # appInput$returnPeriod=returnPeriod
  if(app.sys=="tep"){rciop.log ("DEBUG", paste("...entering getTimeOutput ","!",sep=""), "/util/R/hypeapps-utils.R")}
  if(app.sys=="tep"){rciop.log ("DEBUG", paste("... appInput$timeFileNum: ",as.character(appInput$timeFileNum),sep=""), "/util/R/hypeapps-utils.R")}
  if(app.sys=="tep"){rciop.log ("DEBUG", paste("... appInput$timeFileURL ",appInput$timeFileURL,sep=""), "/util/R/hypeapps-utils.R")}
  
  
  # loop over timeFiles
  nDownLoad=0
  if(appInput$timeFileNum>0){
    for(i in 1:appInput$timeFileNum){
      
      # make subfolder for download:
      targetFolder = paste(appSetup$tmpDir,"/timeData_",as.character(i),sep="")
      dir.create(targetFolder,recursive = T,showWarnings = F)
      
      # get file using opensearch-client
      sysCmd=paste("opensearch-client '",appInput$timeFileURL[i],"' enclosure | ciop-copy -s -U -O ",appSetup$tmpDir,"/timeData_",as.character(i),"/ -", sep="")
      timeFile=system(command = sysCmd,intern = T)
      
      if(file.exists(timeFile)){
        if(app.sys=="tep"){rciop.log ("DEBUG", paste("... timeFile ", timeFile,sep=""), "/util/R/hypeapps-utils.R")}
        
        nDownLoad=nDownLoad+1
        # add timeFile id number in case several files have the same name
        timeFileNew=paste(substr(timeFile,1,nchar(timeFile)-4),"_",as.character(i),".txt",sep="")
        file.copy(from = timeFile, to = timeFileNew,overwrite = T)
        if(nDownLoad==1){
          timeFiles=timeFileNew
        }else{
          timeFiles=c(timeFiles,timeFileNew)
        }
      }
    }
    timeData=list("timeFileNum"=nDownLoad,
                  "timeFile"=timeFiles)
  }else{
    timeData=list("timeFileNum"=0,
                  "timeFile"=NULL)
  }
  return(timeData)
}

## -------------------------------------------------------------------------------
## analyseTimeOutputData - return period analysis on timeOutput data
analyseTimeOutputData<-function(appSetup,appInput,timeData,appDate){
  
  # loop over timeFiles and make the return period analysis
  if(timeData$timeFileNum>0){
    for(i in 1:timeData$timeFileNum){
      timeFileIn  = timeData$timeFile[i]
      timeFileName = strsplit(timeFileIn,split = "/")[[1]]
      timeFileName = timeFileName[length(timeFileName)]
      timeFileOut = paste(appSetup$returnperiodResDir,paste("rp_",appDate,"_",timeFileName,sep=""),sep="/")
      
      if(app.sys=="tep"){rciop.log ("DEBUG", paste("   timeFileIn=",timeFileIn,sep=""), "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste("   timeFileName=",timeFileName,sep=""), "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste("   timeFileOut=",timeFileOut,sep=""), "/util/R/hypeapps-utils.R")}
      
      # analyse file
      returnPeriodMagnitudes(name.in=timeFileIn,name.out=timeFileOut,wl.rp=appInput$returnPeriod,dist="gev")
      
      # save list of output files
      if(i==1){
        timeFileOutAll=timeFileOut
      }else{
        timeFileOutAll=c(timeFileOutAll,timeFileOut)
      }
    } 
  }
  
  appOutput=list("outDir"=appSetup$returnperiodResDir,"files"=timeFileOutAll)
  
  
  return(appOutput)
}


## -------------------------------------------------------------------------------
## DATE2INFODATE - function to return a date as a textstring suitable for info.txt
DATE2INFODATE<-function(datepsx){
  datetxt=as.character(datepsx)
  infodate=paste(substr(datetxt,1,4),substr(datetxt,6,7),substr(datetxt,9,10),sep="")
  return(infodate)
}

## -------------------------------------------------------------------------------
## ID2ZF - function to return gfd zipfile name from a requested issuedate (POSIX)
ID2ZF<-function(idpsx){
  idtxt=as.character(idpsx)
  zfname=paste(substr(idtxt,1,4),substr(idtxt,6,7),substr(idtxt,9,10),".zip",sep="")
  return(zfname)
} 

## -------------------------------------------------------------------------------
## ymd2str - function to return date string "yyyy-mm-dd" from integer year, month day
ymd2str<-function(yr,mn,da){
  dateStr=as.character(yr)
  if(mn>=10){
    dateStr=paste(dateStr,"-",as.character(mn),sep="")
  }else{
    dateStr=paste(dateStr,"-0",as.character(mn),sep="")
  }
  if(da>=10){
    dateStr=paste(dateStr,"-",as.character(da),sep="")
  }else{
    dateStr=paste(dateStr,"-0",as.character(da),sep="")
  }
  return(dateStr)
}


## -------------------------------------------------------------------------------
## ymd2posix - function to return posix date from integer year, month day
ymd2posix<-function(yr,mn,da){
  dateStr=ymd2str(yr,mn,da)
  dateNum=as.POSIXct(dateStr,tz="GMT")
  return(dateNum)
}

## -------------------------------------------------------------------------------
## getGFDzipFromTep - download, archive, hindcast or forecast GFD data from TEP
getGFDzipFromTep<-function(issueDateNum=as.POSIXct("2017-01-01", tz = "GMT"),
                           modelName="niger-hype",productName="hindcast",
                           archive=F,gfdVersion="v1",tmpDir="/tmp"){
  # remote folder
  remoteFolder = paste("https://store.terradue.com//smhi/gfd",modelName,productName,"files",gfdVersion,sep="/")
  
  # zipfile name
  zipFile    = ID2ZF(issueDateNum)
  
  # zipfile URLs (local and remote)
  localFile  = paste(tmpDir,zipFile,sep="/")
  remoteFile = paste(remoteFolder,zipFile,sep="/")
  
  # delete local file if it existst
  if(file.exists(localFile)){file.remove(localFile)}
  
  # system command to download remote file to local path
  sysCmd=paste("curl -o",localFile,remoteFile,sep=" ")
  
  # download by system call
  if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying command >> ",sysCmd,sep=""), "/util/R/hypeapps-utils.R")}
  a=system(sysCmd,intern=T)
  
  # check existance and size of the downloaded file
  if(file.exists(localFile)){
    if(file.size(localFile)>1000){
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" file size OK >> ",as.character(file.size(localFile)),sep=""), "/util/R/hypeapps-utils.R")}
      idataOK=T
    }else{
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" file size NOT OK >> ",as.character(file.size(localFile)),sep=""), "/util/R/hypeapps-utils.R")}
      idataOK=F
    }
  }else{
    idataOK=F
  }
  
  # return list
  resList=list("status"=idataOK,"localFile"=localFile,"remoteFile"=remoteFile)
  return(resList)
}

## -------------------------------------------------------------------------------
## getHindcastForcingData - function to create forcing data for hindcast period
##                          combining archive and downloadable hindcast data
##
##  output = list with forcing files
getHindcastForcingData<-function(startDate,endDate,appSetup,obsFiles,outDir,useRdata=F){
  
  # 1. check if Archive and Hindcasts are needed
  useArchive  = F ; if(startDate<=forcing.archive.end){useArchive  = T}
  useHindcast = F ; if(endDate>forcing.archive.end){useHindcast  = T}
  
  # 2. Archive data - if needed
  if(useArchive){
    
    # log message
    if(app.sys=="tep"){rciop.log ("DEBUG", " getHindcastForcingData: useArchive = TRUE ...", "/util/R/hypeapps-utils.R")}
    
    # create tmpDir/forcing/archive
    archiveDir = paste(appSetup$tmpDir,"/forcing/archive",sep="")
    dir.create(archiveDir,recursive = T,showWarnings = F)
    
    # archive exists
    if(appSetup$forcingArchiveExist){
      # copy forcing files from local archive
      for(i in 1:length(forcing.files)){
        # copy text files, if only archive is needed
        rciop.copy(url = paste(appSetup$forcingArchiveURL,obsFiles[i],sep="/"),
                   target = archiveDir)
      }
      archiveFound = T
    }else{
      # download archive from data catalogue
      sysCmd=paste("curl -o ", appSetup$tmpDir,"/archive.zip https://store.terradue.com//smhi/gfd/niger-hype/hindcast/files/v1/archive.zip",sep="")
      a=system(sysCmd,intern=T)
      # unzip to forcing/archive
      archiveFile=paste(appSetup$tmpDir,"archive.zip",sep="/")
      if(file.exists(archiveFile)){
        unzip(zipfile = archiveFile,overwrite = T,exdir = archiveDir)
        archiveFound = T
      }else{
        return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=T))
        archiveFound = F
      }
    }
  }
  
  # 3. Hindcast data - if needed
  if(useHindcast){
    # read daily hindcasts data from repository to reach the endDate 
    # with as updated data as possible
    
    # first check the issue date of the latest available hindcast
    today.Str  = Sys.Date()
    today.Num   = as.POSIXct(today.Str,tz = "GMT")
    issueDate.Num = today.Num
    
    # try downloading hindcast with issue date equal to the requested
    downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                    modelName=appSetup$modelName,
                                    productName="hindcast",
                                    tmpDir=appSetup$tmpDir)
    if(!downloadInfo$status){
      # try download data issued issueDateNum-1:9
      i=1
      while(i<20 & !downloadInfo$status){
        # search issue data backwards
        issueDate.Num = today.Num-i*86400
        downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                        modelName=appSetup$modelName,
                                        productName="hindcast",
                                        tmpDir=appSetup$tmpDir)
        i=i+1
      }
    }
    # if data downloaded, unzip
    if(downloadInfo$status){
      # save latest issue date
      lastIssueDate.Num=issueDate.Num
      
      #unzip downloaded file
      hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(lastIssueDate.Num),sep="")
      dir.create(hindcastDir,recursive = T,showWarnings = F)
      unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = hindcastDir)
      
      # check first and last date in file
      firstDate.last=substr(system(paste("gawk 'NR==2{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
      lastDate.last=substr(system(paste("gawk 'END{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
      
      # adjust edate to laste date in file, if originally later
      endDate = min(endDate,as.POSIXct(lastDate.last,tz="GMT"))
      
      # download files covering 2017-01-01 until the endDate is within the first month of the file (or until we reach the last file)
      notFound = T
      y1=2017
      m1=1
      d1=1
      
      nIssueDates=0
      
      while(notFound){
        issueDate.Num = ymd2posix(y1,m1,d1)+123*86400
        # try downloading hindcast with issue date equal to the requested
        downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                        modelName=appSetup$modelName,
                                        productName="hindcast",
                                        tmpDir=appSetup$tmpDir)
        if(downloadInfo$status){
          # save issueDate to vector of successful dates
          if(nIssueDates==0){
            issueDates = issueDate.Num
          }else{
            issueDates=c(issueDates,issueDate.Num)
          }
          nIssueDates = nIssueDates+1
          
          #unzip downloaded file
          hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(issueDate.Num),sep="")
          dir.create(hindcastDir,recursive = T,showWarnings = F)
          unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = hindcastDir)
          
          # check first and last date in file
          firstDate=substr(system(paste("gawk 'NR==2{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
          lastDate=substr(system(paste("gawk 'END{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
          
          # next month as integer
          y2=as.integer(substr(firstDate,1,4))
          m2=as.integer(substr(firstDate,6,7))+1
          d2=1
          if(m2>12){
            y2=y2+1
            m2=1
          }
          # next month as posix
          nextMonth=ymd2posix(y2,m2,d2)
          if(endDate<=nextMonth){
            notFound=F
            useLast=F
          }else{
            notFound=T
            useLast=F
            y1=y2
            m1=m2
            d1=d2
          }
        }else{
          notFound=F
          useLast=T
        }
      }
      # add last file if needed
      if(useLast){
        issueDates=c(issueDates,lastIssueDate.Num)
      }
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" issueDates = ", issueDates,sep=","), "/util/R/hypeapps-utils.R")}
      hindcastFound = T
      if(app.sys=="tep"){rciop.log ("DEBUG", " hindcast download ready!", "/util/R/hypeapps-utils.R")}
    }else{
      if(app.sys=="tep"){rciop.log ("DEBUG", " hindcast download failed!", "/util/R/hypeapps-utils.R")}
      hindcastFound = F
    }
  }
  
  # 4. Merge obsfiles and copy to runDir
  for(i in 1:length(obsFiles)){
    if(useArchive){
      if(archiveFound){
        obsData = ReadPTQobs(paste(archiveDir,obsFiles[i],sep="/"))
        iStart  = which(obsData[,1]<=as.POSIXct(startDate,tz="GMT"))
        iStart  = iStart[length(iStart)]
        iEnd    = which(obsData[,1]>=min(c(endDate,forcing.archive.end)))
        iEnd    = iEnd[1]
        obsData = obsData[iStart:iEnd,]
      }
    }
    if(useHindcast){
      if(hindcastFound){
        for(j in 1:length(issueDates)){
          hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(issueDates[j]),sep="")
          newData = ReadPTQobs(paste(hindcastDir,obsFiles[i],sep="/"))
          if(!useArchive & j==1){
            obsData=newData
          }else{
            iOld = which(obsData[,1]<newData[1,1])
            obsData = rbind(obsData[iOld,],newData)
          }
        }
      }
    }
    # final trimming to requested period
    obsData = obsData[which(obsData[,1]>=startDate & obsData[,1]<=endDate),]
    WritePTQobs(x = obsData,paste(outDir,obsFiles[i],sep="/"))
  }
  return(paste(outDir,obsFiles,sep="/"))
}

## -------------------------------------------------------------------------------
## mergeObsfiles - function to merge a set of obsfiles from two directories and write
##                 the result to a third directory given a startdate and enddate
##                 Files in dir1 has priority over dir2 files.
mergeObsFiles<-function(dir1,dir2,outDir,startDate,endDate,obsFiles){
  
  
  # loop over obsFiles
  for(i in 1:length(obsFiles)){
    
    if(i==1){
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... dir1 = ", dir1,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... dir2 = ", dir2,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
    }
    
    obs1 = ReadPTQobs(paste(dir1,obsFiles[i],sep="/"))
    obs2 = ReadPTQobs(paste(dir2,obsFiles[i],sep="/"))
    
    colnames(obs2)=colnames(obs1)
    
    if(i==1){
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... colnames1 = ", colnames(obs1[1:5]), sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... colnames2 = ", colnames(obs2[1:5]),sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
    }
    # clip obs1 and obs2 to the requested start and end dates
    obs1 = obs1[which(obs1[,1]>=startDate),]
    obs1 = obs1[which(obs1[,1]<=endDate),]
    obs2 = obs2[which(obs2[,1]>=startDate),]
    obs2 = obs2[which(obs2[,1]<=endDate),]
    
    # obs2 data before start of obs1
    i2 = which(obs2[,1]<obs1[1,1])
    
    if(length(i2)>=1){
      obs1 = rbind(obs2[i2,],obs1)
    }
    
    # obs2 data after end of obs1
    i2 = which(obs2[,1]>obs1[nrow(obs1),1])
    if(length(i2)>=1){
      obs1 = rbind(obs1,obs2[i2,])
    }
    
    # write merged data to the outDir
    WritePTQobs(x = obs1,paste(outDir,obsFiles[i],sep="/"))
  }
  return(paste(outDir,obsFiles,sep="/"))
}


## -------------------------------------------------------------------------------
## prepare model forcing data and initial conditions
##
## the function returns a list variable with a status code (T/F) and a issuedate
## possibly modified according to the available data on the data source
##
getModelForcing<-function(appSetup,appInput,dataSource="local",hindcast=T){
  
  ## Historical simulation
  if(appSetup$appName=="historical"){  
    ##
    ## A. Simulation and Warmup period (BDATE; CDATE; EDATE)
    ##
    ## CDATE: simulation start date, from application input
    cdate.Str = appInput$cdate
    cdate.Num = as.POSIXct(cdate.Str,tz = "GMT")
    
    # limit cdate to the available state_save files
    cdate.Num = max(c(cdate.Num,appSetup$bdateMin))
    cdate.Str = as.character(cdate.Num)
    
    ## BDATE: start of warmup period
    ##
    if(appInput$assimOn=="on"){
      ##   Assimilation run
      ## 
      ##   make sure BDATE is before the start of the previous rainy season
      ##   which is 1/1 this year if MONTH>6 otherwise 1/1 previous year
      ##
      ##   needed in Niger-River in order to get good ensemble spread
      ##   (later we must can state ensemble inflation to speed up initialization)
      ##
      ##   in any case, it's probabky good to start at the start of hydrological year.
      ##
      cdata.month = as.integer(substr(cdate.Str,start = 6,stop = 7))
      if(cdata.month>6){
        bdate.Str = paste(substr(cdate.Str,1,4),"-01-01",sep="")
        bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
      }else{
        cdata.year = as.integer(substr(cdate.Str,start = 1,stop = 4))
        bdata.year = cdata.year -1
        bdate.Str = paste(as.character(bdata.year),"-01-01",sep="")
        bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
      }
    }else{
      ## Normal simulation (no assimilation)
      ##
      ## start of warmup period, 1 Jan same year as cdate
      bdate.Str = paste(substr(cdate.Str,1,4),"-01-01",sep="")
      bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
    }  
    
    # limit bdate to available state_save files
    bdate.Num = min(c(bdate.Num,appSetup$bdateMax))
    bdate.Num = max(c(bdate.Num,appSetup$bdateMin))
    bdate.Str = as.character(bdate.Num)
    
    ## EDATE: simulation end date, from application input (might change depending on data availability)
    edate.Str = appInput$edate
    edate.Num = as.POSIXct(edate.Str,tz = "GMT")
    
    ##
    ## B. Forcing data (GFD archive + hindcast data from catalogue if needed)
    ##
    
    # B.1 check if Archive and Hindcasts are needed
    useArchive  = F ; if(bdate.Num<=forcing.archive.end){useArchive  = T}
    useHindcast = F ; if(edate.Num>forcing.archive.end){useHindcast  = T}
    
    # B.2 Archive data - if needed
    if(useArchive){
      # log message
      if(app.sys=="tep"){rciop.log ("DEBUG", " useArchive = TRUE ...", "/util/R/hypeapps-utils.R")}
      
      # create tmpDir/forcing/archive
      archiveDir = paste(appSetup$tmpDir,"/forcing/archive",sep="")
      dir.create(archiveDir,recursive = T,showWarnings = F)
      
      # local archive exists (or archive at store.terradue.com/dgustafsson/model ...)
      if(appSetup$forcingArchiveExist){
        # copy forcing files from local archive
        for(i in 1:length(forcing.files)){
          rciop.copy(url = paste(appSetup$forcingArchiveURL,forcing.files[i],sep="/"),
                     target = archiveDir)
        }
        archiveFound = T
      }else{
        # download archive from data catalogue
        sysCmd=paste("curl -o ", appSetup$tmpDir,"/archive.zip https://store.terradue.com//smhi/gfd/niger-hype/hindcast/files/v1/archive.zip",sep="")
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying command >> ",sysCmd,sep=""), "/util/R/hypeapps-utils.R")}
        a=system(sysCmd,intern=T)
        # unzip to forcing/archive
        archiveFile=paste(appSetup$tmpDir,"archive.zip",sep="/")
        if(file.exists(archiveFile)){
          if(app.sys=="tep"){rciop.log ("DEBUG", paste("archiveFile from https://catalogue.terradue.com/hydro-smhi/ = ",archiveFile,sep=""), "/util/R/hypeapps-utils.R")}
          unzip(zipfile = archiveFile,overwrite = T,exdir = archiveDir)
          archiveFound = T
        }else{
          return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=T))
          archiveFound = F
        }
      }
    }
    
    # B.3 Hindcast data - if needed
    if(useHindcast){
      if(app.sys=="tep"){rciop.log ("DEBUG", " useHindcast = TRUE ...", "/util/R/hypeapps-utils.R")}
      # read daily hindcasts data from repository to reach the edate 
      # with as updated data as possible
      
      # first check the issue date of the latest available hindcast
      today.Str  = Sys.Date()
      today.Num   = as.POSIXct(today.Str,tz = "GMT")
      issueDate.Num = today.Num
      
      # try downloading hindcast with issue date equal to the requested
      downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                      modelName=appSetup$modelName,
                                      productName="hindcast",
                                      tmpDir=appSetup$tmpDir)
      if(!downloadInfo$status){
        # try download data issued issueDateNum-1:9
        i=1
        while(i<20 & !downloadInfo$status){
          # search issue data backwards
          issueDate.Num = today.Num-i*86400
          downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                          modelName=appSetup$modelName,
                                          productName="hindcast",
                                          tmpDir=appSetup$tmpDir)
          i=i+1
        }
      }
      # if data downloaded, unzip
      if(downloadInfo$status){
        # save latest issue date
        lastIssueDate.Num=issueDate.Num
        
        #unzip downloaded file
        hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(lastIssueDate.Num),sep="")
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" hindcastDir = ", hindcastDir,sep=""), "/util/R/hypeapps-utils.R")}
        
        dir.create(hindcastDir,recursive = T,showWarnings = F)
        unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = hindcastDir)
        
        # check first and last date in file
        firstDate.last=substr(system(paste("gawk 'NR==2{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
        lastDate.last=substr(system(paste("gawk 'END{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
        
        # adjust edate to laste date in file, if originally later
        edate.Num = min(edate.Num,as.POSIXct(lastDate.last,tz="GMT"))
        
        # download files covering 2017-01-01 until edate or first day in last file
        # is within the first month of the file
        notFound = T
        y1=2017
        m1=1
        d1=1
        
        nIssueDates=0
        
        while(notFound){
          issueDate.Num = ymd2posix(y1,m1,d1)+123*86400
          # try downloading hindcast with issue date equal to the requested
          downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                          modelName=appSetup$modelName,
                                          productName="hindcast",
                                          tmpDir=appSetup$tmpDir)
          if(downloadInfo$status){
            # save issueDate to vector of successful dates
            if(nIssueDates==0){
              issueDates = issueDate.Num
            }else{
              issueDates=c(issueDates,issueDate.Num)
            }
            nIssueDates = nIssueDates+1
            if(app.sys=="tep"){rciop.log ("DEBUG", paste(" issueDates = ", issueDates,sep=""), "/util/R/hypeapps-utils.R")}
            
            #unzip downloaded file
            hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(issueDate.Num),sep="")
            dir.create(hindcastDir,recursive = T,showWarnings = F)
            unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = hindcastDir)
            
            # check first and last date in file
            firstDate=substr(system(paste("gawk 'NR==2{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
            lastDate=substr(system(paste("gawk 'END{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
            
            # next month as integer
            y2=as.integer(substr(firstDate,1,4))
            m2=as.integer(substr(firstDate,6,7))+1
            d2=1
            if(m2>12){
              y2=y2+1
              m2=1
            }
            # next month as posix
            nextMonth=ymd2posix(y2,m2,d2)
            if(edate.Num<=nextMonth){
              notFound=F
              useLast=F
            }else{
              notFound=T
              useLast=F
              y1=y2
              m1=m2
              d1=d2
            }
          }else{
            notFound=F
            useLast=T
          }
          if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... download loop ... ",as.character(issueDate.Num),sep=" "), 
                                        "/util/R/hypeapps-utils.R")}
        }
        # add last file if needed
        if(useLast){
          issueDates=c(issueDates,lastIssueDate.Num)
        }
        if(app.sys=="tep"){rciop.log ("DEBUG", paste(" issueDates = ", issueDates,sep=","), "/util/R/hypeapps-utils.R")}
        hindcastFound = T
        if(app.sys=="tep"){rciop.log ("DEBUG", " hindcast download ready!", "/util/R/hypeapps-utils.R")}
      }else{
        if(app.sys=="tep"){rciop.log ("DEBUG", " hindcast download failed!", "/util/R/hypeapps-utils.R")}
        hindcastFound = F
      }
    }
    
    # B.4 Merge obsfiles and copy to runDir
    if(app.sys=="tep"){rciop.log ("DEBUG", " Starting to merge/copy obs-files ...", "/util/R/hypeapps-utils.R")}
    obsFiles=c("Pobs.txt","Tobs.txt","TMAXobs.txt","TMINobs.txt")
    for(i in 1:length(obsFiles)){
      if(useArchive){
        if(app.sys=="tep"){rciop.log ("DEBUG", " ... trying archive files ...", "/util/R/hypeapps-utils.R")}
        if(archiveFound){
          if(useHindcast){
            if(hindcastFound){
              if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... reading archive ",obsFiles[i],sep=""), "/util/R/hypeapps-utils.R")}
              obsData = ReadPTQobs(paste(archiveDir,obsFiles[i],sep="/"))
              iStart  = which(obsData[,1]<=as.POSIXct(bdate.Num,tz="GMT"))
              iStart  = iStart[length(iStart)]
              iEnd    = which(obsData[,1]>=as.POSIXct("2017-01-01",tz="GMT"))
              iEnd    = iEnd[1]-1
              obsData = obsData[iStart:iEnd,]
            }else{
              file.copy(from = paste(archiveDir,obsFiles[i],sep="/"),to = paste(appSetup$runDir,obsFiles[i],sep="/"))
            }
          }else{
            file.copy(from = paste(archiveDir,obsFiles[i],sep="/"),to = paste(appSetup$runDir,obsFiles[i],sep="/"))
          }
        }
      }
      if(useHindcast){
        if(app.sys=="tep"){rciop.log ("DEBUG", " ... trying hindcast files ...", "/util/R/hypeapps-utils.R")}
        if(hindcastFound){
          if(app.sys=="tep"){rciop.log ("DEBUG", paste(" issueDates = ", issueDates,sep=","), "/util/R/hypeapps-utils.R")}
          for(j in 1:length(issueDates)){
            if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... issueDates[j] ...",issueDates[j],sep=""), "/util/R/hypeapps-utils.R")}
            hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcast/",DATE2INFODATE(issueDates[j]),sep="")
            newData = ReadPTQobs(paste(hindcastDir,obsFiles[i],sep="/"))
            if(!useArchive & j==1){
              obsData=newData
            }else{
              iOld = which(obsData[,1]<newData[1,1])
              obsData = rbind(obsData[iOld,],newData)
            }
          }
          WritePTQobs(x = obsData,paste(appSetup$runDir,obsFiles[i],sep="/"))
        }
      }
    }
    # finally, check bdate, cdate and edate versus content of Pobs.txt in runDir and available state_save files
    firstDate.last=substr(system(paste("gawk 'NR==2{print}'",paste(appSetup$runDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
    lastDate.last=substr(system(paste("gawk 'END{print}'",paste(appSetup$runDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
    
    # bdate - no earlier than the first occurance of 1 Jan at or after first date in file
    bdate.Num = max(bdate.Num,as.POSIXct(firstDate.last,tz="GMT"))
    nextStateDate = which(appSetup$stateDates>=bdate.Num)
    if(length(nextStateDate)>0){
      bdate.Num = appSetup$stateDates[nextStateDate[1]]
    }
    
    # cdate, no earlier than the adjusted bdate
    cdate.Str = appInput$cdate
    cdate.Num = as.POSIXct(cdate.Str,tz = "GMT")
    cdate.Num = max(cdate.Num,bdate.Num)
    
    # edate - no later than last date in obs-file, and not earlier than bdate and cdate
    edate.Num = min(edate.Num,as.POSIXct(lastDate.last,tz="GMT"))
    
    dateError=F
    if(edate.Num<cdate.Num){dateError=T}
    if(edate.Num<bdate.Num){dateError=T}
    if(cdate.Num<bdate.Num){dateError=T}
    
    ##
    ## C. Initial conditions matching the bdate
    ##
    iState = which(appSetup$stateDates==bdate.Num)
    if(length(iState)>0){
      rciop.copy(url = paste(appSetup$stateFilesURL,appSetup$stateFiles[iState],sep="/"), 
                 target = appSetup$runDir)
    }else{
      dateError=T
    }
    
    ##
    ## D. Return list with updated bdate, cdate, and edate
    ##
    if(!dateError){
      return(list("status"=T,
                  "localFile"=paste(appSetup$runDir,"Pobs.txt",sep="/"),
                  "issueDate"=NULL,
                  "archive"=T,
                  "bdate"=bdate.Num,
                  "cdate"=cdate.Num,
                  "edate"=edate.Num,
                  "outstateDate"=NULL,
                  "stateFile"=paste(appSetup$runDir,appSetup$stateFiles[iState],sep="/")))
    }else{
      return(list("status"=F,
                  "localFile"=NULL,
                  "issueDate"=NULL,
                  "archive"=F,
                  "bdate"=NULL,
                  "cdate"=NULL,
                  "edate"=NULL,
                  "outstateDate"=NULL,
                  "stateFile"=NULL))
    }
    
  }else if(appSetup$appName=="forecast"){
    
    ### Forcing data for forecast initialization (hindcast=T) or forecast simulations (hindcast=F)
    
    # set productName
    if(hindcast){
      prod.name="hindcast"
    }else{
      prod.name="forecast"
    }
    
    # requested forecast issue date from user input
    issueDate.Str = appInput$idate
    issueDate.Num = as.POSIXct(appInput$idate,tz = "GMT")
    
    # try downloading GFD forecast/hindcast with issue date equal to the requested
    downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                    modelName=appSetup$modelName,
                                    productName=prod.name,
                                    tmpDir=appSetup$tmpDir)
    if(!downloadInfo$status){
      # try download data issued issueDateNum-1:9
      i=1
      while(i<10 & !downloadInfo$status){
        # search issue data backwards
        issueDate.Num = as.POSIXct(appInput$idate,tz = "GMT")-i*86400
        downloadInfo = getGFDzipFromTep(issueDateNum = issueDate.Num,
                                        modelName=appSetup$modelName,
                                        productName=prod.name,
                                        tmpDir=appSetup$tmpDir)
        i=i+1
      }
    }
    # unzip the downloaded file, or return error
    if(downloadInfo$status){
      
      # write download status to log file
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... download status = ",as.character(downloadInfo$status),sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... remote file = ",downloadInfo$remoteFile,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... local file = ",downloadInfo$localFile,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      
      #unzip downloaded file to temporary directory (hindcast) or runDir (forecast)
      if(hindcast){
        hindcastDir = paste(appSetup$tmpDir,"/forcing/hindcastTemp",sep="")
        dir.create(hindcastDir,recursive = T,showWarnings = F)
        unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = hindcastDir)
        
        # check first and last date in the obsfile:
        # a. get first and last date from file in text string
        firstDate=substr(system(paste("gawk 'NR==2{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
        lastDate=substr(system(paste("gawk 'END{print}'",paste(hindcastDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
      }else{
        unzip(zipfile = downloadInfo$localFile,overwrite = T,exdir = appSetup$runDir)
        firstDate=substr(system(paste("gawk 'NR==2{print}'",paste(appSetup$runDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
        lastDate=substr(system(paste("gawk 'END{print}'",paste(appSetup$runDir,"Pobs.txt",sep="/"),sep=" "),intern=T),1,10)
      }
      
      # b. write to log
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... prod.name = ",prod.name,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... Pobs, first date = ",firstDate,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... Pobs, last date = ",lastDate,sep=" "), 
                                    "/util/R/hypeapps-utils.R")}
      # c. transform text date "yyyy-mm-dd" to Posix date format
      firstDate=as.POSIXct(firstDate,tz = "GMT")
      lastDate=as.POSIXct(lastDate,tz = "GMT")
      
      # adapt bdate, cdate, edate, and outstatedate to forecast issue date and available dates in the file
      if(hindcast){
        cdate=firstDate
        edate=lastDate
        issueDate.Num=lastDate+86400
        outstateDate=issueDate.Num
        
        # set bdate to nearest available state_save file (1st guess, 1 Jan in same year as cdate)
        cdate.Num = cdate
        cdate.Str = as.character(cdate.Num)
        
        # BDATE: start of warmup period - specific if assimilation or not
        if(appInput$assimOn=="on"){
          ##   Assimilation run
          ## 
          ##   make sure BDATE is before the start of the previous rainy season
          ##   which is 1/1 this year if MONTH>6 otherwise 1/1 previous year
          ##
          ##   needed in Niger-River in order to get good ensemble spread
          ##   (later we must can state ensemble inflation to speed up initialization)
          ##
          ##   in any case, it's probably good to start at the previous start of a hydrological year.
          ##
          cdata.month = as.integer(substr(cdate.Str,start = 6,stop = 7))
          if(cdata.month>6){
            bdate.Str = paste(substr(cdate.Str,1,4),"-01-01",sep="")
            bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
          }else{
            cdata.year = as.integer(substr(cdate.Str,start = 1,stop = 4))
            bdata.year = cdata.year -1
            bdate.Str = paste(as.character(bdata.year),"-01-01",sep="")
            bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
          }
        }else{
          ## Normal simulation (no assimilation)
          ##
          ## start of warmup period, 1 Jan same year as cdate
          bdate.Str = paste(substr(cdate.Str,1,4),"-01-01",sep="")
          bdate.Num = as.POSIXct(bdate.Str,tz = "GMT")
        }  
        
        bdate.Num = min(c(bdate.Num,appSetup$bdateMax))
        bdate.Num = max(c(bdate.Num,appSetup$bdateMin))
        bdate.Str = as.character(bdate.Num)
        bdate = bdate.Num
        ##
        ## pad the hindcast forcing files for the period bdate:cdate
        ## using archive and/or additional hindcast data
        ##
        if(bdate.Num<cdate.Num){
          
          # pad with archive dna hindcast data bdate:cdate into hindcastTemp2
          hindcastDir2 = paste(appSetup$tmpDir,"/forcing/hindcastTemp2",sep="")
          dir.create(hindcastDir2,recursive = T,showWarnings = F)
          getHindcastForcingData(bdate.Num,cdate.Num,appSetup,obsFiles=forcing.files,hindcastDir2,useRdata=F)
          
          # merge with the current hindcast data into runDir            
          mergeObsFiles(hindcastDir,hindcastDir2,appSetup$runDir,bdate.Num,edate,obsFiles=forcing.files)
        }
      }else{
        bdate=firstDate # if this is the forecast run, then the bdate=cdate=issuedate=first date in file
        cdate=firstDate
        edate=lastDate
        issueDate.Num=firstDate+86400
        outstateDate=NA
      }
      
      # also make sure the initial statefile is available, either from statefile archive, or from hindcast result folder
      stateFile=paste(appSetup$runDir,"/state_save",DATE2INFODATE(bdate),".txt",sep="")
      if(hindcast){
        iState = which(appSetup$stateDates==bdate.Num)
        if(length(iState)>0){
          rciop.copy(url = paste(appSetup$stateFilesURL,appSetup$stateFiles[iState],sep="/"), 
                     target = appSetup$runDir)
          stateFile = paste(appSetup$runDir,appSetup$stateFiles[iState],sep="/")
          stateError=F
        }else{
          stateError=T
          stateFile=NULL
        }
        if(!stateError){
          if(!file.exists(stateFile)){
            if(app.sys=="tep"){rciop.log ("DEBUG", paste(" statefile missing in inputdata = ",stateFile,sep=" "),
                                          "/util/R/hypeapps-utils.R")}
            return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=F,"stateFile"=NULL))
          }
        }else{
          if(app.sys=="tep"){rciop.log ("DEBUG", " statefile missing in inputdata = stateError",
                                        "/util/R/hypeapps-utils.R")}
          return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=F,"stateFile"=NULL))
        }
      }else{
        #copy statefile from previous hindcast resultdir
        stateFileFromHindcast=paste(appSetup$resDir[1],"/state_save",DATE2INFODATE(bdate),".txt",sep="")
        if(file.exists(stateFileFromHindcast)){
          file.copy(from=stateFileFromHindcast,to=stateFile,overwrite = T)
        }else{
          if(app.sys=="tep"){rciop.log ("DEBUG", paste(" statefile missing in hindcast resultdir = ",stateFileFromHindcast,sep=" "), "/util/R/hypeapps-utils.R")}
          return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=F))
        }
      }
      return(list("status"=T,
                  "localFile"=downloadInfo$localFile,
                  "issueDate"=issueDate.Num,
                  "archive"=F,
                  "bdate"=bdate,
                  "cdate"=cdate,
                  "edate"=edate,
                  "outstateDate"=outstateDate,
                  "stateFile"=stateFile))
    }else{
      return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=F))
    }       
  }else{
    return(list("status"=F,"localFile"=NULL,"issueDate"=NA,"archive"=F))
  }       
}

## -------------------------------------------------------------------------------
## modify some model input files
updateModelInput<-function(appSetup=NULL,appInput=NULL,hindcast=NULL,modelForcing=NULL,xobsInput=NULL){
  
  if(appSetup$appName=="historical"){
    
    # TD. check start and end data versus available dates in forcing data (should be given in modelForcing input list)
    
    ## INFO.TXT ##
    # read info.txt template
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("appSetup$runDir = ",appSetup$runDir,sep=""), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("TMPDIR=",appSetup$tmpDir,sep=""), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("dir(TMPDIR)=",dir(appSetup$tmpDir),sep=" "), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("dir(appSetup$runDir)=",dir(appSetup$runDir),sep=" "), "/util/R/hypeapps-utils.R")}
    
    #select info-file for assimilation
    if(appInput$assimOn=="on"){
      info=readInfo(paste(appSetup$runDir,"info-historical-assimilation.txt",sep="/"))
    }else{
      info=readInfo(paste(appSetup$runDir,"info-historical.txt",sep="/"))
    }
    
    # modify info according to inputs
    # resultdir
    if(info$isResDir){
      info$info.lines[info$resultdir.lineNr]=paste(paste('resultdir', appSetup$resDir, sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('resultdir', appSetup$resDir, sep=" "),"/",sep=""))
    }
    # modeldir
    if(info$isModDir){
      info$info.lines[info$modeldir.lineNr]=paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep=""))
    }
    # bdate,cdate,edate
    info$info.lines[info$bdate.lineNr]=paste('bdate',as.character(modelForcing$bdate),sep=" ")
    info$info.lines[info$cdate.lineNr]=paste('cdate',as.character(modelForcing$cdate),sep=" ")
    info$info.lines[info$edate.lineNr]=paste('edate',as.character(modelForcing$edate),sep=" ")
    
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("  bdate=",as.character(modelForcing$bdate),sep=" "), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("  cdate=",as.character(modelForcing$cdate),sep=" "), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("  edate=",as.character(modelForcing$edate),sep=" "), "/util/R/hypeapps-utils.R")}
    
    
    # output variables (from appInput, xobsInput, and assimVar)
    outputVariables = strsplit(appInput$outvars,split = ",")[[1]] # from appInput
    if(!is.null(xobsInput)){                                      # from xobsInput
      for(i in 1:length(xobsInput$xobsVar)){
        if(!is.na(xobsInput$xobsVar[i])){
          outputVariables = c(outputVariables,xobsInput$xobsVar[i])
        }
      }
    }
    if(appInput$assimOn=="on"){                                   # from assimVar
      for(i in 1:length(appInput$assimVar)){
        if(appInput$assimVar[i]!="9999,9999" & 
           appInput$assimVar[i]!="OPEN,LOOP" & 
           nchar(appInput$assimVar[i])==9){
          assimVariables=strsplit(appInput$assimVar[i],split = ",")
          outputVariables = c(outputVariables,assimVariables)
        }
      }
    }
    outputVariables = tolower(outputVariables)
    outputVariables = unique(outputVariables)
    if(length(outputVariables)>1){
      outVariables=outputVariables[1]
      for(i in 2:length(outputVariables)){
        outVariables=paste(outVariables,outputVariables[i],sep=",")
      }
    }else{
      outVariables=outputVariables
    }
    
    # output basins (from appInput and from xobsInput)
    outputBasins = as.integer(strsplit(appInput$outbasins,split=",")[[1]])
    if(!is.null(xobsInput)){
      for(i in 1:length(xobsInput$xobsSubid)){
        if(!is.na(xobsInput$xobsSubid[i])){
          outputBasins = unique(outputBasins,xobsInput$xobsSubid[i])
        }        
      }
    }
    if(length(outputBasins)>1){
      outBasins=as.character(outputBasins[1])
      for(i in 2:length(outputBasins)){
        outBasins=paste(outBasins,as.character(outputBasins[i]),sep=",")
      }
    }else{
      outBasins=as.character(outputBasins)
    }
    
    # basinoutput
    info$info.lines[info$basinoutput_variable.lineNr]=paste('basinoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    info$info.lines[info$basinoutput_subbasin.lineNr]=paste('basinoutput subbasin',gsub(pattern=",",replacement = " ",outBasins),sep=" ")
    
    # timeoutput
    info$info.lines[info$timeoutput_variable.lineNr]=paste('timeoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    
    # mapoutput
    info$info.lines[info$mapoutput_variable.lineNr]=paste('mapoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    
    # write info file
    writeInfoRes = writeInfo(info$info.lines,filenm = paste(appSetup$runDir,"info.txt",sep="/"))
    
    # select AssimInfo.txt for assimilation simulation
    if(appInput$assimOn=="on"){
      # copy AssimInfo depending on assimilation variables
      if(appInput$assimVar=="AOWL,WCOM"|appInput$assimVar=="aowl,wcom"){
        file.copy(from = paste(appSetup$runDir,"AssimInfo-AOWL.txt",sep="/"),to = paste(appSetup$runDir,"AssimInfo.txt",sep="/"))
        if(app.sys=="tep"){rciop.log ("DEBUG", paste("  assimilation run:","AOWL vs WCOM",sep=" "), "/util/R/hypeapps-utils.R")}
        
      }
      if(appInput$assimVar=="OPEN,LOOP"|appInput$assimVar=="open,loop"){
        file.copy(from = paste(appSetup$runDir,"AssimInfo-Openloop.txt",sep="/"),to = paste(appSetup$runDir,"AssimInfo.txt",sep="/"))
        if(app.sys=="tep"){rciop.log ("DEBUG", paste("  assimilation run:","open-loop",sep=" "), "/util/R/hypeapps-utils.R")}
      }
    }
    
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("dir(appSetup$runDir)=",dir(appSetup$runDir),sep=" "), "/util/R/hypeapps-utils.R")}
    
    # return updated list of output variables and subbasins
    modelInput = list("outvars"=outVariables,"outbasins"=outBasins)
    
    return(modelInput)
    
    
  }else if(appSetup$appName=="forecast"){
    
    # read template info for hindcast or forecast simulation
    if(hindcast){
      info=readInfo(paste(appSetup$runDir,"info-hindcast.txt",sep="/"))
      dirNum=1
    }else{
      info=readInfo(paste(appSetup$runDir,"info-forecast.txt",sep="/"))
      dirNum=2
    }
    
    # update resultdir
    if(info$isResDir){
      info$info.lines[info$resultdir.lineNr]=paste(paste('resultdir', appSetup$resDir[dirNum], sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('resultdir', appSetup$resDir[dirNum], sep=" "),"/",sep=""))
    }
    
    # update modeldir
    if(info$isModDir){
      info$info.lines[info$modeldir.lineNr]=paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep=""))
    }
    
    # bdate,cdate,edate
    info$info.lines[info$bdate.lineNr]=paste('bdate',DATE2INFODATE(modelForcing$bdate),sep=" ")
    info$info.lines[info$cdate.lineNr]=paste('cdate',DATE2INFODATE(modelForcing$cdate),sep=" ")
    info$info.lines[info$edate.lineNr]=paste('edate',DATE2INFODATE(modelForcing$edate),sep=" ")
    
    # output variables (from appInput, xobsInput, and assimVar)
    outputVariables = strsplit(appInput$outvars,split = ",")[[1]] # from appInput
    if(!is.null(xobsInput)){                                      # from xobsInput
      for(i in 1:length(xobsInput$xobsVar)){
        if(!is.na(xobsInput$xobsVar[i])){
          outputVariables = c(outputVariables,xobsInput$xobsVar[i])
        }
      }
    }
    if(appInput$assimOn=="on"){                                   # from assimVar
      for(i in 1:length(appInput$assimVar)){
        if(appInput$assimVar[i]!="9999,9999" & 
           appInput$assimVar[i]!="OPEN,LOOP" & 
           nchar(appInput$assimVar[i])==9){
          assimVariables=strsplit(appInput$assimVar[i],split = ",")
          outputVariables = c(outputVariables,assimVariables)
        }
      }
    }
    outputVariables = tolower(outputVariables)
    outputVariables = unique(outputVariables)
    if(length(outputVariables)>1){
      outVariables=outputVariables[1]
      for(i in 2:length(outputVariables)){
        outVariables=paste(outVariables,outputVariables[i],sep=",")
      }
    }else{
      outVariables=outputVariables
    }
    
    # output basins (from appInput and from xobsInput)
    outputBasins = as.integer(strsplit(appInput$outbasins,split=",")[[1]])
    if(!is.null(xobsInput)){
      for(i in 1:length(xobsInput$xobsSubid)){
        if(!is.na(xobsInput$xobsSubid[i])){
          outputBasins = unique(outputBasins,xobsInput$xobsSubid[i])
        }        
      }
    }
    if(length(outputBasins)>1){
      outBasins=as.character(outputBasins[1])
      for(i in 2:length(outputBasins)){
        outBasins=paste(outBasins,as.character(outputBasins[i]),sep=",")
      }
    }else{
      outBasins=as.character(outputBasins)
    }
    
    
    # basinoutput
    info$info.lines[info$basinoutput_variable.lineNr]=paste('basinoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    info$info.lines[info$basinoutput_subbasin.lineNr]=paste('basinoutput subbasin',gsub(pattern=",",replacement = " ",outBasins),sep=" ")
    
    # timeoutput
    info$info.lines[info$timeoutput_variable.lineNr]=paste('timeoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    
    # mapoutput
    info$info.lines[info$mapoutput_variable.lineNr]=paste('mapoutput variable',gsub(pattern=",",replacement = " ",outVariables),sep=" ")
    
    # outstatedate
    if(hindcast){
      info$info.lines[info$outstatedate.lineNr]=paste('outstatedate',DATE2INFODATE(modelForcing$issueDate),sep=" ")
    }
    
    # remove existing XobsFile if existing
    if(is.null(xobsInput) & file.exists(paste(appSetup$runDir,"Xobs.txt",sep="/"))){
      file.remove(paste(appSetup$runDir,"Xobs.txt",sep="/"))
    }
    
    # write info file
    return(writeInfo(info$info.lines,filenm = paste(appSetup$runDir,"info.txt",sep="/")))
    
  }
}

## -------------------------------------------------------------------------------
## prepare application outputs
prepareHypeAppsOutput<-function(appSetup=NULL,appInput=NULL,modelInput=NULL,modelForcing = NULL,runRes=NULL,appDate=NULL){
  
  # Create folder for data to be published
  outDir = paste(appSetup$tmpDir,'output',sep="/")
  dir.create(outDir,recursive = T,showWarnings = F)
  
  ## output file prefixes, to order the results better
  prefix.img =paste("001","_",appDate,sep="")
  prefix.csv =paste("002","_",appDate,sep="")
  prefix.bas =paste("003","_",appDate,sep="")
  prefix.map =paste("004","_",appDate,sep="")
  prefix.tim =paste("005","_",appDate,sep="")
  prefix.oth =paste("006","_",appDate,sep="")
  prefix.log =paste("000","_",appDate,sep="")
  prefix.wl.txt = paste("004","_",appDate,sep="")
  prefix.wl.png = paste("001","_",appDate,sep="")
  
  ## get hype2csv file from its URL
  if(!is.null(appSetup$hype2csvURL)){
    rciop.copy(url = appSetup$hype2csvURL, target = appSetup$tmpDir)
    if(file.exists(paste(appSetup$tmpDir,appSetup$hype2csvFile,sep="/"))){
      hype2csvFile = paste(appSetup$tmpDir,appSetup$hype2csvFile,sep="/")
      hype2csvExists=T
    }else{
      hype2csvFile = "9999"
      hype2csvExists=F
    }
  }else{
    hype2csvFile = "9999"
    hype2csvExists=F    
  }
  
  ## Post-process requested outputs (copy some files...)
  if(appSetup$appName=="historical"){
    
    # if runRes=1 (error)
    if(runRes==1){
      # if existing, copy log-file in model working folder
      hyssLogFile = dir(path = appSetup$runDir , pattern =".log")
      if(length(hyssLogFile)>0){
        if(app.sys=="tep"){
          #          res <- rciop.copy(paste(appSetup$runDir,hyssLogFile[1],sep="/"), outDir, uncompress=TRUE)
          res <- file.copy(from = paste(appSetup$runDir,hyssLogFile[1],sep="/"), 
                           to = paste(outDir,paste(prefix.log,hyssLogFile[1],sep="_"),sep="/"), 
                           overwrite = T)
        }
      }
      
    }else{
      
      # copy log-file in model working folder
      hyssLogFile = dir(path = appSetup$runDir , pattern =".log")
      if(length(hyssLogFile)>0){
        if(app.sys=="tep"){
          #          res <- rciop.copy(paste(appSetup$runDir,hyssLogFile[1],sep="/"), outDir, uncompress=TRUE)
          res <- file.copy(from = paste(appSetup$runDir,hyssLogFile[1],sep="/"), 
                           to = paste(outDir,paste(prefix.log,hyssLogFile[1],sep="_"),sep="/"), 
                           overwrite = T)
        }
      }
      
      # list files in model results folder
      timeFiles   = dir(path = appSetup$resDir , pattern ="time")
      mapFiles    = dir(path = appSetup$resDir , pattern ="map")
      subassFiles = dir(path = appSetup$resDir , pattern ="subass")
      simassFile  = dir(path = appSetup$resDir , pattern ="simass")
      allFiles    = dir(path = appSetup$resDir , pattern =".txt")
      
      # copy time files
      if(length(timeFiles)>0){
        for(i in 1:length(timeFiles)){
          if(app.sys=="tep"){
            #            rciop.copy(paste(appSetup$resDir,timeFiles[i],sep="/"), outDir, uncompress=TRUE)
            if(appInput$assimOn=="off"){
              file.copy(from = paste(appSetup$resDir,timeFiles[i],sep="/"),
                        to = paste(outDir,paste(prefix.tim,timeFiles[i],sep="_"),sep="/"),
                        overwrite = T)
            }else{
              # check if file extension indicates "max" or "min"
              nchars = nchar(timeFiles[i])
              last4  = substr(timeFiles[i],nchars-7,nchars-4)
              if(last4=="_002"){
                fileOUT = paste(substr(timeFiles[i],1,nchars-7),"min.txt",sep="")
              }else if(last4=="_003"){
                fileOUT = paste(substr(timeFiles[i],1,nchars-7),"max.txt",sep="")
              }else{
                fileOUT = paste(substr(timeFiles[i],1,nchars-4),"_mean.txt",sep="")
              }
              file.copy(from = paste(appSetup$resDir,timeFiles[i],sep="/"),
                        to = paste(outDir,paste(prefix.tim,fileOUT,sep="_"),sep="/"),
                        overwrite = T)
            }
          }
        }
      }
      
      # copy map files
      if(length(mapFiles)>0){
        for(i in 1:length(mapFiles)){
          if(app.sys=="tep"){
            #            rciop.copy(paste(appSetup$resDir,mapFiles[i],sep="/"), outDir, uncompress=TRUE)
            if(appInput$assimOn=="off"){
              file.copy(from = paste(appSetup$resDir,mapFiles[i],sep="/"),
                        to = paste(outDir,paste(prefix.map,mapFiles[i],sep="_"),sep="/"),
                        overwrite = T)
            }else{
              # check if file extension indicates "max" or "min"
              nchars = nchar(mapFiles[i])
              last4  = substr(mapFiles[i],nchars-7,nchars-4)
              if(last4=="_002"){
                fileOUT = paste(substr(mapFiles[i],1,nchars-7),"min.txt",sep="")
              }else if(last4=="_003"){
                fileOUT = paste(substr(mapFiles[i],1,nchars-7),"max.txt",sep="")
              }else{
                fileOUT = paste(substr(mapFiles[i],1,nchars-4),"_mean.txt",sep="")
              }
              file.copy(from = paste(appSetup$resDir,mapFiles[i],sep="/"),
                        to = paste(outDir,paste(prefix.map,fileOUT,sep="_"),sep="/"),
                        overwrite = T)
            }
            
          }
        }
      }
      # copy subass files
      if(length(subassFiles)>0){
        for(i in 1:length(subassFiles)){
          if(app.sys=="tep"){
            #            rciop.copy(paste(appSetup$resDir,subassFiles[i],sep="/"), outDir, uncompress=TRUE)
            file.copy(from = paste(appSetup$resDir,subassFiles[i],sep="/"),
                      to = paste(outDir,paste(prefix.oth,subassFiles[i],sep="_"),sep="/"),
                      overwrite = T)
          }
        }
      }
      # copy simass files
      if(length(simassFile)>0){
        if(app.sys=="tep"){
          #          rciop.copy(paste(appSetup$resDir,simassFile[1],sep="/"), outDir, uncompress=TRUE)
          file.copy(from = paste(appSetup$resDir,simassFile[1],sep="/"),
                    to = paste(outDir,paste(prefix.oth,simassFile[1],sep="_"),sep="/"),
                    overwrite = T)
        }
      }
      
      # basin outputfiles
      outbasins = strsplit(modelInput$outbasins,split = ",")[[1]]
      zeroString="0000000000000000000000000000000"
      basinFiles=NULL
      if(length(outbasins)>0){
        for(i in 1:length(outbasins)){
          outFile=paste(outbasins[i],".txt",sep="")
          ni=nchar(outFile)
          for(j in 1:length(allFiles)){
            nj=nchar(allFiles[j])
            if(nj>=ni){
              if(substr(allFiles[j],nj-ni+1,nj)==outFile){
                if(nj>ni){
                  if(substr(allFiles[j],1,nj-ni)==substr(zeroString,1,nj-ni)){
                    if(app.sys=="tep"){
                      #                      rciop.copy(paste(appSetup$resDir,allFiles[j],sep="/"), outDir, uncompress=TRUE)
                      if(appInput$assimOn=="off"){
                        # Normal run
                        file.copy(from = paste(appSetup$resDir,allFiles[j],sep="/"),
                                  to = paste(outDir,paste(prefix.bas,allFiles[j],sep="_"),sep="/"),
                                  overwrite = T)
                      }else{
                        # Assimilation run
                        nchars = nchar(allFiles[j])
                        # min (002)
                        fileIN = paste(substr(allFiles[j],1,nchars-4),"_002.txt",sep="")
                        fileOUT = paste(substr(allFiles[j],1,nchars-4),"_min.txt",sep="")
                        file.copy(from = paste(appSetup$resDir,fileIN,sep="/"),
                                  to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                  overwrite = T)
                        # max (003)
                        fileIN = paste(substr(allFiles[j],1,nchars-4),"_003.txt",sep="")
                        fileOUT = paste(substr(allFiles[j],1,nchars-4),"_max.txt",sep="")
                        file.copy(from = paste(appSetup$resDir,fileIN,sep="/"),
                                  to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                  overwrite = T)
                        # mean (no extension)
                        fileOUT = paste(substr(allFiles[j],1,nchars-4),"_mean.txt",sep="")
                        file.copy(from = paste(appSetup$resDir,allFiles[j],sep="/"),
                                  to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                  overwrite = T)
                      }
                      basinFiles=c(basinFiles,allFiles[j])
                    }
                  }
                }else{
                  if(app.sys=="tep"){
                    #                    rciop.copy(paste(appSetup$resDir,allFiles[j],sep="/"), outDir, uncompress=TRUE)
                    if(appInput$assimOn=="off"){
                      # Normal run
                      file.copy(from = paste(appSetup$resDir,allFiles[j],sep="/"),
                                to = paste(outDir,paste(prefix.bas,allFiles[j],sep="_"),sep="/"),
                                overwrite = T)
                    }else{
                      # Assimilation run
                      nchars = nchar(allFiles[j])
                      # min (002)
                      fileIN = paste(substr(allFiles[j],1,nchars-4),"_002.txt",sep="")
                      fileOUT = paste(substr(allFiles[j],1,nchars-4),"_min.txt",sep="")
                      file.copy(from = paste(appSetup$resDir,fileIN,sep="/"),
                                to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                overwrite = T)
                      # max (003)
                      fileIN = paste(substr(allFiles[j],1,nchars-4),"_003.txt",sep="")
                      fileOUT = paste(substr(allFiles[j],1,nchars-4),"_max.txt",sep="")
                      file.copy(from = paste(appSetup$resDir,fileIN,sep="/"),
                                to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                overwrite = T)
                      # mean (no extension)
                      fileOUT = paste(substr(allFiles[j],1,nchars-4),"_mean.txt",sep="")
                      file.copy(from = paste(appSetup$resDir,allFiles[j],sep="/"),
                                to = paste(outDir,paste(prefix.bas,fileOUT,sep="_"),sep="/"),
                                overwrite = T)
                    }
                    basinFiles=c(basinFiles,allFiles[j])
                  }
                }
              }
            }
          }
        }
      }
      # transform basinoutput files to csv format
      if(length(basinFiles)>0){
        if(hype2csvExists){
          for(i in 1:length(basinFiles)){
            resCsv = basinfiles2csv(hypeFile=paste(appSetup$resDir,basinFiles[i],sep="/"),
                                    csvFile=paste(outDir,paste(prefix.csv,"_",substr(basinFiles[i],1,nchar(basinFiles[i])-3),"csv",sep=""),sep="/"),
                                    hype2csvFile=hype2csvFile,assimOn=appInput$assimOn)
          }
        }
      }
      
      # plot content of basinoutput files
      # ---------------------------------
      # plotting is made through system call to a separate plotting script, see support issue #5720
      #
      # short story: a) basic plotting is disabled since X11 is not available
      #              b) alternative plotting device is enabled by Cairo package.
      #              c) Cairo package require jpeglib version 9 which is in 
      #                 unresolvable conflict with some dependancies from package rgdal
      #
      if(length(basinFiles)>0){
        for(i in 1:length(basinFiles)){
          if(appInput$assimOn=="off"){
            syscmd = paste(app.rscript4plotting,"--vanilla --slave --quite",app.plotscriptBasinOutput,
                           appSetup$resDir,outDir,basinFiles[i],appSetup$modelName,
                           hype2csvFile,prefix.img,appInput$assimOn,"0","99999999",sep=" ")
          }else{
            nSimObsTemp=length(appInput$assimVar)
            if(nSimObsTemp>0){
              nSimObs=0
              for(j in 1:nSimObsTemp){
                if(appInput$assimVar[j]!="OPEN,LOOP"){
                  nSimObs = nSimObs + 1
                  if(nSimObs==1){
                    obsSimVar = paste(substr(appInput$assimVar[j],1,4),substr(appInput$assimVar[j],6,9),sep="")
                  }else{
                    obsSimVar = c(obsSimVar,paste(substr(appInput$assimVar[j],1,4),substr(appInput$assimVar[j],6,9),sep=""))
                  }
                }
              }
              if(nSimObs==0){
                obsSimVar="99999999"
              }
              nSimObs=as.character(nSimObs)
            }else{
              nSimObs="0"
              obsSimVar="99999999"
            }
            syscmd = paste(app.rscript4plotting,"--vanilla --slave --quite",app.plotscriptBasinOutput,
                           appSetup$resDir,outDir,basinFiles[i],appSetup$modelName,
                           hype2csvFile,prefix.img,appInput$assimOn,nSimObs,obsSimVar,sep=" ")
          }
          if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying plot script:  ",syscmd,sep=""), "/util/R/hypeapps-utils.R")}
          plotres = system(command = syscmd,intern = T)
          if(app.sys=="tep"){rciop.log ("DEBUG", paste(" plot result:  ",plotres,sep=""), "/util/R/hypeapps-utils.R")}
        }
      }
      # plot map output files
      if(length(mapFiles)>0){
        # cdate and edate text strings for plot header
        cdateTXT = as.character(modelForcing$cdate,format="%Y%m%d")
        edateTXT = as.character(modelForcing$edate,format="%Y%m%d")
        for(i in 1:length(mapFiles)){
          if(app.sys=="tep"){
            if(appInput$assimOn=="off"){
              # call mapoutput plot function
              syscmd = paste(app.rscript4plotting,"--vanilla --slave --quite",app.plotscriptMapOutput,
                             appSetup$resDir,outDir,mapFiles[i],appSetup$modelName,
                             appSetup$shapefileRdata,prefix.img,cdateTXT,edateTXT,sep=" ")
              if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying map output plot script:  ",syscmd,sep=""), "/util/R/hypeapps-utils.R")}
              plotres = system(command = syscmd,intern = T)
              if(app.sys=="tep"){rciop.log ("DEBUG", paste(" map output plot result:  ",plotres,sep=""), "/util/R/hypeapps-utils.R")}
            }else{
              # skip file extension that indicates "max" or "min"
              nchars = nchar(mapFiles[i])
              last4  = substr(mapFiles[i],nchars-7,nchars-4)
              if(last4!="_002" & last4!="_003"){
                # call mapoutput plot function
                syscmd = paste(app.rscript4plotting,"--vanilla --slave --quite",app.plotscriptMapOutput,
                               appSetup$resDir,outDir,mapFiles[i],appSetup$modelName,
                               appSetup$shapefileRdata,prefix.img,cdateTXT,edateTXT,sep=" ")
                if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying map output plot script:  ",syscmd,sep=""), "/util/R/hypeapps-utils.R")}
                plotres = system(command = syscmd,intern = T)
                if(app.sys=="tep"){rciop.log ("DEBUG", paste(" map output plot result:  ",plotres,sep=""), "/util/R/hypeapps-utils.R")}
              }
            }
          }
        }
      }
    }
    ## list all files in output folder:
    outFiles=dir(outDir,all.files = F,full.names = T,recursive = F)
    
  }else if(appSetup$appName=="forecast"){
    
    # add hindcast and forecast subfolders to outDir
    outDir = paste(appSetup$tmpDir,'output/hindcast',sep="/")
    outDir = c(outDir,paste(appSetup$tmpDir,'output/forecast',sep="/"))
    dir.create(outDir[1],recursive = T,showWarnings = F)
    dir.create(outDir[2],recursive = T,showWarnings = F)
    
    # loop over hindcast and forecast
    for(k in 1:2){
      
      if(k==1){
        prodTag="hindcast"
      }else{
        prodTag="forecast"
      }
        
      # copy log-files from rundir to outdirs (only when k==1)
      if(k==1){
        hyssLogFile = dir(path = appSetup$runDir, pattern =".log")
        if(app.sys=="tep"){
          if(length(hyssLogFile)>=0){
            for(j in 1:length(hyssLogFile)){
              file.copy(from = paste(appSetup$runDir,hyssLogFile[j],sep="/"), 
                        to = paste(outDir[k],paste(prefix.log,hyssLogFile[j],sep="_"),sep="/"))
            }
          }
        }
      }
      
      # list files in result folder
      timeFiles   = dir(path = appSetup$resDir[k] , pattern ="time")
      mapFiles    = dir(path = appSetup$resDir[k] , pattern ="map")
      subassFiles = dir(path = appSetup$resDir[k] , pattern ="subass")
      simassFile  = dir(path = appSetup$resDir[k] , pattern ="simass")
      allFiles    = dir(path = appSetup$resDir[k] , pattern =".txt")
      
      # copy time files
      if(length(timeFiles)>0){
        for(i in 1:length(timeFiles)){
          if(app.sys=="tep"){
            file.copy(from = paste(appSetup$resDir[k],timeFiles[i],sep="/"), 
                      to = paste(outDir[k],paste(prefix.tim,prodTag,timeFiles[i],sep="_"),sep="/"))
          }
        }
      }
      # copy map files
      if(length(mapFiles)>0){
        for(i in 1:length(mapFiles)){
          if(app.sys=="tep"){
            file.copy(from = paste(appSetup$resDir[k],mapFiles[i],sep="/"), 
                      to = paste(outDir[k],paste(prefix.map,prodTag,mapFiles[i],sep="_"),sep="/"))
          }
        }
      }
      # copy subass files
      if(length(subassFiles)>0){
        for(i in 1:length(subassFiles)){
          if(app.sys=="tep"){
            file.copy(from = paste(appSetup$resDir[k],subassFiles[i],sep="/"), 
                      to = paste(outDir[k],paste(prefix.oth,prodTag,subassFiles[i],sep="_"),sep="/"))
          }
        }
      }
      # copy simass files
      if(length(simassFile)>0){
        if(app.sys=="tep"){
          file.copy(from = paste(appSetup$resDir[k],simassFile[i],sep="/"), 
                    to = paste(outDir[k],paste(prefix.oth,prodTag,simassFile[i],sep="_"),sep="/"))
        }
      }
      
      # basin outputfiles
      outbasins = strsplit(appInput$outbasins,split = ",")[[1]]
      zeroString="0000000000000000000000000000000"
      basinFiles=NULL
      if(length(outbasins)>0){
        for(i in 1:length(outbasins)){
          outFile=paste(outbasins[i],".txt",sep="")
          ni=nchar(outFile)
          for(j in 1:length(allFiles)){
            nj=nchar(allFiles[j])
            if(nj>=ni){
              if(substr(allFiles[j],nj-ni+1,nj)==outFile){
                if(nj>ni){
                  if(substr(allFiles[j],1,nj-ni)==substr(zeroString,1,nj-ni)){
                    if(app.sys=="tep"){
                      file.copy(from = paste(appSetup$resDir[k],allFiles[j],sep="/"), 
                                to = paste(outDir[k],paste(prefix.bas,prodTag,allFiles[j],sep="_"),sep="/"))
                      basinFiles=c(basinFiles,allFiles[j])
                    }
                  }
                }else{
                  if(app.sys=="tep"){
                    file.copy(from = paste(appSetup$resDir[k],allFiles[j],sep="/"), 
                              to = paste(outDir[k],paste(prefix.bas,prodTag,allFiles[j],sep="_"),sep="/"))
                    basinFiles=c(basinFiles,allFiles[j])
                  }
                }
              }
            }
          }
        }
      }
      # transform basinoutput files to csv format
      if(length(basinFiles)>0){
        for(i in 1:length(basinFiles)){
          resCsv = basinfiles2csv(hypeFile=paste(appSetup$resDir[k],basinFiles[i],sep="/"),
                                  csvFile=paste(outDir[k],paste(prefix.csv,"_",prodTag,"_",substr(basinFiles[i],1,nchar(basinFiles[i])-3),"csv",sep=""),sep="/"),
                                  hype2csvFile=hype2csvFile,assimOn=appInput$assimOn)
        }
      }
      
      # HINDCAST and FORECAST standard plots (basin outputs and maps)
      
      # FORECAST special COUT plots
      if(k==2){
        # make warning level plots only if return period level file exists
        if(!is.null(appSetup$rpFileCOUT)){

          rpFile=appSetup$rpFileCOUT
          
          # plot forecast hydrographs for selected subbasins
          if(length(basinFiles)>0){
            for(i in 1:length(basinFiles)){
              syscmd = paste(app.rscript4plotting,
                             "--vanilla --slave --quite",
                             app.plotscriptForecastBasin,
                             appSetup$resDir[1],
                             appSetup$resDir[2],
                             outDir[2],
                             basinFiles[i],
                             appSetup$modelName,
                             hype2csvFile,
                             rpFile,
                             paste(prefix.img,"_forecast",sep=""),
                             sep=" ")
              if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying forecast basin plot script:  ",syscmd,sep=""), "/util/R/hypeapps-utils.R")}
              plotres = system(command = syscmd,intern = T)
              if(app.sys=="tep"){rciop.log ("DEBUG", paste(" plot result:  ",plotres,sep=""), "/util/R/hypeapps-utils.R")}
            }
          }
          
          # plot forecast warning level maps
#           args         = commandArgs(trailingOnly=TRUE)
#           fileDir      = args[1]   # folder with timeOutput files (inputs to this script)
#           plotDir      = args[2]   # folder where the output files should be written (warning levels by text and by plot)
#           name.hypeout = args[3]   # timeOutput filename (input) 
#           name.retlev  = args[4]   # return period level file (input, full path) 
#           name.wl.txt  = args[5]   # warning level text output filename
#           name.wl.png  = args[6]   # warning level png  output filename
#           rdataFile    = args[7]   # path to rdata file with subbasin spatial polygons data frame
#           modelNameIN  = args[8]
          
        
          name.hypeout = "timeCOUT.txt"
          if(file.exists(paste(appSetup$resDir[2],name.hypeout,sep="/"))){
            name.retlev  = appSetup$rpFileCOUT
            name.wl.txt  = paste(prefix.wl.txt,"_forecast_mapWarningLevel.txt",sep="")
            name.wl.png  = paste(prefix.wl.png,"_forecast_mapWarningLevel.png",sep="")
            rdataFile    = appSetup$shapefileRdata
          
            syscmd = paste(app.rscript4plotting,"--vanilla --slave --quite",
                           app.plotscriptWarningLevelMap,
                           appSetup$resDir[2],
                           outDir[2],
                           name.hypeout,
                           name.retlev,
                           name.wl.txt,
                           name.wl.png,
                           rdataFile,
                           appSetup$modelName,
                           sep=" ")
            if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying warning level map plot script:  ",syscmd,sep=""), "/util/R/hypeapps-utils.R")}
            plotres = system(command = syscmd,intern = T)
            if(app.sys=="tep"){rciop.log ("DEBUG", paste(" plot result:  ",plotres,sep=""), "/util/R/hypeapps-utils.R")}
          }
        }
      }  
    }
    #outDir = paste(appSetup$tmpDir,'output',sep="/")
    
    ## list all files in output folders:
    outFiles=dir(outDir[1],all.files = F,full.names = T,recursive = F)
    outFiles=c(outFiles,dir(outDir[2],all.files = F,full.names = T,recursive = F))
    
  }
#  ## return outDir
#  return(outDir)

  ## return outFiles
  return(outFiles)
  
}

# functions for application logfile that will be published as part of application results
appLogOpen<-function(appName,tmpDir,appDate,prefix=NULL){
  fileName=paste(appDate,"_","hypeapps-",appName,".log",sep="")
  if(!is.null(prefix)){
    fileName = paste(prefix,"_",fileName,sep="")
  }
  fileName = paste(tmpDir,"/",fileName,sep="")
  fileConn<-file(fileName,open="wt")
  writeLines(paste("hypeapps-",appName," starting, ",as.character(date()),sep=""),fileConn)
  return(list("fileName"=fileName,"fileConn"=fileConn))
}
appLogWrite<-function(logText,fileConn){
  writeLines(paste(logText,", ",as.character(date()),sep=""),fileConn)
  return(0)
}
appLogClose<-function(appName,fileConn){
  writeLines(paste("hypeapps-",appName," ending, ",as.character(date()),sep=""),fileConn)
  close(fileConn)
  return(0)
}

# internal log succesful sourcing of file
if(app.sys=="tep"){rciop.log ("DEBUG", paste("all functions sourced"), "/util/R/hypeapps-utils.R")}
