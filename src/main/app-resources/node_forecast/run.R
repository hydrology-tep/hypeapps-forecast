#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-forecast/src/main/app-resources/node_forecast/run.R

# Copyright 2017 SMHI
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

# Application 1: "Niger-HYPE 10 day forecast" (hypeapps-forecast)
# Author:         David Gustafsson, SMHI
# Version:        2017-09-21

# Workflow overview:
# ------------------
# 1 Initialization          (load environmental variables, libraries, utility functions, etc)
# 2 Application inputs      (read all user inputs)
# 3 Application setup       (setup working folders, copy model parameters)
# 4 Hindcast input          (hindcast forcing, initial state, info)
# 5 Run hindcast            (run the hindcast model)
# 6 Forecast input          (forecast forcing, initial state, info)
# 7 Run forecast            (run the forecast model)
# 8 Output                  (pepare and publish output data)
# 9 End of workflow

#################################################################################
## 1 - Initialization
## ------------------------------------------------------------------------------
## set application name
app.name = "forecast"
## ------------------------------------------------------------------------------
## flag which environment is used, if not set
if(!exists("app.sys")){
  app.sys ="tep"
}
## ------------------------------------------------------------------------------
## load rciop package and set working directory to TMPDIR when running on TEP 
if(app.sys=="tep"){
  library("rciop")
  
  rciop.log ("DEBUG", " *** hypeapps-forecast *** TEP hydrological modelling applications ***", "/node_forecast/run.R")
  rciop.log ("DEBUG", " rciop library loaded", "/node_forecast/run.R")
  
  setwd(TMPDIR)
  rciop.log("DEBUG", paste(" R session working directory set to ",TMPDIR,sep=""), "/node_forecast/run.R")
}
## ------------------------------------------------------------------------------
## load hypeapps environment and additional R utility functions
if(app.sys=="tep"){
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/hypeapps-environment.R",sep="/"))
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/hypeapps-utils.R", sep="/"))
  rciop.log ("DEBUG", paste(" libraries loaded and utilities sourced"), "/node_forecast/run.R")
}else if(app.sys=="win"){
  source("application/util/R/hypeapps-environment.R")  
  source("application/util/R/hypeapps-utils.R")
}
#################################################################################
## 2 - Application user inputs
## ------------------------------------------------------------------------------
## application input parameters
app.input <- getHypeAppInput(appName = app.name)

if(app.sys=="tep"){rciop.log ("DEBUG", paste(" hypeapps inputs and parameters read"), "/node_forecast/run.R")}

#################################################################################
## 3 - Application setup
## ------------------------------------------------------------------------------
## Prepare basic model setup (static input files and hype model executable copied to working folder)
app.setup <- getHypeAppSetup(modelName = model.name,
                             modelBin  = model.bin,
                             tmpDir    = app.tmp_path,
                             appDir    = app.app_path,
                             appName   = app.name,
                             appInput  = app.input,
                             modelFilesURL = model.files.url,
                             forcingArchiveURL = forcing.archive.url,
                             stateFilesURL = state.files.url,
                             stateFilesIN = state.files)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("HypeApp setup read"), "/node_forecast/run.R")}

#################################################################################
## 4 - Hindcast input data
## ------------------------------------------------------------------------------
## forcing data
hindcast.forcing <- getModelForcing(appSetup   = app.setup,
                                    appInput  = app.input,
                                    dataSource = forcing.data.source,
                                    hindcast   = T)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("hindcast forcing set"), "/node_forecast/run.R")}

## ------------------------------------------------------------------------------
## get Xobs input file(s) from open catalogue
xobs.data <- getXobsData(appInput = app.input,
                         appSetup = app.setup)
if(app.sys=="tep"){rciop.log ("DEBUG", paste("xobs data downloaded from catalogue"), "/node_historical/run.R")}

## ------------------------------------------------------------------------------
## read downloaded Xobs input file(s) - merge into one Xobs.txt in the model run folder
xobs.input <- readXobsData(appSetup = app.setup,
                           xobsData = xobs.data)
if(app.sys=="tep"){rciop.log ("DEBUG", paste("xobs data merged to model rundir"), "/node_historical/run.R")}

## ------------------------------------------------------------------------------
## modify some model files based on input parameters
hindcast.input <- updateModelInput(appSetup = app.setup, appInput = app.input, 
                                   hindcast = T, modelForcing = hindcast.forcing, xobsInput = xobs.input)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("hindcast inputs modified"), "/node_forecast/run.R")}

#################################################################################
## 5 - Run hindcast
## ------------------------------------------------------------------------------
##  run hindcast
if(hindcast.input==0){
  if(app.sys=="tep"){rciop.log ("DEBUG", " ...starting hindcast model run", "/node_forecast/run.R")}
  hindcast.run = system(command = app.setup$runCommand,intern = T)
}
if(app.sys=="tep"){rciop.log ("DEBUG", " ...hindcast model run ready", "/node_forecast/run.R")}

#################################################################################
## 6 - Forecast input data
## ------------------------------------------------------------------------------
## forcing data
forecast.forcing <- getModelForcing(appSetup   = app.setup,
                                    appInput  = app.input,
                                    dataSource = forcing.data.source,
                                    hindcast   = F)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("...forecast forcing set"), "/node_forecast/run.R")}

## ------------------------------------------------------------------------------
## modify some model files based on input parameters
forecast.input <- updateModelInput(appSetup = app.setup, appInput = app.input, 
                                   hindcast = F, modelForcing = forecast.forcing, xobsInput = NULL)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("...forecast inputs modified"), "/node_forecast/run.R")}

#################################################################################
## 7 - Run forecast
## ------------------------------------------------------------------------------
##  run model
if(forecast.input==0){
  if(app.sys=="tep"){rciop.log ("DEBUG", " ...starting forecast model run", "/node_forecast/run.R")}
  forecast.run = system(command = app.setup$runCommand,intern = T)
}
if(app.sys=="tep"){rciop.log ("DEBUG", " ...forecast model run ready", "/node_forecast/run.R")}

#################################################################################
## 8 - Output
## ------------------------------------------------------------------------------
## post-process output data
#if(attr(forecast.run,"status")==1){
#  rciop.publish(paste(app.setup$runDir,"/*",sep=""),recursive=TRUE,metalink=TRUE) 
#}else{
app.outdir <- prepareHypeAppsOutput(appSetup  = app.setup, appInput = app.input, 
                                    modelInput = forecast.input, modelForcing = forecast.forcing,
                                    runRes = attr(forecast.run,"status"))

## ------------------------------------------------------------------------------
## publish postprocessed results
if(app.sys=="tep"){
  for(k in 1:length(app.outdir)){
    rciop.publish(path=paste(app.outdir[k],"/*",sep=""), recursive=FALSE, metalink=TRUE)
  }
}
#}
#################################################################################
## 9 - End of workflow
## ------------------------------------------------------------------------------
## exit with appropriate status code
q(save="no", status = 0)
