#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-environment.R
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
# hypeapps-environment.R: Paths and other settings for h-TEP hydrological modelling applications
# Author:                 David Gustafsson, SMHI
# Version:                2017-09-22
#

## set system flag, if not set
if(!exists("app.sys")){
  app.sys ="tep"
}

## load rciop library
if(app.sys=="tep"){
  library("rciop")
}

if(app.sys=="tep"){rciop.log ("DEBUG", paste("rciop package loaded"), "/util/R/hypeapps-environment.R")}

## Folder settings
if(app.sys=="tep"){
  app.app_path = Sys.getenv("_CIOP_APPLICATION_PATH")
  app.tmp_path = TMPDIR
}else if(app.sys=="win"){
  app.app_path = paste(getwd(),'application',sep="/")
  app.tmp_path = paste(getwd(),'tmp',sep="/")
}else{
  app.app_path = ""
  app.tmp_path = ""
}
if(app.sys=="tep"){rciop.log ("DEBUG", paste("application paths set"), "/util/R/hypeapps-environment.R")}

## settings for plotting with Cairo (TEP only)
if(app.sys=="tep"){
  app.rscript4plotting = "/opt/anaconda/envs/cairo-env/bin/Rscript"
}else if(app.sys=="win"){
  app.rscript4plotting = "Rscript.exe"
}else{
  app.rscript4plotting = "Rscript"
}
app.plotscriptBasinOutput    = paste(app.app_path,"util/R/hypeapps-plot-basinoutput.R",sep="/")
app.plotscriptMapOutput      = paste(app.app_path,"util/R/hypeapps-plot-mapoutput.R",sep="/")
app.plotscriptForecastBasin  = paste(app.app_path,"util/R/hypeapps-plot-forecast-basin.R",sep="/")
app.plotscriptWarningLevelMap= paste(app.app_path,"util/R/hypeapps-plot-warninglevel-map.R",sep="/")

## Model settings and data access information
if(app.sys=="tep"){source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/hypeapps-model-settings.R",sep="/"))}
if(app.sys=="tep"){rciop.log ("DEBUG", paste("model and data access settings sourced"), "/util/R/hypeapps-environment.R")}

## Flag that this file has been sourced
app.envset=T
