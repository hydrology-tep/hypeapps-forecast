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
# hypeapps-model-settings.R: Model name, URL to store.terradue.com etc.
# Author:                    David Gustafsson, SMHI
# Version:                   2019-05-16

## [DO NOT EDIT] set all required variables to NULL
model.name      = NULL
model.bin       = NULL
model.files.url = NULL
forcing.files   = NULL
forcing.data.source = NULL
forcing.archive.url = NULL
forcing.archive.start = NULL
forcing.archive.end = NULL
state.files.url = NULL
state.files = NULL
hype2csv.url = NULL
hype2csv.file = NULL
shapefile.url   = NULL
shapefile.layer = NULL
shapefile.ext   = c(".shp",".prj",".dbf",".shx")

## [EDIT HERE] Model name and executable binary
model.name = "niger-hype"                                            # model name (used for model input file folder, shapefile, in plot legends etc)
model.bin  = "hype_assimilation-5.x.0.exe"                           # name of model binary executable

## [EDIT HERE] Data store settings (model files, forcing archive)
model.files.url = "https://store.terradue.com/dgustafsson/hype-files/model/niger-hype/v2.23" # model files root index

## [EDIT HERE] Forcing data settings
forcing.files = c("Pobs.txt","Tobs.txt","TMINobs.txt","TMAXobs.txt") # list of necessary forcing data files

forcing.data.source =  "hydro-smhi"  # could be "local" or "hydro-smhi"

forcing.archive.url  = "https://store.terradue.com/dgustafsson/hype-files/model/niger-hype/v2.23/forcingarchive"
forcing.archive.start   = as.POSIXct("1979-01-01",tz="GMT")                    # first date in the forcing arhive data
forcing.archive.end     = as.POSIXct("2018-12-31",tz="GMT")                    # last date in the forcing arhive data

# [EDIT HERE] Statefiles settings
state.files.url = "https://store.terradue.com/dgustafsson/hype-files/model/niger-hype/v2.23/statefiles"
state.files="state_save19790101.txt"
for(i in 1980:2019){
  state.files=c(state.files,paste("state_save",as.character(i),"0101.txt",sep=""))
}

# [EDIT HERE] hype2csv URL
hype2csv.url = "https://store.terradue.com/dgustafsson/hype-files/model/niger-hype/v2.23/hype2csv/niger-hype2csv.txt"
hype2csv.file = "niger-hype2csv.txt"

# [EDIT HERE] sub-basin shapefile URL (shapefile.url should point to shapefile [model.name].shp, and in the same folder should be .dbf, .prj and .shx)
shapefile.url   = "https://store.terradue.com/dgustafsson/hype-files/model/niger-hype/v2.23/shapefiles"
shapefile.layer = "niger-hype"
shapefile.ext   = c(".shp",".prj",".dbf",".shx") 

# log message
if(app.sys=="tep"){rciop.log ("DEBUG", paste("model and data access settings set"), "/util/R/hypeapps-model-settings.R")}
