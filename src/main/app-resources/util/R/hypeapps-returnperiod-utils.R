#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-returnperiod-utils.R
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
# hypeapps-returnperiod-utils.R: Return Period Analysis (river discharge)
# Author:                        David Gustafsson, SMHI
# Version:                       2017-07-11

# ------------------------------------------------------
# load necessary libraries and source code
# ------------------------------------------------------

# library lmomco for statistical analysis
library(lmomco)

# hypeapp-hype-utils.R to import HYPE model outputs
if(app.sys=="tep"){
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-hype-utils.R",sep="/"))
}else if(app.sys=="win"){
  source(paste("application","util/R/hypeapps-hype-utils.R",sep="/"))
}else{
  stop("app.sys not set (allowed vaues tep or win)")
}


# ------------------------------------------------------------
# FrequencyAnalysis  - Fits a given extreme value distribution
#                     to an extreme value series
#   INPUT:
#   series       A vector representing an extreme value series (e.g., annual maximum flood)
#   distribution A three-character name of the extreme value distribution (see ?dist.list())
#   nep          A vector of non-exceedance probabilities
#
#   OUTPUT: 
#   A list object containing: (1) distribution information and (2) output
#  (quantile estimates at various non-exceedance probabilities)
# ------------------------------------------------------------
FrequencyAnalysis <- function( series, distribution, nep = nonexceeds() ) {
  library(lmomco)
  distribution <- tolower(distribution)
  transformed <- FALSE
  
  # add log Pearson Type 3 to list of distributions supported
  # by lmomco package
  base.dist <- c('lp3', dist.list())
  
  if( any(distribution %in% base.dist) ) {
    
    # log transform series 
    if( distribution == 'lp3' ) {
      series <- log10(series)
      transformed <- TRUE
      distribution <- 'pe3'
    }
    
    # compute L-moments
    samLmom <- lmom.ub(series)
    
    # estimate distribution parameters
    distPar <- lmom2par(samLmom, type = distribution)
    
    # compute quantiles for nonexceedances
    quant <- par2qua(f = nep, para = distPar)
    
    if( distribution == 'pe3' & transformed ) {
      distribution <- 'lp3'
      quant <- 10^quant
    }
    
    # return result as list object
    return(
      list(
        distribution = list(
          name = distribution,
          logTransformed = transformed,
          parameters = distPar),
        output = data.frame(nep = nep, rp = prob2T(nep), estimate = quant) 
      ) )
    
  } else {
    stop(
      sprintf('Distribution \'%s\' not recognized!', distribution))
  }
}


# ------------------------------------------------------
# annualMaximum - function to identify the annual maximum 
#                 value in a time series data frame
# ------------------------------------------------------
annualMaximum<-function(ts.data){
  # years
  years = as.character(ts.data$DATE)
  years = as.integer(substr(years,1,4))
  
  # initiate output data frame
  output=data.frame("year"=unique(years))
  output$max=NA
  
  # loop through years and find max value
  for(i in 1:nrow(output)){
    iyear=which(years==output$year[i] & !is.na(ts.data[,2]))
    if(length(iyear)>0){
      output$max[i]=max(ts.data[iyear,2])
    }
  }
  
  # end, return
  return(output)
}

# -------------------------------------------------------------------------
# returnPeriodMagnitudes - Calculate return levels for selected return periods
#                          (return-period magnitudes)
# INPUT:
# name.in    Text string with path to input file (timeXXXX.txt)
# name.out   Text string with path to output file (retlev.txt)
# wl.rp      Vector with return periods (years), default 2, 5, 30 years
# dist       Name of distrbution ("gev")
# -------------------------------------------------------------------------
returnPeriodMagnitudes<-function(name.in=NULL,name.out="retlev.txt",wl.rp=c(2,5,30),dist="gev"){
  
  if(file.exists(name.in)){
    # read timefile
    hypeout <- ReadTimeOutput(filename = name.in,dt.format = "%Y-%m-%d")

    # log some file properties
#    if(app.sys=="tep"){rciop.log ("DEBUG", paste("   name.in=",name.in,sep=""), "/util/R/hypeapps-returnperiod-utils.R")}
#    if(app.sys=="tep"){rciop.log ("DEBUG", paste("   ncol(hypeoutput)=",as.character(ncol(hypeout)),sep=""), "/util/R/hypeapps-returnperiod-utils.R")}
#    if(app.sys=="tep"){rciop.log ("DEBUG", paste("   nrow(hypeoutput)=",as.character(nrow(hypeout)),sep=""), "/util/R/hypeapps-returnperiod-utils.R")}
#    if(app.sys=="tep"){rciop.log ("DEBUG", paste("   colnames(hypeoutput)=",colnames(hypeout),sep=""), "/util/R/hypeapps-returnperiod-utils.R")}
    
    # change column names
    colnames(hypeout) <- sub("X","",colnames(hypeout))

    # Define object to store return levels
    retlev<-as.data.frame(matrix(NA,nrow=ncol(hypeout)-1,ncol=1+length(wl.rp)))
    colnames(retlev)<-c("SUBID",paste("RP",wl.rp,sep=""))
    retlev[,"SUBID"]<-colnames(hypeout)[-1]
  
    # Calculat annual maximum for each year and subid
    annmax<-aggregate(hypeout[,-1],by=list(year=as.numeric(format(hypeout[,"DATE"],format="%Y"))),FUN="max")
    
   # if(app.sys=="tep"){rciop.log ("DEBUG", paste("   annmax=",as.character(annmax),sep=""), "/util/R/hypeapps-returnperiod-utils.R")}
    
    
    # Fit the distribution and derive return levels for each subbasin
    for(j in 1:nrow(retlev)){ # j<-2
      # fit frequency distribution
      suppressWarnings(  # gives error and warns if the FrequencyAnalysis fails (e.g. in very dry areas)
        fa.sim <- tryCatch(FrequencyAnalysis(series=annmax[,retlev[j,"SUBID"]], distribution=dist), error = function(e) {NULL})
      )
      if(!is.null(fa.sim)) {
        # extract return levels, note the loop is there in case some RP is chosen that is not included in the fa.sim$output$rp, then it interpolates between the nearest (davids method)
        wl.sim <- rep(0,length(wl.rp)) 
        rp.sim <- as.numeric(fa.sim$output$rp)
        for(i in 1:length(wl.rp)){ #i<-1
          ilow=which(rp.sim<=wl.rp[i])
          ihig=which(rp.sim>=wl.rp[i])
          if(ilow[length(ilow)]==ihig[1]){ # if the RP is in the estimate list
            wl.sim[i]=fa.sim$output$estimate[ihig[1]]
          }else{
            wl.sim[i]=fa.sim$output$estimate[ilow[length(ilow)]] + (wl.rp[i]-rp.sim[ilow[length(ilow)]]) * (fa.sim$output$estimate[ihig[1]]-fa.sim$output$estimate[ilow[length(ilow)]]) / (rp.sim[ihig[1]]-rp.sim[ilow[length(ilow)]])
          }
        }
        
        retlev[j,1+1:length(wl.rp)] <- wl.sim  # insert the derived return levels
        rm(wl.sim,rp.sim)  # Clean
      }
      
      rm(fa.sim)
    }
  
    # Write return period magnitude output file
    write.table(retlev,name.out,quote = F,row.names=F)  # note it writes with 11|12 decimals so reading in the object will have slighlty different values...
  
    # return ok signal
    return(0)
  }else{
    return(-1)
  }
}
