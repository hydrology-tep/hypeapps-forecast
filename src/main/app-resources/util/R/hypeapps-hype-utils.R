#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hype-utils.R
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
# hype-utils.R: HYPE model input/output utility functions (HYPEtools light)
# Author:       David Gustafsson, SMHI
# Version:      2017-03-31

# Functions:
#
#   writeInfo
#   readInfoLine
#   readInfoDate
#   readInfo
#   WriteXobs
#   ReadBasinOutput
#   ReadTimeOutput
#   writeHype2CSV

# dependencies
library(data.table)

#-------------------------------------------
#functions to read/write and modify Info.txt
#-------------------------------------------
#write Info.txt
writeInfo<-function(x,filenm="Info.txt"){
  #x is a string vector with one record for each line in Info.txt
  write(x,filenm)
  return(0)
}

#read data from a line in Info.txt
readInfoLine<-function(thisline,skipChar){
  #thisline is a string for one line in info.txt
  #skipChar is the number of characters to skip before reading the data
  thisline=substr(thisline,skipChar+1,nchar(thisline))
  #look for TAB or space
  if(regexpr("\t",thisline)[[1]]>=1){
    #TAB delimited
    while(regexpr("\t",thisline)[[1]]==1){
      thisline=substr(thisline,2,nchar(thisline))
    }
    infoval=thisline
  }else{
    #SPACE delimited
    while(regexpr(" ",thisline)[[1]]==1){
      thisline=substr(thisline,2,nchar(thisline))
    }
    infoval=thisline
  }
  return(infoval)
}

#read a date from a line in info
readInfoDate<-function(thisline){
  #the date is assumed to be given in the format YYYY-mm-dd
  #function could be further generalized...
  iSep = regexpr("-",thisline)[[1]]
  infodate=as.Date(substr(thisline,iSep-4,iSep+5))
}

#main function to read info.txt into a useful list for easy modification
readInfo<-function(filenm = "Info.txt"){ 
  #
  # updated to HYPE_DA_v5.x.0 (2017-07-14)
  #
  #read info file into a list called info
  info.lines = readLines(filenm)
  info=list("info.lines"=info.lines)
  
  #some default flags
  info$isResDir=F
  info$isModDir=F
  
  #loop over lines and identify some relevant data
  noutstate = 0
  for(i in 1:length(info$info.lines)){
    thisline = info$info.lines[i]
    #bdate        2012-10-01
    if(substr(thisline,1,5)=="bdate"){
      #iSep = regexpr("-",thisline)[[1]]
      #info$bdate=substr(thisline,iSep-4,iSep+5)
      info$bdate=readInfoDate(thisline)
      info$bdate.lineNr=i
    }
    #cdate      	2013-01-10
    if(substr(thisline,1,5)=="cdate"){
      info$cdate=readInfoDate(thisline)
      info$cdate.lineNr=i
    }
    #edate      	2013-07-31
    if(substr(thisline,1,5)=="edate"){
      info$edate=readInfoDate(thisline)
      info$edate.lineNr=i
    }
    #resultdir	
    if(substr(thisline,1,9)=="resultdir"){
      info$resultdir=readInfoLine(thisline,9)
      info$resultdir.lineNr=i
      info$isResDir=T
    }
    #modeldir	
    if(substr(thisline,1,8)=="modeldir"){
      info$modeldir=readInfoLine(thisline,8)
      info$modeldir.lineNr=i
      info$isModDir=T
    }
    #instate	n
    if(substr(thisline,1,7)=="instate"){
      info$instate=readInfoLine(thisline,7)
      info$instate.lineNr=i
    }
    #outstatedate	2010-09-30
    if(substr(thisline,1,12)=="outstatedate"){
      if(noutstate==0){
        noutstate=1
        info$outstatedate=readInfoDate(thisline)
        info$outstatedate.lineNr=i
      }else{
        noutstate=noutstate+1
        info$outstatedate=c(info$outstatedate,readInfoDate(thisline))
        info$outstatedate.lineNr=c(info$outstatedate.lineNr,i)
      }
    }
    #assimilation          y/n
    if(substr(thisline,1,12)=="assimilation"){
      info$assimilation=readInfoLine(thisline,12)
      info$assimilation.lineNr=i
    }
    # basinoutput variable
    if(substr(thisline,1,20)=="basinoutput variable"){
      info$basinoutput_variable=readInfoLine(thisline,20)
      info$basinoutput_variable.lineNr=i
    }
    # basinoutput subbasin
    if(substr(thisline,1,20)=="basinoutput subbasin"){
      info$basinoutput_subbasin=readInfoLine(thisline,20)
      info$basinoutput_subbasin.lineNr=i
    }
    # timeoutput variables
    if(substr(thisline,1,19)=="timeoutput variable"){
      info$timeoutput_variable=readInfoLine(thisline,19)
      info$timeoutput_variable.lineNr=i
    }
    
    # mapoutput variables
    if(substr(thisline,1,18)=="mapoutput variable"){
      info$mapoutput_variable=readInfoLine(thisline,18)
      info$mapoutput_variable.lineNr=i
    }
    
  }
  info$noutstate = noutstate      
  return(info)
}


## ------------------------------------------------------------
## HYPE import/export functions from HYPEtools
## ------------------------------------------------------------
##
## WriteXobs --------------------------------------------------
##
WriteXobs <- function(x, filename = "Xobs.txt", append = F, comment = NA, variable = NA, subid = NA, 
                      lastDate = NA, timestep = "d") {
  if (!append) {
    # no appending, just write to a new file
    
    # Export of comment (if it exists) and header information, uses a cascade of conditional checks
    # which will stop the export with informative error messages if inconsistencies are found
    
    # remove existing export file. this makes it possible to have a consistent 'append' connection open below, instead of
    # overwriting first and then reopening in append mode...
    invisible(suppressWarnings(file.remove(filename)))
    
    # create and open a connection to which the following writeLines() will be written
    fcon <- file(description = filename, open ="at")
    
    ## export comment line
    if (is.na(comment)){
      # comment argument is empty
      if(!is.null(attr(x, which = "comment"))) {
        # comment attribute exists, export
        writeLines(paste(attr(x, which = "comment"), collapse = "\t"), con = fcon)
      } else {
        # comment attribute does not exist, export the following string
        writeLines("Exported from R", con = fcon)
      }
    } else {
      # export comment argument
      writeLines(comment, con = fcon)
    }
    
    ## export variable line
    if (length(variable) == 1 & is.na(variable[1])) {
      # variable argument is empty
      if(!is.null(attr(x, which = "variable"))) {
        # attribute variable exists
        if (length(attr(x, which = "variable")) == ncol(x) - 1) {
          # attribute and export dataframe match in length, export attribute with padded 'name' string and newline
          tmp <- paste(c("x", attr(x, which = "variable")), collapse = "\t")
          writeLines(tmp, con = fcon)
        } else {
          # mismatch in length, stop with error
          close(fcon)
          stop("Length of attribute 'variable' does not match number of variables in export object.\n 
               Check consistency, e.g. with attr(x, 'variable') and names(x)[-1].")
        }
        } else {
          # attribute variable does not exist, stop with error
          close(fcon)
          stop("'variable' argument not given and 'variable' attribute not existing in export object.")
      }
    } else {
      # export the variable argument with padded 'name' string and newline, if length matches no. of observation data cols in x
      if (length(variable) == ncol(x) - 1) {
        tmp <- paste(c("x", variable), collapse ="\t")
        writeLines(tmp, con = fcon)
      } else {
        # mismatch in length, stop with error
        close(fcon)
        stop("Length of argument 'variable' does not match number of variables in export object.")
      }
    }
    
    ## export subid line
    if (length(subid == 1) & is.na(subid[1])) {
      # subid argument is empty
      if(!is.null(attr(x, which = "subid"))) {
        # attribute subid exists
        if (length(attr(x, which = "subid")) == ncol(x) - 1) {
          # attribute and export dataframe match in length, export attribute with padded 0 and newline
          tmp <- paste(as.character(c(0, attr(x, which = "subid"))), collapse = "\t")
          writeLines(tmp, con = fcon)
        } else {
          # mismatch in length, stop with error
          close(fcon)
          stop("Length of attribute 'subid' does not match number of variables in export object.\n 
               Check consistency, e.g. with attr(x, 'subid') and names(x)[-1].")
        }
        } else {
          # attribute subid does not exist, stop with error
          close(fcon)
          stop("'subid' argument not given and 'subid' attribute not existing in export object.")
      }
    } else {
      # export the subid argument with padded 0 and newline, if length matches no. of observation data cols in x
      if (length(subid) == ncol(x) - 1) {
        tmp <- paste(as.character(c(0, subid)), collapse = "\t")
        writeLines(tmp, con = fcon)
      } else {
        # mismatch in length, stop with error
        close(fcon)
        stop("Length of argument 'subid' does not match number of variables in export object.")
      }
    }
    
    close(fcon)
    
    
    
  } else {
    # export will be appended to existing file
    
    # check if file to export to exists, stop otherwise
    stopifnot(file.exists(filename))
    
    # read variable names from existing file into two vectors, to be compared against the export data for consistency
    tmp <- readLines(filename,n=3)
    existingVar <- strsplit(tmp[2], split = "\t")[[1]][-1]
    existingSbd <- as.integer(strsplit(tmp[3], split = "\t")[[1]][-1])
    
    # first consistency check: number of columns identical?
    if (length(existingVar) != ncol(x) - 1) {
      stop("Inconsistent number of data columns between export object and existing file.")
    }
    
    ## second consistency check: variables and SUBIDs identical and in the same order?
    # select export data to compare with, either from function argument or from attribute of x
    if (!is.na(variable[1])) {
      exportVar <- variable
    } else {
      if (!is.null(attr(x, "variable"))) {
        exportVar <- attr(x, "variable")
      } else {
        stop("'variable' argument not given and 'variable' attribute not existing in export object.")
      }
    }
    
    if (!is.na(subid[1])) {
      exportSbd <- subid
    } else {
      if (!is.null(attr(x, "subid"))) {
        exportSbd <- attr(x, "subid")
      } else {
        stop("'subid' argument not given and 'subid' attribute not existing in export object.")
      }
    }
    
    # check consistency, will fail if at least one column is different for either subid or variable
    if (any(!(exportVar == existingVar), !(exportSbd == existingSbd))) {
      stop("Inconsistent variable names or SUBIDs.")
    }
    
    ## third consistency check: is last date in existing Xobs earlier than the first export row?
    ## split into hour and day cases
    if (timestep == "h" | timestep == "hourly") {
      tdiff <- difftime(x[1,1], lastDate, units = "hours")
      if (tdiff < 1) {
        stop("Time series in existing and new Xobs overlap or difference is smaller than one hour.")
      }
      ## export '-9999' lines if gap is larger than one hour
      # create date vector
      dpad <- seq(from = lastDate, to = x[1,1], length = tdiff + 1)
      # cut off end and start dates
      dpad <- dpad[-c(1, length(dpad))]
      # create data frame from a matrix (it is more straightforward to span up an empty matrix and then convert to df)
      pad <- cbind(format(dpad, format = "%Y-%m-%d %H:%M"), as.data.frame(matrix(data = -9999, nrow = length(dpad), 
                                                                                 ncol = ncol(x) - 1)))
      write.table(pad, file = filename, col.names = F, sep = "\t", append = T, na = "-9999", row.names = F, quote = F)
    }
    
    if (timestep == "d" | timestep == "daily") {
      tdiff <- difftime(x[1,1], lastDate, units = "days")
      if (tdiff < 1) {
        stop("Time series in existing and new Xobs overlap or difference is smaller than one day.")
      }
      ## export '-9999' lines if gap is larger than one day
      # create date vector
      dpad <- seq(from = lastDate, to = x[1,1], length = tdiff + 1)
      # cut off end and start dates
      dpad <- dpad[-c(1, length(dpad))]
      # create data frame from a matrix (it is more straightforward to span up an empty matrix and then convert to df)
      pad <- cbind(format(dpad, format = "%Y-%m-%d"), as.data.frame(matrix(data = -9999, nrow = length(dpad), 
                                                                           ncol = ncol(x) - 1)))
      write.table(pad, file = filename, col.names = F, sep = "\t", append = T, na = "-9999", row.names = F, quote = F)
    }
  }
  
  # Export of the dataframe, format date-times to HYPE requirements first
  if (timestep == "d" | timestep == "daily") {
    x[,1] <- format(x[,1], format = "%Y-%m-%d")
  }
  if (timestep == "h" | timestep == "hourly") {
    x[,1] <- format(x[,1], format = "%Y-%m-%d %H:%M")
  }
  write.table(x, file = filename, col.names = F, sep = "\t", append = T, na = "-9999", row.names = F, quote = F)
  
  
}

##
## ReadBasinOutput --------------------------------------------------
##
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

## 
## ReadTimeOutput
##
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

###
### MERGE XOBS DATA FRAMES
###
MergeXobs <- function(x, y, comment = "") {
  
#   # check time step with in both inputs and if they are identical
#   # requires equidistant time steps
#   if (!any(class(x) == "HypeXobs")) {
#     warning("'x' not of class HypeXobs.")
#     x.tstep <- difftime(x[2, 1], x[1, 1])
#   } else {
#     x.tstep <- attr(x, "timestep")
#   }
#   if (!any(class(y) == "HypeXobs")) {
#     warning("'y' not of class HypeXobs.")
#     y.tstep <- difftime(x[2, 1], x[1, 1])
#   } else {
#     y.tstep <- attr(y, "timestep")
#   }
#   
#   if (x.tstep != y.tstep) {
#     stop("Time step lengths in 'x' and 'y' differ.")
#   }
#   
#   # check if there are any duplicated dates in x or y
#   if (anyDuplicated(x[, 1]) || anyDuplicated(y[, 1])) {
#     stop("Duplicated dates in either 'x' or 'y'.")
#   }
#   
  # Create data frame with common time period column for both input xobs in appropriate time steps 
  date.min <- min(min(x[, 1]), min(y[, 1]))
  date.max <- max(max(x[, 1]), max(y[, 1]))
  res <- data.frame(date = seq(date.min, date.max, by = x.tstep))
  
  
  # merge x and y with new time axis individually
  res1 <- merge(res, x, by = 1, all = TRUE)
#  attr(res1, "comment") <- attr(x, "comment")
  attr(res1, "variable") <- attr(x, "variable")
  attr(res1, "subid") <- attr(x, "subid")
#  attr(res1, "class") <- attr(x, "class")
#  attr(res1, "timestep") <- attr(x, "timestep")
  
  res2 <- merge(res, y, by = 1, all = TRUE)
#  attr(res2, "comment") <- attr(y, "comment")
  attr(res2, "variable") <- attr(y, "variable")
  attr(res2, "subid") <- attr(y, "subid")
#  attr(res2, "class") <- attr(y, "class")
#  attr(res2, "timestep") <- attr(y, "timestep")
  
  
  # extract variable names from the merged data and match common column where observations have to be merged
  names1 <- paste0(attr(res1, "variable"), "_", attr(res1, "subid"))
  names2 <- paste0(attr(res2, "variable"), "_", attr(res2, "subid"))
  common.cols <- match(names1, names2)
  
  # conditional: common columns exist, merge them
  if(length(na.omit(common.cols)) > 0) {
    cat("Common columns found, merging.\n")
    cat(paste0("Common column indices in 'x': ", paste(which(!is.na(common.cols)) + 1, collapse = " "), "\n"))
    cat(paste0("Common column indices in 'y': ", paste(as.integer(na.omit(common.cols)) + 1, collapse = " "), "\n"))
    # columns to merge, res1
    te1 <- res1[, c(TRUE, !is.na(common.cols))]
    # columns to merge, res2
    te2 <- res2[, c(1, as.integer(na.omit(common.cols)) + 1)]
    
    # fill observations from xobs without precedence into the one with precedence, if no obs exist there
    # mapply this to all identified columns (te1 and te2 are ALWAYS of the same length, therefore mapply is safe)
    te3 <- as.data.frame(mapply(function(x, y) {ifelse(!is.na(x), x, y)}, te1, te2))
    
    # update columns in data source with precedence
    res1[, c(FALSE, !is.na(common.cols))] <- te3[, -1]
    
    # remove columns from data source without precedence
    res2 <- res2[, -(as.integer(na.omit(common.cols)) + 1)]
    
  }
  
  # combine the results, catch special case where all columns are common and only a date vector is left in res2
  if(is.data.frame(res2)) {
    res <- suppressWarnings(cbind(res, res1[, -1], res2[, -1]))
  } else {
    res <- suppressWarnings(cbind(res, res1[, -1]))
  }
  
  # update comment attribute, conditional on function argument value
  if (comment == "") {
    comment <- paste0("!Created by MergeXobs. Original comments: ", 
                      attr(x,"comment"), " (x); ", attr(y,"comment"), " (y)")
  }
  
  # reconstruct other HypeXobs attributes from res1 and res2
  res <- HypeXobs(x = res, comment = comment, 
                  variable = c(attr(res1, "variable"), attr(res2, "variable")), 
                  subid = c(attr(res1, "subid"), attr(res2, "subid")))
  
  return(res)
}

##
## ReadPTQobs
##
ReadPTQobs <- function (filename, dt.format = "%Y-%m-%d", nrows = -1) {
  
  ## import ptqobs file header, extract attribute
  # import
  xattr <- readLines(filename,n = 1)
  # extract obsids
  sbd <- as.integer(strsplit(xattr, split = "\t")[[1]][-1])
  
  # read the data
  x <- fread(filename,  na.strings = "-9999", sep = "\t", header = T, data.table = F, nrows = nrows)
  #colClasses = c("NA", rep("numeric", length(sbd))))
  
  attr(x, which = "obsid") <- sbd
  
  # date conversion 
  xd <- as.POSIXct(strptime(x[, 1], format = dt.format), tz = "GMT")
  x[, 1] <- xd
  #x[, 1] <- tryCatch(na.fail(xd), error = function(e) {
  #  print("Date/time conversion attempt led to introduction of NAs, date/times returned as strings"); return(x[, 1])
  # })
  
  return(x)
}

##
## WritePTQobs
##
WritePTQobs <- function (x, filename, dt.format = "%Y-%m-%d", digits = 3, nsmall = 1, obsid = NULL) {
  
  ## check if consistent header information is available, obsid arguments take precedence before attribute
  if(!is.null(obsid)) {
    if (length(obsid) == ncol(x) - 1) {
      header <- c("DATE", obsid)
    } else {
      stop("Length of function argument 'obsid' does not match number of obsid columns in export object.")
    }
  } else if (!is.null(attr(x, which = "obsid"))) {
    if (length(attr(x, which = "obsid")) == ncol(x) - 1) {
      header <- c("DATE", attr(x, which = "obsid"))
    } else {
      stop("Length of attribute 'obsid' does not match number of obsid columns in export object.")
    }
  } else {
    stop("No information available from 'obsid' argument or 'obsid' attribute to construct export header.")
  }
  
  # date conversion, conditional on that the date column is a posix class
  if (any(class(x[, 1]) == "POSIXct")) {
    x[, 1] <- format(x[, 1], format = dt.format)
  } else {
    warning("First column in export data frame is not of class 'POSIXct', will be exported unchanged.")
  }
  
  # convert NAs to -9999, needed because format() below does not allow for automatic replacement of NA strings 
  x[is.na(x)] <- -9999
  
  # export
  write.table(format(x, digits = digits, nsmall = nsmall, scientific = F, drop0trailing = T, trim = T), file = filename, 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = header)
  
}


###
### ReadXobs
###
ReadXobs <- function (filename = "Xobs.txt", dt.format="%Y-%m-%d", nrows = -1L) {
  
  ## import xobs file header, extract attributes
  # import (3-row header)
  xattr <- readLines(filename,n=3)
  # 1st row, comment
  # split string elements along tabs, returns list of character vectors
  cmt <- strsplit(xattr[1], split = "\t")
  # remove empty strings (excel export artefacts)
  cmt <- sapply(cmt, function(x) {te <- nchar(x);te <- ifelse(te == 0, F, T);x[te]})
  # 2nd row, HYPE variable IDs
  hype.var <- toupper(strsplit(xattr[2], split = "\t")[[1]][-1])
  # 3rd row, SUBIDs
  sbd <- as.integer(strsplit(xattr[3], split = "\t")[[1]][-1])
  
  
  # read the data, skip header and comment rows, force numeric data (automatic column classes can be integer)
  xobs <- fread(filename,  na.strings = "-9999", skip = 3, sep = "\t", header = F, data.table = F, nrows = nrows, 
                colClasses = c("NA", rep("numeric", length(sbd))))
  
  # update header, composite of variable and subid
  names(xobs) <- c("DATE", paste(hype.var, sbd, sep = "_"))
  
  # warn if duplicate columns found, throw useful msg
  if (length(names(xobs)) != length(unique(names(xobs)))) {
    warning(paste0("Duplicated variable-SUBID combination(s) in file: ", paste(names(xobs)[duplicated(names(xobs))], collapse = " ")))
    duplifree <- FALSE
  } else {
    duplifree <- TRUE
  }
  
  # date conversion 
  xd <- as.POSIXct(strptime(xobs[, 1], format = dt.format), tz = "GMT")
#  xobs[, 1] <- tryCatch(na.fail(xd), error = function(e) {
#    cat("Date/time conversion attempt led to introduction of NAs, date/times returned as strings.\nImported as data frame, not as 'HypeXobs' object.\n"); return(xobs[, 1])})
  
  
  # if date conversion worked and time steps are HYPE-conform (need at least 2 time steps), make returned object class HypeXobs
#  if(!is.character(xobs[, 1]) && duplifree && nrow(xobs) > 1) {
    
    # create HypeXobs object, can fail if multi-day time steps in imported table
#    xobs <- tryCatch(HypeXobs(x = xobs, comment = cmt, variable = hype.var, subid = sbd), 
#                     error = function(e) {cat("Longer-than-daily time steps not allowed in HypeXobs objects.\n"); return(xobs)})
    
    # update with additional attributes if HypeXobs class assignment failed
#    if (!any(class(xobs) == "HypeXobs")) {
      attr(xobs, which = "comment") <- cmt
      attr(xobs, which = "variable") <- hype.var
      attr(xobs, which = "subid") <- sbd
#      warning("Imported as data frame, not as 'HypeXobs' object.")
#    }
    
#  } else {
#    # update with additional attributes if not a HypeXobs object
#    attr(xobs, which = "comment") <- cmt
#    attr(xobs, which = "variable") <- hype.var
#    attr(xobs, which = "subid") <- sbd
#    
#    warning("Imported as data frame, not as 'HypeXobs' object.")
#  }
  
  return(xobs)
}
