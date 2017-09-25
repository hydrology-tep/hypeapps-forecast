#!/opt/anaconda/envs/cairo-env/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-plot-mapoutput.R
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
# hypeapps-plot-basinoutput.R: Script to plot HYPE mapoutput data, used for TEP Hydrology.
# Author:                      David Gustafsson, SMHI
# Version:                     2017-09-21
# -------------------------------------------------------------------
# 1 - Argument inputs
# -------------------------------------------------------------------
args         = commandArgs(trailingOnly=TRUE)
fileDir      = args[1]   # folder with mapoutput files
plotDir      = args[2]   # folder where the plot files should be written
mapFile      = args[3]   # mapoutput filename
modelName    = args[4]   # modelName (used for plot main title)
rdataFile    = args[5]   # path to rdata file with subbasin spatial polygons data frame
prefix.png   = args[6]   # png filename prefix
cdate        = args[7]   # start of simulation
edate        = args[8]   # end of simulation
 
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
library(sp)         # sp for reading Rdata file

pngORjpg = 0        # 0=> png, 1=>jpg

# -------------------------------------------------------------------
# 4- Local functions to read and plot HYPE map output file
# -------------------------------------------------------------------

### colors from HYPEtools
### ---------------------
ColNitr <- colorRampPalette(c("#fff5a8", "#6b0601"))
ColPhos <- colorRampPalette(c("#dcf5e9", "#226633"))
ColPrec <- colorRampPalette(c("#e0e7e8", "#00508c"))
ColTemp <- colorRampPalette(c("#0000ff", "#0080ff", "#80ffff", "#f0f0f0", "#ffff00", "#ff8000", "#ff0000"))
ColQ <- colorRampPalette(c("#ede7ff", "#2300ff"))
ColDiffTemp <- colorRampPalette(c("#2892c7", "#e9f0e8", "#e81515"))
ColDiffGeneric <- colorRampPalette(c("#e81515", "#e9f0e8", "#2892c7"))
ColBlues <- colorRampPalette(c("#0a0a96", "#a3a3db"))
ColReds <- colorRampPalette(c("#f77497", "#670101"))
ColGreens <- colorRampPalette(c("#04eb04", "#004400"))
ColYOB <- colorRampPalette(c("#ffe851", "#da531d", "#5b1e00"))
ColPurples <- colorRampPalette(c("#da62ed", "#300275"))

### PlotMapOutput from HYPEtools
### ----------------------------
PlotMapOutput <- function(x, map, map.subid.column = 1, var.name = "", map.adj = 0, plot.legend = T, 
                          legend.pos = "right", legend.title = NULL, legend.outer = F, legend.inset = c(0, 0), 
                          col.ramp.fun = "auto", col.breaks = NULL, plot.scale = T, plot.arrow = T, 
                          par.cex = 1, par.mar = rep(0, 4) + .1, add = FALSE, restore.par = FALSE,main.title=NULL,par.mai=NULL) {
  
  # input argument checks
#  stopifnot(is.data.frame(x), dim(x)[2] == 2, class(map)=="SpatialPolygonsDataFrame", 
#            is.null(col.breaks) || is.numeric(col.breaks))
#  stopifnot(map.adj %in% c(0, .5, 1))
#  stopifnot(legend.pos %in% c("bottomright", "right", "topright", "topleft", "left", "bottomleft"))
#  if (length(col.breaks) == 1) {
#    col.breaks <- range(x[, 2], na.rm = T)
#    warning("Just one value in user-provided argument 'col.breaks', set to range of 'x[, 2]'.")
#  }
#  if (!is.null(col.breaks) && (min(col.breaks) > min(x[, 2], na.rm = T) || max(col.breaks) < max(x[, 2], na.rm = T))) {
#    warning("Range of user-provided argument 'col.breaks' does not cover range of 'x[, 2]. 
#            Areas outside range will be excluded from plot.")
#  }
  
  # add y to legend inset if not provided by user
  if (length(legend.inset) == 1) {
    legend.inset[2] <- 0
  }
  
  # save current state of par() variables which are altered below, for restoring on function exit
  par.mar0 <- par("mar")
  par.xaxs <- par("xaxs")
  par.yaxs <- par("yaxs")
  par.lend <- par("lend")
  par.xpd <- par("xpd")
  par.cex0 <- par("cex")
  par.mai0 <- par("mai")
  if(is.null(par.mai)){
    par.mai = par.mai0
  }
  if (restore.par) {
    on.exit(par(mar = par.mar0, xaxs = par.xaxs, yaxs = par.yaxs, lend = par.lend, xpd = par.xpd, cex = par.cex0, mai = par.mai0))
  }
  
  # data preparation and conditional assignment of color ramp functions and break point vectors 
  # to internal variables crfun and cbrks
  
  if (is.function(col.ramp.fun)) {
    # Case 1: a color ramp palette function is supplied
    crfun <- col.ramp.fun
    if (!is.null(col.breaks)) {
      cbrks <- col.breaks
    } else {
      cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
    }
  } else if (is.character(col.ramp.fun)) {
    # Case 2: no color ramp palette function is supplied and one of the predefined is requested
    # First treat the specific palette function strings, then "auto" requests, and last error handling for all other strings.
    # Specific palettes get a generic class break points if not provided with another by the user
    # THIS CODE IS REPETITIVE, COULD BE STREAMLINED BY BREAKING OUT cbrks ASSIGNMENT
    if (col.ramp.fun == "ColNitr") {
      crfun <- ColNitr
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else if (col.ramp.fun == "ColPhos") {
      crfun <- ColPhos
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else if (col.ramp.fun == "ColTemp") {
      crfun <- ColTemp
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else if (col.ramp.fun == "ColPrec") {
      crfun <- ColPrec
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else if (col.ramp.fun == "ColQ") {
      crfun <- ColQ
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else if (col.ramp.fun == "ColDiffTemp") {
      crfun <- ColDiffTemp
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        cbrks <- c(ifelse(min(x[,2]) < 7.5, min(x[,2]) - 1, 30), -7.5, -5, -2.5, 1, 0, 1, 2.5, 5, 7.5, ifelse(max(x[,2]) > 7.5, max(x[,2]) + 1, 30))
        #cbrks <- quantile(x[, 2], probs = seq(0, 1, .1))
      }
      
    } else if (col.ramp.fun == "ColDiffGeneric") {
      crfun <- ColDiffGeneric
      if (!is.null(col.breaks)) {
        cbrks <- col.breaks
      } else {
        # create a break point sequence which is centered around zero, with class withs based on equal intervals of the log-scaled
        # variable distribution
        cbrks <- c(rev(exp(seq(0, log(max(abs(range(x[,2]))) + 1), length.out = 5)) * -1), exp(seq(0, log(max(abs(range(x[,2]))) + 1), length.out = 5)))
        #cbrks <- quantile(x[, 2], probs = seq(0, 1, .1))
      }
      
    } else if (col.ramp.fun == "auto") {
      # Here follows a limited set of pre-defined color ramps and break point vectors for select HYPE variables, and
      # at the end a generic "catch the rest" treatment for undefined variables
      if (toupper(var.name) == "CCTN") {
        crfun <- ColNitr
        cbrks <- c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, ifelse(max(x[,2]) > 5000, max(x[,2]) + 1, 10000))
        if (is.null(legend.title)) {
          legend.title <- expression(paste("Total N (", mu, "g l"^"-1", ")"))
        }
      } else if (toupper(var.name) == "CCTP") {
        crfun <- ColPhos
        cbrks <- c(0, 5, 10, 25, 50, 100, 150, 200, 250, ifelse(max(x[,2]) > 250, max(x[,2]) + 1, 1000))
        if (is.null(legend.title)) {
          legend.title <- expression(paste("Total P (", mu, "g l"^"-1", ")"))
        }
      } else if (toupper(var.name) == "COUT") {
        crfun <- ColQ
        cbrks <- c(0, .5, 1, 5, 10, 50, 100, 500, ifelse(max(x[,2]) > 500, max(x[,2]) + 1, 2000))
        if (is.null(legend.title)) {
          legend.title <- expression(paste("Q (m"^3, "s"^"-1", ")"))
        }
      } else if (toupper(var.name) == "TEMP") {
        crfun <- ColTemp
        cbrks <- c(ifelse(min(x[,2]) < -7.5, min(x[,2]) - 1, -30), -7.5, -5, -2.5, 1, 0, 1, 2.5, 5, 7.5, ifelse(max(x[,2]) > 7.5, max(x[,2]) + 1, 30))
        if (is.null(legend.title)) {
          legend.title <- expression(paste("Air Temp. ("*degree, "C)"))
        }
      } else {
        crfun <- ColDiffGeneric
        cbrks <- quantile(x[, 2], probs = seq(0, 1, .1), na.rm = T)
      }
    } else {
      # Error treatment for all other strings
      stop("Invalid 'col.ramp.fun' argument. Neither a function nor a recognised character string.")
    }
  } else {
    # Error treatment for all other types of user input
    stop("Invalid 'col.ramp.fun' argument. Neither a function nor a character string.")
  }
  
  
  # in variables with large numbers of "0" values, the lower 10%-percentiles can be repeatedly "0", which leads to an error with cut,
  # so cbrks is shortened to unique values (this affects only the automatic quantile-based breaks)
  # if just one value remains (or was requested by user), replace crbks by minmax-based range (this also resolves unexpected behaviour
  # with single-value cbrks in 'cut' below).
  cbrks <- unique(cbrks)
  if (length(cbrks) == 1) {
    cbrks <- range(cbrks) + c(-1, 1)
  }
  # discretise the modeled values in x into classed groups, add to x as new column (of type factor)
  x[, 3] <- cut(x[, 2], breaks = cbrks, include.lowest = T)
  # replace the factor levels with color codes using the color ramp function assigned above
  levels(x[, 3]) <- crfun(length(cbrks) - 1)
  # convert to character to make it conform to plotting requirements below
  x[, 3] <- as.character(x[, 3])
  # give it a name
  names(x)[3] <- "color"
  
  # add x to subid map table (in data slot, indicated by @), merge by SUBID
  map@data <- data.frame(map@data, x[match(map@data[, map.subid.column], x[,1]),])
  
  # update legend title if none was provided by user or "auto" selection
  if (is.null(legend.title)) {
    legend.title <- toupper(var.name)
  }
  #x11(width = 4.5, height = 9)
  # par settings: lend set to square line endings because the legend below works with very thick lines 
  # instead of boxes (a box size limitation work-around); xpd set to allow for plotting a legend on the margins
  if (!add) {
    par(mar = par.mar, xaxs = "i", yaxs = "i", lend = 1, xpd = T, cex = par.cex, mai = par.mai)
    #plot.window(xlim = 0:1, ylim = 0:1)
    frame()
  } else {
    par(lend = 1, xpd = T, cex = par.cex)
  }
  
  
  ## the positioning of all plot elements works with three scales for the device's plot region: 
  ## inches, fraction, and map coordinates
  
  # plot width (inches)
  p.in.wd <- par("pin")[1]
  
  # legend position (fraction if 'add' is FALSE, otherwise already in map coordinates) 
  leg.fr.pos <- legend(legend.pos, legend = rep(NA, length(cbrks) - 1),
                       col = crfun(length(cbrks) - 1), lty = 1, lwd = 14,  bty = "n", title = legend.title, plot = F)
  # legend width (fraction if 'add' is FALSE, otherwise already in map coordinates) 
  leg.fr.wd <- leg.fr.pos$rect$w
  # legend box element height (fraction), with workaround for single-class maps
  if (length(leg.fr.pos$text$y) == 1) {
    te <- legend(legend.pos, legend = rep(NA, length(cbrks)),
                 col = crfun(length(cbrks)), lty = 1, lwd = 14,  bty = "n", title = legend.title, plot = F)
    legbx.fr.ht <- diff(c(te$text$y[length(cbrks)], te$text$y[length(cbrks) - 1]))
  } else {
    legbx.fr.ht <- diff(c(leg.fr.pos$text$y[length(cbrks) - 1], leg.fr.pos$text$y[length(cbrks) - 2]))
  }
  
  
  ## prepare legend annotation
  
  # formatted annotation text (to be placed between legend boxes which is not possible with legend() directly)
  ann.txt <- signif(cbrks, digits = 2)
  # conditional: remove outer break points
  if (!legend.outer) {
    ann.txt[c(1, length(ann.txt))] <- ""
  }
  # annotation width (inches)
  ann.in.wd <- max(strwidth(ann.txt, "inches"))
  # legend inset required to accomodate text annotation, and scalebar (always below legend)
  leg.inset <- c(ann.in.wd/p.in.wd, if(legend.pos %in% c("bottomright", "bottomleft")) {0.1} else {0})
  
  # conditional on legend placement side (legend annotation always right of color boxes)
  if (legend.pos %in% c("bottomright", "right", "topright")) {
    
    # update legend inset
    legend.inset <- legend.inset + leg.inset
    ## annotation positions (fraction if 'add' is FALSE, otherwise already in map coordinates)
    # inset scaling factor, used if 'add' is TRUE, otherwise 1 (explicitly because usr does not get updated directly when set)
    if (add) {
      f.inset.x <- par("usr")[2] - par("usr")[1]
      f.inset.y <- par("usr")[4] - par("usr")[3]
    } else {
      f.inset.x <- 1
      f.inset.y <- 1
    }
    ann.fr.x <- rep(leg.fr.pos$text$x[1], length(ann.txt)) - legend.inset[1] * f.inset.x - 0.01
    if (legend.pos == "bottomright") {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks))) + legend.inset[2] * f.inset.y
    } else if (legend.pos == "right") {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks)))
    } else {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks))) - legend.inset[2] * f.inset.y
    }
    
  } else {
    # left side legend
    # update legend inset
    legend.inset[2] <- legend.inset[2] + leg.inset[2]
    ## annotation positions (fraction if 'add' is FALSE, otherwise already in map coordinates)
    # inset scaling factor, used if 'add' is TRUE, otherwise 1 (explicitly because usr does not get updated directly when set)
    if (add) {
      f.inset.x <- par("usr")[2] - par("usr")[1]
      f.inset.y <- par("usr")[4] - par("usr")[3]
    } else {
      f.inset.x <- 1
      f.inset.y <- 1
    }
    ann.fr.x <- rep(leg.fr.pos$text$x[1], length(ann.txt)) + legend.inset[1] * f.inset.x - 0.01
    if (legend.pos == "bottomleft") {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks))) + legend.inset[2] * f.inset.y
    } else if (legend.pos == "left") {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks)))
    } else {
      ann.fr.y <- rev(seq(from = leg.fr.pos$text$y[length(cbrks) - 1] - legbx.fr.ht/2, by = legbx.fr.ht, length.out = length(cbrks))) - legend.inset[2] * f.inset.y
    }
  }
  
  
  ## calculate coordinates for map positioning
  
  # map coordinates,unprojected maps need a workaround with dummy map to calculate map side ratio
  if (is.projected(map)) {
    bbx <- bbox(map)
    # map side ratio (h/w)
    msr <- apply(bbx, 1, diff)[2] / apply(bbx, 1, diff)[1]
    # plot area side ratio (h/w)
    psr <- par("pin")[2] / par("pin")[1]
  } else {
    bbx <- bbox(map)
    # set user coordinates using a dummy plot (no fast way with Spatial polygons plot, therefore construct with SpatialPoints map)
    if(!is.null(main.title)){
      plot(SpatialPoints(coordinates(map), proj4string = CRS(proj4string(map))), col = NULL, xlim = bbx[1, ], ylim = bbx[2, ],main=main.title)
    }else{
      plot(SpatialPoints(coordinates(map), proj4string = CRS(proj4string(map))), col = NULL, xlim = bbx[1, ], ylim = bbx[2, ])
    }
    # create a map side ratio based on the device region in user coordinates and the map bounding box
    p.range.x <- diff(par("usr")[1:2])
    p.range.y <- diff(par("usr")[3:4])
    m.range.x <- diff(bbox(map)[1, ])
    m.range.y <- diff(bbox(map)[2, ])
    # map side ratio (h/w)
    msr <- m.range.y / m.range.x
    # plot area side ratio (h/w)
    psr <- p.range.y / p.range.x
  }
  
  
  # define plot limits, depending on (a) map and plot ratios (plot will be centered if left to automatic) and (b) user choice
  if (msr > psr) {
    # map is smaller than plot window in x direction, map can be moved left or right
    if (map.adj == 0) {
      pylim <- as.numeric(bbx[2, ])
      pxlim <- c(bbx[1, 1], bbx[1, 1] + diff(pylim)/psr)
    } else if (map.adj == .5) {
      pylim <- as.numeric(bbx[2, ])
      pxlim <- c(mean(as.numeric(bbx[1, ])) - diff(pylim)/psr/2, mean(as.numeric(bbx[1, ])) + diff(pylim)/psr/2)
    } else {
      pylim <- as.numeric(bbx[2, ])
      pxlim <- c(bbx[1, 2] - diff(pylim)/psr, bbx[1, 2])
    }
  } else {
    # map is smaller than plot window in y direction, map can be moved up or down
    if (map.adj == 0) {
      pxlim <- as.numeric(bbx[1, ])
      pylim <- c(bbx[2, 1], bbx[2, 1] + diff(pxlim)*psr)
    } else if (map.adj == .5) {
      pxlim <- as.numeric(bbx[1, ])
      pylim <- c(mean(as.numeric(bbx[2, ])) - diff(pxlim)*psr/2, mean(as.numeric(bbx[2, ])) + diff(pxlim)*psr/2)
    } else {
      pxlim <- as.numeric(bbx[1, ])
      pylim <- c(bbx[2, 2] - diff(pxlim)*psr, bbx[2, 2])
    }
  }
  
  
  ## plot the map and add legend using the positioning information derived above
  
  # map, plot in current frame if not added because a new frame was already created above for calculating all the coordinates
  if (!add) {
    par(new = TRUE)
  }
  if(!is.null(main.title)){
    plot(map, col = map$color, border = NA, ylim = pylim, xlim = pxlim, add = add, main=main.title)
  }else{
    plot(map, col = map$color, border = NA, ylim = pylim, xlim = pxlim, add = add)
  }
  # legend
  if (plot.legend) {
    legend(legend.pos, legend = rep(NA, length(cbrks) - 1), inset = legend.inset, 
           col = crfun(length(cbrks) - 1), lty = 1, lwd = 14,  bty = "n", title = legend.title)
    # convert annotation positioning to map coordinates, only if 'add' is FALSE
    # then plot annotation text
    if (!add) {
      ann.mc.x <- ann.fr.x * diff(pxlim) + pxlim[1]
      ann.mc.y <- ann.fr.y * diff(pylim) + pylim[1]
      text(x = ann.mc.x, y = ann.mc.y, labels = ann.txt, adj = c(0, .5), cex = 0.8)
    } else {
      text(x = ann.fr.x, y = ann.fr.y, labels = ann.txt, adj = c(0, .5), cex = 0.8)
    }
  }
  
  
  ## scale position (reference point: lower left corner), also used as reference point for north arrow
  ## conditional on 'add'
  
  if (add) {
    
    # x position conditional on legend placement side
    if (legend.pos %in% c("bottomright", "right", "topright")) {
      lx <- par("usr")[2] - signif(diff(par("usr")[1:2])/4, 0) - legend.inset[1] * diff(par("usr")[1:2])
    } else {
      lx <- par("usr")[1] + (legend.inset[1] + 0.02) * diff(par("usr")[1:2])
    }
    
    # y position conditional legend placement position (leg.fr.pos here is already in map coordinates)
    if (legend.pos %in% c("bottomright", "bottomleft")) {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h + legend.inset[2]*f.inset.y/2)
    } else if (legend.pos %in% c("right", "left")) {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h + (legend.inset[2]/2 - .1) * f.inset.y)
    } else {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h - (legend.inset[2]/2 - .1) * f.inset.y)
    }
  } else {
    
    # x position conditional on legend placement side
    if (legend.pos %in% c("bottomright", "right", "topright")) {
      lx <- pxlim[2] - signif(diff(bbx[1,])/4, 0) - legend.inset[1] * diff(pxlim)
    } else {
      lx <- pxlim[1] + (legend.inset[1] + 0.02) * diff(pxlim)
    }
    
    # y position conditional legend placement position
    if (legend.pos %in% c("bottomright", "bottomleft")) {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h + legend.inset[2]/2) * diff(pylim) + pylim[1]
    } else if (legend.pos %in% c("right", "left")) {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h + legend.inset[2]/2 - .1) * diff(pylim) + pylim[1]
    } else {
      ly <- (leg.fr.pos$rect$top - leg.fr.pos$rect$h - legend.inset[2]/2 - .1) * diff(pylim) + pylim[1]
    }
  }
  
  if (plot.scale) {
    if (!is.projected(map)) {
      warning("Scale bar meaningless with un-projected maps. Set 'plot.scale = F' to remove it.")
    }
    if (!add) {
      ldistance <- signif(diff(bbx[1,])/4, 0)
    } else {
      ldistance <- signif(diff(par("usr")[1:2])/4, 0)
    }
    .Scalebar(x = lx, 
              y = ly, 
              distance = ldistance, 
              scale = 0.001, t.cex = 0.8)
  }
  
  if (plot.arrow) {
    
    if (add) {
      nlen <- diff(par("usr")[1:2])/70
      # north arrow x position conditional on side where legend is plotted
      if (legend.pos %in% c("bottomright", "right", "topright")) {
        nx <- lx - 0.02 * diff(par("usr")[1:2])
      } else {
        nx <- lx + signif(diff(par("usr")[1:2])/4, 0) + 0.055 * diff(par("usr")[1:2])
      }
    } else {
      nlen <- diff(bbx[1,])/70
      # north arrow x position conditional on side where legend is plotted
      if (legend.pos %in% c("bottomright", "right", "topright")) {
        nx <- lx - 0.02 * diff(pxlim)
      } else {
        nx <- lx + signif(diff(bbx[1,])/4, 0) + 0.055 * diff(pxlim)
      }
    }
    
    .NorthArrow(xb = nx, 
                yb = ly, 
                len = nlen, cex.lab = .8)
  }
  
  
  # invisible unless assigned: return map with added data and color codes
  invisible(map)
}


### ReadMapOutput from HYPEtools
### ----------------------------
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
# 5- function to write a world file to geotag an image file
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

# --------------------------------------------------------------------------
# 5- plot content of the map output as png-file to output dir
# --------------------------------------------------------------------------

# next lines used only for development
#   fileDir      = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test"   # folder with mapoutput files
#   plotDir      = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test"   # folder where the plot files should be written
#   mapFile      = "004_mapCOUT.txt"   # mapoutput filename
#   modelName    = "Niger-hype"   # modelName (used for plot main title)
#   rdataFile    = "D:/TEP/ModellingService/sandbox/home/dgustafsson/hype-files/test/shp.Rdata"   # path to rdata file with subbasin spatial polygons data frame
#   prefix.png   = "001"   # png filename prefix
#   cdate=as.POSIXct("1979-01-01",tz="GMT")
#   edate=as.POSIXct("1980-01-01",tz="GMT")

# echo arguments to the TEP log file'
rciop.log ("DEBUG", paste(" plot-mapoutput, fileDir: ",fileDir,sep=""), "/util/R/hypeapps-plot-mapoutput.R")
rciop.log ("DEBUG", paste(" plot-mapoutput, plotDir: ",plotDir,sep=""), "/util/R/hypeapps-plot-mapoutput.R")
rciop.log ("DEBUG", paste(" plot-mapoutput, mapFile: ",mapFile,sep=""), "/util/R/hypeapps-plot-mapoutput.R")
rciop.log ("DEBUG", paste(" plot-mapoutput, rdataFile: ",rdataFile,sep=""), "/util/R/hypeapps-plot-mapoutput.R")

# read map file
mapData = ReadMapOutput(filename = paste(fileDir,mapFile,sep="/"))

# load subbasin spatial points data frame (shapefileData)
load(rdataFile)

# # next line used only for development
#shapefileData = shp

# variable name from MapOutput filename
varName = substr(mapFile,nchar(mapFile)-7,nchar(mapFile)-4)
# Main title
mainTitle=paste(toupper(substr(modelName,1,1)),tolower(substr(modelName,2,nchar(modelName))),sep="")
# hype to HYPE or Hype to HYPE
testHYPE = regexpr(pattern ="hype",mainTitle)
if(testHYPE[1]>0){
  substr(mainTitle,testHYPE[1],testHYPE[1]+3)<-"HYPE"
}else{
  testHYPE = regexpr(pattern ="ype",mainTitle)
  if(testHYPE[1]>0){
    substr(mainTitle,testHYPE[1],testHYPE[1]+2)<-"YPE"
  }
}

if(varName=="COUT"){
  # discharge
  mainTitle=paste(mainTitle," - Mean river discharge \n(",cdate,"-",edate,")",sep="")
  varPlot="m3/s"
}else if(varName=="CPRC" | varName=="PREC"){
  # precipitation
  mainTitle=paste(mainTitle," - Mean annual precipitation \n(",cdate,"-",edate,")",sep="")
  varPlot="mm/yr"
}else if(varName=="CRUN"){
  # runoff
  mainTitle=paste(mainTitle," - Mean annual runoff \n(",cdate,"-",edate,")",sep="")
  varPlot="mm/yr"
}else if(varName=="EVAP"){
  # evapotranspiration
  mainTitle=paste(mainTitle," - Mean annual evapotraspiration \n(",cdate,"-",edate,")",sep="")
  varPlot="mm/yr"
}else if(varName=="CTMP" | varName=="TEMP"){
  # temperature
  mainTitle=paste(mainTitle," - Mean air temperature \n(",cdate,"-",edate,")",sep="")
  varPlot="degree Celcius"
}else{
  mainTitle=paste(mainTitle," - Mean ", varName,"\n(",cdate,"-",edate,")",sep=", ")
  varPlot=varName
}

# graphics file name
plotFileName = paste(plotDir,prefix.png,sep="/")
plotFileName = paste(plotFileName,"_",substr(mapFile,1,nchar(mapFile)-3),sep="")
if(pngORjpg==1){
  plotFileName = paste(plotFileName,"jpg",sep="")
}else{
  plotFileName = paste(plotFileName,"png",sep="")
}


# Plot dimensions and World file
ydiff = bbox(shapefileData)[2,2]-bbox(shapefileData)[2,1]
xdiff = bbox(shapefileData)[1,2]-bbox(shapefileData)[1,1]

width = round(xdiff/(1/60),digits=0)
height = round(width * ydiff/xdiff*1.03,digits=0)


cx = bbox(shapefileData)[1,2] - 0.5 * (bbox(shapefileData)[1,2]-bbox(shapefileData)[1,1])
cy = bbox(shapefileData)[2,2] - 0.5 * (bbox(shapefileData)[2,2]-bbox(shapefileData)[2,1])
wfres = writeWorldFile(paste(plotFileName,"w",sep=""), pxWidth=round(width,digits=0), pxHeight=height,
                       degWidth=xdiff,degHeight=ydiff, lonBasin=cx, latBasin=cy, plotPos="center")


# initiate jpeg or png file for plotting using Cairo graphics device
if(pngORjpg==1){
  CairoJPEG(filename = plotFileName, width = width, height = height, units = "px",bg = "white")
}else{
  CairoPNG(filename = plotFileName, width = width, height = height, units = "px",bg = "white")
}

# mapplot
PlotMapOutput(x = mapData, 
             map = shapefileData,
             col.ramp.fun = "ColQ",
             plot.scale = F,
             plot.arrow = F,
             legend.title = varPlot,
             legend.pos = "topleft",
             par.mar = c(0,0,0,0),par.cex = 3,par.mai=c(0,0,0,0))
# text
legend("bottomleft",legend = mainTitle,bty="n")
# Close plot
dev.off()



