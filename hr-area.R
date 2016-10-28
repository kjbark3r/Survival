######################################################
## Migration - Volume Overlap Calculations ##
########  NSERP - Kristin Barker - June 2016  ########
######################################################

##SET WD
####Work computer, personal laptop, or external hard drive
wd_workcomp <- "C:\\Users\\kristin.barker\\Documents\\GitHub\\Survival"
wd_laptop <- "C:\\Users\\kjbark3r\\Documents\\GitHub\\Survival"
wd_external <- "E:\\Kristins\\Survival\\"

if (file.exists(wd_workcomp)) {
  setwd(wd_workcomp)
} else {
  if(file.exists(wd_laptop)) {
    setwd(wd_laptop)
  } else {
    if(file.exists(wd_external)) {
      setwd(wd_external)
    } else {
      cat("Are you SURE you got that file path right?\n")
    }
  }
}
rm(wd_workcomp, wd_laptop, wd_external)

##LOAD PACKAGES
library(sp) #for kernel centroid estimate
library(adehabitatHR) #for kernel centroid estimate
library(raster) #prob don't need this one here
library(rgdal) #for latlong/stateplane conversions
library(gsubfn)
library(maptools) #for writeSpatialShape
library(dplyr) #for joins

###########################################################################################
#SET UP DATA
###########################################################################################

#GPS DATA FROM COLLARS, PLUS MIGRATION SEASON/YEAR
locs <- read.csv("../Migration/HRoverlap/collardata-locsonly-equalsampling.csv", as.is = TRUE, header = TRUE)
locs$Date <- as.Date(locs$Date, "%Y-%m-%d")

# delineate winter and summer HRs
locs$MigHR <- ifelse(between(locs$Date, as.Date("2014-01-01"), as.Date("2014-03-15")), "Winter", 
                     ifelse(between(locs$Date, as.Date("2014-07-01"), as.Date("2014-08-31")), "Summer", 
                            ifelse(between(locs$Date, as.Date("2015-01-01"), as.Date("2015-03-15")), "Winter", 
                                   ifelse(between(locs$Date, as.Date("2015-07-01"), as.Date("2015-08-31")), "Summer",
                                          ifelse(between(locs$Date, as.Date("2016-01-01"), as.Date("2016-03-15")), "Winter",
                                                 ifelse(NA))))))
locs$IndivYr <- ifelse(locs$Date < "2015-01-01", 
                       paste(locs$AnimalID, "-14", sep=""),
                       paste(locs$AnimalID, "-15", sep=""))

#LIST OF ANIMALS TO RUN
#indivyrs <- as.data.frame(unique(locs$IndivYr))
#numelk <- nrow(indivyrs)

#DEFINE PROJECTIONS
latlong <- CRS("+init=epsg:4326")
stateplane <- CRS("+init=epsg:2818")

###########################################################################################
#CALCULATE HR AREA PER INDIVIDUAL
###########################################################################################

win <- subset(locs, MigHR == "Winter")
xy <- data.frame("x"=win$Long,"y"=win$Lat)
xy.spdf.ll <- SpatialPointsDataFrame(xy, win, proj4string = latlong)
xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
kuds <- kernelUD(xy.spdf.sp[,23], h = "href", same4all = FALSE)
hrs <- getverticeshr(kuds)
win.a <- sapply(slot(hrs, "polygons"), slot, "area")
median(win.a)

sum <- subset(locs, MigHR == "Summer")
xy <- data.frame("x"=sum$Long,"y"=sum$Lat)
xy.spdf.ll <- SpatialPointsDataFrame(xy, sum, proj4string = latlong)
xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
kuds <- kernelUD(xy.spdf.sp[,23], h = "href", same4all = FALSE)
hrs <- getverticeshr(kuds)
sum.a <- sapply(slot(hrs, "polygons"), slot, "area")
median(sum.a)

median(win.a)/median(sum.a)
# well, that was a bust
# turns out summer HR is larger
# and ratio of winter:summer is 0.65