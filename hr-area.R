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
#SET UP DATA  ####
###########################################################################################

#GPS DATA FROM COLLARS, PLUS MIGRATION SEASON/YEAR
locs <- read.csv("../ElkDatabase/collardata-locsonly-equalsampling.csv", as.is = TRUE, header = TRUE)
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
#CALCULATE HR AREA PER INDIVIDUAL ####
###########################################################################################

win <- subset(locs, MigHR == "Winter")
xy <- data.frame("x"=win$Long,"y"=win$Lat)
xy.spdf.ll <- SpatialPointsDataFrame(xy, win, proj4string = latlong)
xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
kuds <- kernelUD(xy.spdf.sp[,23], h = "href", same4all = FALSE)
hrs <- getverticeshr(kuds)
win.a <- sapply(slot(hrs, "polygons"), slot, "area")
mean(win.a); median(win.a)

sum <- subset(locs, MigHR == "Summer")
xy <- data.frame("x"=sum$Long,"y"=sum$Lat)
xy.spdf.ll <- SpatialPointsDataFrame(xy, sum, proj4string = latlong)
xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
kuds <- kernelUD(xy.spdf.sp[,23], h = "href", same4all = FALSE)
hrs <- getverticeshr(kuds)
sum.a <- sapply(slot(hrs, "polygons"), slot, "area")
mean(sum.a); median(sum.a)

median(win.a)/median(sum.a)
# well, that was a bust
# turns out summer HR is larger
# and ratio of winter:summer is 0.65

###########################################################################################
#DOUBLE-CHECK YOU WERE CONSISTENT IN DOING OVERLAP OF *WINTER* BY SUMMER    ####
###########################################################################################
###########################################################################################
#tweaked original overlapcalc code
###########################################################################################
# spoiler alert: yes
# also that didn't matter and was a waste of time to look at, oops
##########################

#GPS DATA FROM COLLARS, PLUS MIGRATION SEASON/YEAR
locs <- read.csv("../Migration/HRoverlap/collardata-locsonly-equalsampling.csv", as.is = TRUE, header = TRUE)
locs$Date <- as.Date(locs$Date, "%Y-%m-%d")

# delineate winter and summer HRs
locs$MigHR <- ifelse(between(locs$Date, as.Date("2014-01-01"), as.Date("2014-03-15")), "Winter 2014", 
                     ifelse(between(locs$Date, as.Date("2014-07-01"), as.Date("2014-08-31")), "Summer 2014", 
                            ifelse(between(locs$Date, as.Date("2015-01-01"), as.Date("2015-03-15")), "Winter 2015", 
                                   ifelse(between(locs$Date, as.Date("2015-07-01"), as.Date("2015-08-31")), "Summer 2015",
                                          ifelse(between(locs$Date, as.Date("2016-01-01"), as.Date("2016-03-15")), "Winter 2016",
                                                 ifelse(NA))))))
locs$IndivYr <- ifelse(locs$Date < "2015-01-01", 
                       paste(locs$AnimalID, "-14", sep=""),
                       paste(locs$AnimalID, "-15", sep=""))

######################
#SPRING 2014 MIGRATION

  #subset individual and seasonal locations
  temp_dat_spr14 <- subset(locs, MigHR == "Winter 2014" | MigHR == "Summer 2014")
  
  #Get xy points, write to dataframe, to spatial data frame, to stateplane projection
  xy <- data.frame("x"=temp_dat_spr14$Long,"y"=temp_dat_spr14$Lat)
  xy.spdf.ll <- SpatialPointsDataFrame(xy, temp_dat_spr14, proj4string = latlong)
  xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
  
  #calculate area overlap and volume intersection
  vol.spr <- kerneloverlap(xy.spdf.sp[,22], method = "VI", percent = 95, conditional = FALSE)
  

######################
#FALL 2014 MIGRATION

  #subset individual and seasonal locations
  temp_dat_fall14 <- subset(locs, MigHR == "Summer 2014" | MigHR == "Winter 2015")
  
  #Get xy points, write to dataframe, to spatial data frame, to stateplane projection
  xy <- data.frame("x"=temp_dat_fall14$Long,"y"=temp_dat_fall14$Lat)
  xy.spdf.ll <- SpatialPointsDataFrame(xy, temp_dat_fall14, proj4string = latlong)
  xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
  
  #calculate area overlap and volume intersection
  vol.fall <- kerneloverlap(xy.spdf.sp[,22], method = "VI", percent = 95, conditional = FALSE)

###########################################################################################
# CHECKING OUT PHR, PROBABILITY OF OVERLAP    ####
###########################################################################################

# testing storing new data in a spatialpointsdf
xy.spdf.sp@data$test <- "Woot"
test <- as.data.frame(xy.spdf.sp)
# damn i'm good
  # at simplistic tasks
      
#GPS DATA FROM COLLARS, PLUS MIGRATION SEASON/YEAR
locs <- read.csv("../Migration/HRoverlap/collardata-locsonly-equalsampling.csv", as.is = TRUE, header = TRUE)
locs$Date <- as.Date(locs$Date, "%Y-%m-%d")

# delineate winter and summer HRs
locs$MigHR <- ifelse(between(locs$Date, as.Date("2014-01-01"), as.Date("2014-03-15")), "Winter 2014", 
                     ifelse(between(locs$Date, as.Date("2014-07-01"), as.Date("2014-08-31")), "Summer 2014", 
                            ifelse(between(locs$Date, as.Date("2015-01-01"), as.Date("2015-03-15")), "Winter 2015", 
                                   ifelse(between(locs$Date, as.Date("2015-07-01"), as.Date("2015-08-31")), "Summer 2015",
                                          ifelse(between(locs$Date, as.Date("2016-01-01"), as.Date("2016-03-15")), "Winter 2016",
                                                 ifelse(NA))))))
locs$IndivYr <- ifelse(locs$Date < "2015-01-01", 
                       paste(locs$AnimalID, "-14", sep=""),
                       paste(locs$AnimalID, "-15", sep=""))


#DEFINE PROJECTIONS
latlong <- CRS("+init=epsg:4326")
stateplane <- CRS("+init=epsg:2818")

#LISTS OF ANIMALS TO RUN
#Because code to automate this takes too long to run on my subpar computer
list.spr14 <- read.csv("../Migration/HRoverlap/spr14.csv", header = TRUE)
numelk.spr14 <- nrow(list.spr14)
list.fall14 <- read.csv("../Migration/HRoverlap/fall14.csv", header = TRUE)
numelk.fall14 <- nrow(list.fall14)
list.spr15 <- read.csv("../Migration/HRoverlap/spr15.csv", header = TRUE)
numelk.spr15 <- nrow(list.spr15)
list.fall15 <- read.csv("../Migration/HRoverlap/fall15.csv", header = TRUE)
numelk.fall15 <- nrow(list.fall15)


######################
#SPRING 2014 MIGRATION

  #set up df for storing
  spr14 <- data.frame(matrix(ncol = 2, nrow = numelk.spr14)) #create df wo NAs
  colnames(spr14) <- c("AnimalID", "PHR")

for(i in 1:numelk.spr14) {
  elk <- list.spr14[i,]
  
  #subset individual and seasonal locations
  temp_dat_spr14 <- subset(locs, AnimalID == elk) 
  temp_dat_spr14 <- subset(temp_dat_spr14, MigHR == "Winter 2014" | MigHR == "Summer 2014")

  #Get xy points, write to dataframe, to spatial data frame, to stateplane projection
  xy <- data.frame("x"=temp_dat_spr14$Long,"y"=temp_dat_spr14$Lat)
  xy.spdf.ll <- SpatialPointsDataFrame(xy, temp_dat_spr14, proj4string = latlong)
  xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
  
  #calculate area overlap and volume intersection
  phr <- kerneloverlap(xy.spdf.sp[,22], method = "PHR", percent = 95, conditional = FALSE)
  
spr14[i,1] <- elk
spr14[i,2] <- phr[2,1]
}
spr14 <- spr14 %>%
  mutate(IndivYr = paste(AnimalID, "-14", sep="")) %>%
  select(-AnimalID)

######################
#SPRING 2015 MIGRATION

  #set up df for storing
  spr15 <- data.frame(matrix(ncol = 2, nrow = numelk.spr15)) #create df wo NAs
  colnames(spr15) <- c("AnimalID", "PHR")

for(i in 1:numelk.spr15) {
  elk <- list.spr15[i,]
  
  #subset individual and seasonal locations
  temp_dat_spr15 <- subset(locs, AnimalID == elk) 
  temp_dat_spr15 <- subset(temp_dat_spr15, MigHR == "Winter 2015" | MigHR == "Summer 2015")

  #Get xy points, write to dataframe, to spatial data frame, to stateplane projection
  xy <- data.frame("x"=temp_dat_spr15$Long,"y"=temp_dat_spr15$Lat)
  xy.spdf.ll <- SpatialPointsDataFrame(xy, temp_dat_spr15, proj4string = latlong)
  xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
  
  #calculate area overlap and volume intersection
  phr <- kerneloverlap(xy.spdf.sp[,22], method = "PHR", percent = 95, conditional = FALSE)
  
spr15[i,1] <- elk
spr15[i,2] <- phr[2,1]
}
spr15 <- spr15 %>%
  mutate(IndivYr = paste(AnimalID, "-15", sep="")) %>%
  select(-AnimalID)


######################
#save and export

phrdata <- bind_rows(spr14, spr15)
write.csv(phrdata, file="phr.csv", row.names=FALSE)
