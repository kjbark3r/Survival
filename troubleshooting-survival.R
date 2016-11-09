##############################
## misc code related to
## survival consequences of varying mig behavs
### kjb nov 2016
##############################

## SETUP #### 

## WD
wd_workcomp <- "C:\\Users\\kristin.barker\\Documents\\GitHub\\Survival"
wd_laptop <- "C:\\Users\\kjbark3r\\Documents\\GitHub\\Survival"
if (file.exists(wd_workcomp)) {
  setwd(wd_workcomp)
} else {
  if(file.exists(wd_laptop)) {
    setwd(wd_laptop)
    } else {
      cat("Are you SURE you got that file path right?\n")
    }
}
rm(wd_workcomp, wd_laptop)

## libraries
library(dplyr)
library(adehabitatHR) # needed?
library(raster)

## read in data (already subsetted to females only)
rawdata <- read.csv("../Migration/HRoverlap/volumeintersection.csv")
phrdata <- read.csv("phr.csv")
lookdata <- read.csv("migstatus-prelimlook.csv")
lookdata <- lookdata %>% select(c(AnimalID, Status))

## DISCRETIZING ####

## categorize individuals (and rank, for funzies) 
mig <- rawdata %>%
  right_join(phrdata, by = "IndivYr") %>%
  right_join(lookdata, by = "AnimalID") %>%
  select(-FallVI) %>%
  rename(VI = SprVI) %>%
  transform(Rank = rank(-VI, ties.method = "random")) %>%
  transform(Mig = ifelse(PHR >= 0.5, "Resident",
                  ifelse(PHR == 0, "Migrant", 
                         "Intermediate"))) %>%
  transform(PHRrank = rank(-PHR, ties.method="random"))
            
## visualize
plot(VI ~ Rank, col = Mig, data=mig)
plot(PHR ~ PHRrank, col = Mig, data=mig)
  #oh hell, everybody's intermediate... rethink
hist(mig$VI)

length(which(mig$Mig == "Migrant")) # feel good about this one
length(which(mig$Mig == "Resident")) # nope nope nope
length(which(mig$Mig == "Intermediate"))

## use ratio of avg summer:winter HR area
## to determine residency proportion cutoff

# read in winter and summer HR tifs from previous analysis
  # time frames of these are slightly different than what you're doing now
  # i don't think it matters, but can rerun with new dates if needed
kde.w <- list.files("../zOld/Migration", pattern = ".*-[Ww].*\\.tif$")
kde.s <- list.files("../zOld/Migration", pattern = ".*-[Ss].*\\.tif$")

test <- raster("../zOld/Migration/KDE140040-S15.tif")
# test2 <- getverticeshr(test) # newp
plot(test) # ooooh pretty
# area(test) # newp
# test2 <- kernelUD(test) # newp
# test2 <- rasterToPolygons(test, fun = function(x) {x < 0.95*max(x)}) # r-splosion
# test2 <- SpatialPointsDataFrame(test) # newp
# test2 <- getverticeshr("../zOld/Migration/KDE140040-S15.tif") # newp
# test2 <- kernel.area(SOMETHING, percent = 95) # nevermind, x has to be estUD
# feck.. fine, i'll start over with the locs...

pts.w <- list.files("../zOld/Migration", pattern = "Pts.*-[Ww].*\\.shp$")
pts.s <- list.files("../zOld/Migration", pattern = "Pts.*-[Ss].*\\.shp$")


## VI on 50% (core area) rather than 95% UD ####

#GPS DATA FROM COLLARS, PLUS MIGRATION SEASON/YEAR
locs <- read.csv("../ElkDatabase/collardata-locsonly-equalsampling.csv", as.is = TRUE, header = TRUE)
locs$Date <- as.Date(locs$Date, "%Y-%m-%d")

# summer = july 15 - aug 31 (match biomass sampling/gdm estimation time pd)
# winter = starts feb 26 2014 bc day after last capture
  # other yrs set to be same length of time as summer timeframe
locs$MigHR <- ifelse(between(locs$Date, as.Date("2014-02-26"), as.Date("2014-03-31")), "Winter 2014", 
                     ifelse(between(locs$Date, as.Date("2014-07-15"), as.Date("2014-08-31")), "Summer 2014", 
                            ifelse(between(locs$Date, as.Date("2015-02-15"), as.Date("2015-03-31")), "Winter 2015", 
                                   ifelse(between(locs$Date, as.Date("2015-07-15"), as.Date("2015-08-31")), "Summer 2015",
                                          ifelse(between(locs$Date, as.Date("2016-02-15"), as.Date("2016-03-31")), "Winter 2016",
                                                 ifelse(NA))))))

# projections
latlong <- CRS("+init=epsg:4326")
stateplane <- CRS("+init=epsg:2818")

# lists of animals per searon
list.spr14 <- read.csv("../Migration/HRoverlap/spr14.csv", header = TRUE)
numelk.spr14 <- nrow(list.spr14)
list.spr15 <- read.csv("../Migration/HRoverlap/spr15.csv", header = TRUE)
numelk.spr15 <- nrow(list.spr15)

######################
#SPRING 2014 MIGRATION

spr14 <- data.frame(matrix(ncol = 2, nrow = numelk.spr14)) #create df wo NAs
colnames(spr14) <- c("AnimalID", "SprVI50")

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
  vol <- kerneloverlap(xy.spdf.sp[,20], method = "VI", percent = 50, conditional = FALSE)
  
  #store results
  spr14[[i,1]] <- elk
  spr14[[i,2]] <- vol[2,1]
}    
spr14$IndivYr <- paste(spr14$AnimalID, "-14", sep="")

######################
#SPRING 2015 MIGRATION

spr15 <- data.frame(matrix(ncol = 2, nrow = numelk.spr15)) #create df wo NAs
colnames(spr15) <- c("AnimalID", "SprVI50")

for(i in 1:numelk.spr15) {
  elk <- list.spr15[i,]
  
  #subset individual and seasonal locations
  temp_dat_spr15 <- subset(locs, AnimalID == elk) 
  temp_dat_spr15 <- subset(temp_dat_spr15, MigHR == "Winter 2015" | MigHR == "Summer 2015")
  
  #Get xy points, write points to dataframe, to spatial data frame, to stateplane projection
  xy <- data.frame("x"=temp_dat_spr15$Long,"y"=temp_dat_spr15$Lat)
  xy.spdf.ll <- SpatialPointsDataFrame(xy, temp_dat_spr15, proj4string = latlong)
  xy.spdf.sp <- spTransform(xy.spdf.ll,stateplane)
  
  #calculate area overlap and volume intersection
  vol <- kerneloverlap(xy.spdf.sp[,20], method = "VI", percent = 50, conditional = FALSE)
  
  #store results
  spr15[[i,1]] <- elk
  spr15[[i,2]] <- vol[2,1]
}    
spr15$IndivYr <- paste(spr15$AnimalID, "-15", sep="")

sex <- distinct(select(locs, AnimalID, Sex))
vi50 <- bind_rows(spr14, spr15) %>%
  left_join(sex, by = "AnimalID") %>%
  filter(Sex == "Female")

write.csv(vi50, file = "volumeintersection50.csv", row.names=F)

######################################################################################
###### MESSUPS AND DELETIONS ####
#######################################################################################

# can apply work on a list? 
  # no, have to make it df
testfcn <- function(x) {paste(x, "BLARG", sep="")}
testfcn("hey")
kde.w <- as.data.frame(list.files("../zOld/Migration", pattern = ".*-[Ww].*\\.tif$"))
apply(kde.w, 2, testfcn)