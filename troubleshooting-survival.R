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



## DISCRETIZING ####

## categorize individuals (and rank, for funzies) 
mig <- rawdata %>%
  select(-FallVI) %>%
  rename(VI = SprVI) %>%
  transform(Rank = rank(-VI, ties.method = "random")) %>%
  transform(Mig = ifelse(VI >= 0.5, "Resident",
                  ifelse(VI == 0, "Migrant", 
                         "Intermediate")))
            
## visualize
plot(VI ~ Rank, col = Mig, data=mig)
  #oh hell, everybody's intermediate... rethink

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


# next steps:
  # read in above as tifs
  # calculate average area of each group
  # calculate proportion diff

######################################################################################
###### MESSUPS AND DELETIONS ####
#######################################################################################

# can apply work on a list? 
  # no, have to make it df
testfcn <- function(x) {paste(x, "BLARG", sep="")}
testfcn("hey")
kde.w <- as.data.frame(list.files("../zOld/Migration", pattern = ".*-[Ww].*\\.tif$"))
apply(kde.w, 2, testfcn)