######################################################
#### SURVIVAL consequences of migratory behaviors ####
########  NSERP - Kristin Barker - Nov 2016   ########
######################################################

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

## read in data (already subsetted to females only)
rawdata <- read.csv("../Migration/HRoverlap/volumeintersection.csv")

## categorize individuals (and rank, for funzies) 
mig <- rawdata %>%
  dplyr::select(-FallVI) %>%
  rename(VI = SprVI) %>%
  transform(Rank = rank(-VI, ties.method = "average")) %>%
  transform(Mig = ifelse(VI >= 0.35, "Resident",
                  ifelse(VI == 0, "Migrant", 
                         "Intermediate")))
            
## visualize
plot(VI ~ Rank, col = Mig, data=mig)
  #oh hell, everybody's intermediate... rethink
length(which(mig$Mig == "Migrant")) # still feel ok about this one
length(which(mig$Mig == "Resident")) # nope nope nope
length(which(mig$Mig == "Intermediate"))
