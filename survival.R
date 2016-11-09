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
vi50 <- read.csv("volumeintersection50.csv")

## categorize individuals 
#  and rank by 50% UD volume overlap 
mig <- rawdata %>%
  dplyr::select(-FallVI) %>%
  rename(VI = SprVI) %>%
  right_join(vi50, by = "IndivYr") %>%
  select(-c(Sex, AnimalID.x, AnimalID.y)) %>%
  rename(VI50 = SprVI50) %>%
  transform(Rank = rank(-VI50, ties.method = "average")) 

# determine appropriate cutoff
hist(mig$VI50)
median(mig$VI50)
hist((mig$VI50)^(1/2))

migstatus <- mig %>%
  transform(Mig = ifelse(VI50 >= 0.25, "Resident",
                  ifelse(VI50 < 0.001, "Migrant", 
                         "Intermediate")))

par(mfrow=c(2,1))
hist(migstatus$VI)
hist(migstatus$VI50)

plot(VI50 ~ Rank, col = Mig, data=migstatus)
  #oh hell, everybody's intermediate... rethink
length(which(migstatus$Mig == "Migrant")) # still feel ok about this one
length(which(migstatus$Mig == "Resident")) # nope nope nope
length(which(migstatus$Mig == "Intermediate"))
