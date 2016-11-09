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

## read in data; combine and clean
vi95 <- read.csv("../Migration/HRoverlap/volumeintersection.csv")
vi50 <- read.csv("../Migration/HRoverlap/volumeintersection50.csv")
look <- read.csv("../Migration/MethodComparison/migstatus-prelimlook.csv") %>%
  select(AnimalID, Status)

mig <- vi95 %>%
  rename(VI95 = SprVI) %>%
  select(-c(AnimalID, Sex)) %>% 
  left_join(vi50, by = "IndivYr") %>%
  rename(VI50 = SprVI) %>%
  select(IndivYr, AnimalID, VI95, VI50) %>%
  left_join(look, by = "AnimalID") %>%
  rename(Look = Status)


## categorize individuals into migratory bins

# migrant
length(which(mig$Look == "Migrant"))
length(which(mig$VI95 == 0)); length(which(mig$VI50 == 0))
  # define migrant using 95%, not 50%

#resident
length(which(mig$Look == "Resident"))
length(which(mig$VI95 > 0.1)); length(which(mig$VI95 > 0.2))
  # ok,somewhere in the ballpark or 10-20% VI
length(which(mig$VI50 > 0))
  # ooh i like this idea
  # resident is anyone whose CORES (50% UD) overlap at all
  # migrant is anyone whose HRs (95% UD) never overlap
  # intermediate is everyone else


mig <- transform(mig, 
       MigStatus = ifelse(VI50 > 0, "Resident",
                   ifelse(VI95 == 0, "Migrant",
                          "Intermediate")))

length(which(mig$MigStatus == "Migrant")) 
length(which(mig$MigStatus == "Resident")) 
length(which(mig$MigStatus == "Intermediate"))
# me likey

write.csv(mig, file = "migstatus.csv", row.names=FALSE)
