####Survival Analysis for 2014-5 Sapphire adult elk
require(lubridate)
require(tidyr)
require(dplyr)
require(survival)
require(ggplot2)
require(survMisc)
require(data.table)
require(eha)
require(RODBC)

#setwd("C:/Users/CF2752/Documents/Sapphire Elk/Progress Reports/FinalReport/Adult Survival")
setwd("E:\\NSERP_DeVoe\\Sapphire Elk\\Data\\SurvivalAnalysis\\Adult Survival")

channel <- odbcDriverConnect("Driver={Microsoft Access Driver (*.mdb, *.accdb)};
                                   dbq=E:/NSERP_DeVoe/Sapphire Elk/Data/GPSdata/SapphireElkProject_ElkDatabase.accdb")

#####################################################################################
#Prep the data - Data from Mortality table in GPS database
#####################################################################################

###Read in and format file of survival information/mortality table (1 row per individual)
adultelk <- sqlQuery(channel, paste("select * from MortalityInfo"))
adultelk <- adultelk %>%
  rename(AnimalID = `Animal ID`) %>%
  select(AnimalID, Gender, AgeCapture, CaptureDate, Fate, Cause1, MortDate) %>%
  mutate(Gender = factor(trimws(Gender)), #trim off space present in gender column ("Female ")
         MortDate = mdy(MortDate),  # convert dates
         CaptureDate = as.Date(CaptureDate, format = "%Y-%m-%d")) %>%
  filter(Fate != "Unknown") # remove unknowns

#####################################################################################
#Formatting/structuring data for the biological year
#####################################################################################

#Enter bio year begin dates
bioYear1 <- mdy("06/01/2014") # June 1, 2014 - May 31, 2015
bioYear2 <- mdy("06/01/2015") # June 1, 2015 - May 31, 2016
bioYear3 <- mdy("06/01/2016") # June 1, 2016 +  # this category will be filtered out, but needs to be used to properly code exit timing
lastDay <- as.numeric(mdy("02/19/2016") - bioYear2) # Last day animals checked (as far as I know): February 19, 2016 (Day 263 from beginning of bioYear2)

#Adding rows for each individual surviving into subsequent bioyears, first create column for each bioYear
adultelk <- adultelk %>% 
  mutate(bioYear1 = ifelse(CaptureDate < bioYear2, "2014-2015", NA), # in this dataset, every elk part of bioYear1
         bioYear2 = ifelse(is.na(MortDate) | MortDate >= bioYear2, "2015-2016", NA), # only for those that made it to bioYear2
         bioYear3 = ifelse(is.na(MortDate) | MortDate >= bioYear3, "2016-2017", NA)) # only for those that made it to bioYear3

#Using gather to create rows for each bioYear surviving to 
adultelk <- gather(adultelk, bioYear, bioYear.yrs, bioYear1:bioYear3, na.rm = TRUE)

#Calculate entry and exit days from biological year (Day 1 = June 1 , 2014, Day 1 = June 1, 2015)
adultelk <- adultelk %>%
  mutate(enter.bio = ifelse(CaptureDate <= bioYear1, 1, # captured before bioYear 1 (and survived to begining of bioYear1), start day 1
                            ifelse(CaptureDate > bioYear1 & CaptureDate < bioYear2 & bioYear == "bioYear1", CaptureDate - bioYear1, # captured within bioYear1, capture date minus bioYear1
                                   ifelse(CaptureDate > bioYear2, CaptureDate - bioYear2, 1)))) %>% # captured within bioYear2, capture date minus bioYear2, otherwise enterred into bioYear2 on day 1
         #currently, exit.bio only calculating through bioYear2, except for known mortalities in bioYear3
  mutate(exit.bio = ifelse((is.na(MortDate) | MortDate > bioYear2) & bioYear == "bioYear1", 365, # survived bioYear1, day 365
                           ifelse(is.na(MortDate) & bioYear == "bioYear2", lastDay,  # survived bioYear2, last day checked
                                  ifelse(MortDate > bioYear3 & bioYear == "bioYear2", lastDay, # survived bioYear2, putting as last day of bioYear2 # if analyzing into bioYear3, change lastday to 365
                                         ifelse(MortDate >= bioYear1 & MortDate < bioYear2, MortDate - bioYear1, # died in bioYear1, end date minus bioYear1
                                                ifelse((MortDate >= bioYear2 & MortDate < bioYear3) & bioYear == "bioYear2", MortDate - bioYear2, # died in bioYear2, end date minus bioYear2
                                                       ifelse(MortDate >= bioYear3 & bioYear == "bioYear3", MortDate - bioYear3, # died in bioYear3, end date minus bioYear3
                                                              "Check code!"))))))) %>%
         #Event indication (1=dead, 0=alive)
  mutate(event.bio = ifelse((is.na(MortDate) | MortDate > bioYear2) & bioYear == "bioYear1", 0, # survived bioYear1, = 0
                            ifelse((is.na(MortDate) | MortDate > bioYear3) & bioYear == "bioYear2", 0,  # survived bioYear2,= 0,
                                   ifelse(is.na(MortDate) & bioYear == "bioYear3", 0,  # survived bioYear3 = 0
                                          ifelse(MortDate >= bioYear1 & MortDate < bioYear2, 1, # died in bioYear1, = 1
                                                 ifelse(MortDate >= bioYear2 & MortDate < bioYear3 & bioYear == "bioYear2", 1,  # died in bioYear2, = 1
                                                        ifelse(MortDate >= bioYear3 & bioYear == "bioYear3", 1, # died in bioYear3, = 1
                                                               "Check code!")))))))
#Converting and formatting
adultelk <- adultelk %>%
  mutate(enter.bio = as.integer(enter.bio),
         exit.bio = as.integer(exit.bio),
         event.bio = as.integer(event.bio),
         bioYear.yrs = factor(bioYear.yrs))

##Capture numbers
length(unique(adultelk$AnimalID)) # total 
length(unique(adultelk$AnimalID[adultelk$Gender == "Female"])) # females
length(unique(adultelk$AnimalID[adultelk$Gender == "Male"])) # males

#Filter out bioYear3
adultelk <- filter(adultelk, bioYear != "bioYear3")

#####################################################################################
#Look at patterns of survival over all 3 years
#####################################################################################
attach(adultelk)
head(adultelk)

par(mfrow = c(1,2))
plot(Surv(enter.bio,exit.bio,event.bio), fn = "cum")
plot(Surv(enter.bio,exit.bio,event.bio), fn = "surv")
detach(adultelk)
##Look at female survival
adultelk1 <- subset(adultelk, Gender=="Female")
attach(adultelk1)
par(mfrow = c(2,2))
plot(Surv(enter.bio,exit.bio,event.bio), fn = "cum", main = "Cumulative hazard function females")
plot(Surv(enter.bio,exit.bio,event.bio), fn = "surv", main = "Survivor function females")
detach(adultelk1)
##Look at patterns of male survival
adultelk2 <- subset(adultelk, Gender=="Male")
attach(adultelk2)
plot(Surv(enter.bio,exit.bio,event.bio), fn = "cum", main = "Cumulative hazard function males")
plot(Surv(enter.bio,exit.bio,event.bio), fn = "surv", main = "Survivor function males")
detach(adultelk2)


#####################################################################################
#Model surival 
#####################################################################################
###Cow Annual KM estimates - ageGroup, rstudyArea, Year
annual.all<-survfit(Surv(enter.bio,exit.bio,event.bio)~1,conf.type="log-log",data=adultelk1)  #overall annual cow survival
summary(annual.all)

annual.Year<-survfit(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs,conf.type="log-log",data=adultelk1)
summary(annual.Year)

###Bull Annual KM estimates - ageGroup, rstudyArea, Year
annual.all<-survfit(Surv(enter.bio,exit.bio,event.bio)~1,conf.type="log-log",data=adultelk2)  #overall annual cow survival
summary(annual.all)

annual.Year<-survfit(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs,conf.type="log-log",data=adultelk2)
summary(annual.Year)



#####################################################################################
#Log-Rank tests for diff in year
#####################################################################################
###Annual KM estimates - Sex
summary(coxph(Surv(enter.bio,exit.bio,event.bio)~Gender, data=adultelk))

###Cow Annual KM estimates - Year
summary(coxph(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs, data=droplevels(adultelk1)))

###Bull Annual KM estimates - Year
summary(coxph(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs, data=droplevels(adultelk2))) 




############################################################################################################################
###KM plot graphics - survMisc##############################################################################################

## set the number of digits for plots to 2 with this function fmt
fmt <- function(){
    function(x) format(x,nsmall = 2,scientific = FALSE)
}

##load survMisc package (must be installed first) and run this code update 
require(survMisc)
autoplot.survfit <- function(object, ...,
xLab="Days since Jun 1",
yLab="Survival probability",
title="Marks show times with censoring",
titTextSize=15,
axisTitSize=15,
axisLabSize=10,
survLineSize=0.5,
type=c("single", "CI", "fill"),
palette=c("Dark2", "Set2", "Accent", "Paired",
"Pastel1", "Pastel2", "Set1", "Set3"),
jitter=c("none", "noEvents", "all"),
legend=TRUE,
legLabs=NULL,
legTitle="Strata",
legTextSize=10,
legSize=0.75,
alpha=0.05,
CIline=10,
fillLineSize=0.05,
pVal=FALSE,
sigP=1,
pX=0.1,
pY=0.1,
timeTicks=c("major", "minor"),
tabTitle="Number at risk by time",
tabTitTextSize=15,
tabLegTextSize=5,
nRiskSize=5){
stopifnot(inherits(object, "survfit"))
if(!is.null(legLabs) &! length(object$strata)==0) stopifnot(
length(legLabs)==length(object$strata))
### generate data to plot
### declare variables (for R CMD check)
### st1 is vector for strata identification
surv <- n.risk <- n.censor <- n.event <- upper <- lower <- NULL
.SD <- st1 <- stNames <- st <- s1 <- minT <- l <- maxT <- u <- NULL
### change names for strata to legLabs if required
if(is.null(legLabs)){
stNames <- names(object$strata)
} else {
stNames <- legLabs
}
### if only one strata (intercept only model)
if (is.null(object$strata)) {
if(is.null(legLabs)) {
st1 <- as.factor(rep(1, length(object$time)))
} else {
stopifnot(length(legLabs)==1)
st1 <- as.factor(rep(legLabs, length(object$time)))
}
} else {
### add vector for one strata according to number of rows of strata
st1 <- unlist(sapply(1:length(object$strata),
function (i) rep(stNames[i], object$strata[i]) ))
}
### create data.table with data from survfit
### add column for strata
### (using data.table here as avoids duplication when adding rows later)
### also rename strata as 'st' to avoid calling survival::function
dt1 <- data.table(time=object$time,
n.risk=object$n.risk,
n.event=object$n.event,
n.censor=object$n.censor,
surv=object$surv,
upper=object$upper,
lower=object$lower,
st=as.factor(st1))
### make two rows for each stratum
### for time=0 to time=time of first event
dt2 <- rbindlist(list(dt1[, .SD[1, ], by=st,allow.cartesian = TRUE],
dt1[, .SD[1, ], by=st,allow.cartesian = TRUE]))
### set n.event and n.censored to zero
dt2[, c("n.event", "n.censor") := list(0), by=st,allow.cartesian = TRUE]
### set surv, upper and lower to one
dt2[, c("surv", "upper", "lower") := list(1), by=st,allow.cartesian = TRUE]
### set first time to zero
dt2[seq(length(unique(dt2$st))), "time" := (0L),allow.cartesian = TRUE ]
### reorder to allow binding
setcolorder(dt2, names(dt1))
dt1 <- rbindlist(list(dt2, dt1))
### jitter
jitter <- match.arg(jitter)
### for groups with no events add random no.to survival (by strata)
if (jitter=="noEvents") {
### add column to indicate no. events by group
dt1[, s1 := sum(n.event), by=list(st),allow.cartesian = TRUE]
dt1[s1==0, surv := surv+(runif(1, 0.01, 0.05)), by=st,allow.cartesian = TRUE]
}
if(jitter=="all"){
### for groups with no events add random no.to survival (by strata)
dt1[, surv := surv+(runif(1, 0.01, 0.05)), by=st,allow.cartesian = TRUE]
}
###
dt1 <- dt1[order(st),allow.cartesian = TRUE]
### plot single lines only
g1 <- ggplot(data=dt1, aes(group=st, colour=st, fill=st)) +
geom_step(aes(x=time, y=surv), direction="hv", size=survLineSize)
###
type <- match.arg(type)
if (type=="CI"){
g1 <- g1 +
geom_step(aes(x=time, y=upper),
direction="hv", linetype=CIline, alpha=alpha) +
geom_step(aes(x=time, y=lower),
direction="hv", linetype=CIline, alpha=alpha)
}
if (type=="fill"){
### copy dt1 to work allow further work
dt2 <- dt1[, list(l=unique(lower),
u=unique(upper),
minT=as.numeric(min(time)),
time=as.numeric(time)
), by=list(surv, st),allow.cartesian = TRUE]
### make max. time column
dt2[, "maxT" := c(minT[2:length(minT)], NA), by=st,allow.cartesian = TRUE]
### merge columns
dt1 <- merge(dt1, dt2, by=c("time", "surv", "st"), all.y=TRUE,allow.cartesian=TRUE)
dt1 <- dt1[order(st),allow.cartesian = TRUE]
### add shading
g1 <- g1 + geom_rect(data=dt1, aes(x=NULL, y=NULL,
ymax=surv, ymin=l,
xmax=maxT, xmin=minT,
colour=st, group=st, fill=st),
alpha=alpha, size=fillLineSize) +
geom_rect(data=dt1, aes(x=NULL, y=NULL,
ymax=u, ymin=surv,
xmax=maxT, xmin=minT,
colour=st, group=st, fill=st),
alpha=alpha, size=fillLineSize)
}
### palette
### use palette Dark2 for prominent shades
### (suitable for colorblind)
### use palette Set2 for lighter shades as large fill area
palette <- match.arg(palette)
if(type=="fill"){
g1 <- g1 + scale_fill_brewer(type="qual", palette=palette,
guide=guide_legend(
keywidth=legSize,
keyheight=legSize))
}
g1 <- g1 + scale_colour_brewer(type="qual", palette=palette,
guide=guide_legend(
keywidth=legSize,
keyheight=legSize))
### scales
g1 <- g1 +
scale_y_continuous(yLab) +
ggtitle(title)
### times to show
timeTicks <- match.arg(timeTicks)
### use marks from existing plot
if(timeTicks=="major"){
times1 <- ggplot_build(g1)$panel$ranges[[1]]$x.major_source
} else {
times1 <- ggplot_build(g1)$panel$ranges[[1]]$x.minor_source
}
### x axis
g1 <- g1 +
scale_x_continuous(name=xLab,
breaks=times1)
### font sizes
g1 <- g1 +
theme(title = element_text(size=titTextSize),
legend.text=element_text(size=legTextSize),
legend.title=element_text(size=legTextSize),
axis.text = element_text(size = axisLabSize),
axis.title = element_text(size = axisTitSize)
)
### legend title
if(type=="fill"){
g1 <- g1 + labs(group=legTitle, colour=legTitle, fill=legTitle)
} else {
g1 <- g1 + labs(colour=legTitle)
}
### remove legend if required
if(!legend) g1 <- g1 + theme(legend.position = "none")
### p value for log-rank test (only if >=2 groups)
if(pVal & !is.null(object$strata)) {
sd1 <- survival::survdiff(eval(object$call$formula),
data=eval(object$call$data))
p1 <- stats::pchisq(sd1$chisq,
length(sd1$n) - 1,
lower.tail=FALSE)
p1txt <- ifelse(p1 < 0.0001,
"Log-rank test \n p < 0.0001",
paste("Log-rank test \n p =", signif(p1, sigP))
)
g1 <- g1 + annotate("text",
x = pX * max(dt1$time),
y = pY,
label = p1txt,
size = element_text(size=legTextSize))
}
### data for table
dt3 <- data.table(
time = summary(object, times = times1, extend = TRUE)$time,
n.risk = summary(object, times = times1, extend = TRUE)$n.risk
)
### if intercept-only model
if (is.null(object$strata)) {
dt3[, "st" := as.factor(rep(1, length(times1))),allow.cartesian = TRUE]
} else {
dt3[, "st" := summary(object, times=times1, extend=TRUE)$strata,allow.cartesian = TRUE]
}
### change names of strata to legend labels
if(!is.null(legLabs)) dt3[, "st" := factor(st, labels=legLabs),allow.cartesian = TRUE ]
### table
### reverse here to plot in same order as in main plot
g2 <- ggplot(data=dt3, aes(x=time, y=rev(st), shape=rev(st))) +
geom_point(size=0) +
geom_text(aes(label=n.risk), colour=1, size=nRiskSize) +
scale_x_continuous(name=xLab, limits=c(0, max(object$time)),
breaks=times1) +
### reverse here to plot in same order as in main plot
scale_y_discrete(name=legTitle, breaks=as.character(levels(dt3$st)),
labels=rev(levels(dt3$st))) +
ggtitle(tabTitle) +
theme(axis.text = element_text(size=axisLabSize),
axis.title = element_text(size=axisTitSize),
plot.title = element_text(size=tabTitTextSize),
legend.title = element_text(size=tabLegTextSize),
legend.text = element_text(size=tabLegTextSize)
) +
guides(shape = guide_legend(title=legTitle,
keywidht=legSize,
keyheight=legSize))
### remove legend
if(!legend) g2 <- g2 + theme(legend.position = "none")
res <- list("table"=g2,
"plot"=g1)
class(res) <- c("tableAndPlot", "list")
return(res)
}


#########################################################################################################
#### MAKE KAPLAN MEIER PLOTS IN GGPLOTS FOR ADULT FEMALES AND MALES                  
#########################################################################################################

#Set hunting seasons
arch.beg <- as.numeric(mdy("09/07/2014")-mdy("05/31/2014")) # Begin. Archery Season
rifl.beg <- as.numeric(mdy("10/26/2014")-mdy("05/31/2014")) # Begin. Rifle Season
rifl.end <- as.numeric(mdy("12/01/2014")-mdy("05/31/2014")) # End Rifle Season

##Make KM plot for cows by year ######################
reg.out1<-survfit(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs,conf.type="log-log",data=adultelk1)
summary(reg.out1)

p2 <- autoplot(reg.out1,
               title="",legTitle="",legLabs=c("2014-2015","2015-2016"),
               type="fill",survLineSize=2, palette="Set1",
               alpha=0.2)$plot
p4<-p2 + theme_bw() +
	 xlab("Days since Jun 1") + theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=12)) +
  	 ylab("Adult female survival probability") + theme(axis.text.y=element_text(size=12, color="black")) + theme(axis.title.y=element_text(size=12)) +
	geom_vline(xintercept=c(arch.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.end), linetype="dotted",size=1)+
	scale_y_continuous(breaks=seq(0,1,by=.2),limits=c(0,1)) +
  scale_x_continuous(breaks=c(seq(0,300,by=100),365)) +
 	theme(axis.title=element_text(size=12)) + #this sets the size of the axis lables
 	theme(legend.text=element_text(size=10), legend.title=element_text(size=10)) + # this sets the size for the legend title and text 
  annotate(geom="text", x=arch.beg-8, y=0, label="Archery", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.beg-8, y=0, label="Rifle", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.end-8, y=0, label="End Rifle", hjust=0, angle=90, size=3)
p4
##save file as PNG type and export file to working directory
ggsave(file="CowYearFIG.png",width=7,height=4)


########
##Make KM plot for bulls by year ######################
reg.out2<-survfit(Surv(enter.bio,exit.bio,event.bio)~bioYear.yrs,conf.type="log-log",data=adultelk2)
summary(reg.out2)

p2 <- autoplot(reg.out2,
               title="",legTitle="",legLabs=c("2014-2015","2015-2016"),
               type="fill",survLineSize=2, palette="Set1",
               alpha=0.2)$plot
p4<-p2 + theme_bw() +
	 xlab("Days since Jun 1") + theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=12)) +
  	 ylab("Adult male survival probability") + theme(axis.text.y=element_text(size=12, color="black")) + theme(axis.title.y=element_text(size=12)) +
	geom_vline(xintercept=c(arch.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.end), linetype="dotted",size=1)+
	scale_y_continuous(breaks=seq(0,1,by=.2),limits=c(0,1)) +
  scale_x_continuous(breaks=c(seq(0,300,by=100),365)) +
  theme(axis.title=element_text(size=12)) + #this sets the size of the axis lables
 	theme(legend.text=element_text(size=10), legend.title=element_text(size=10)) + # this sets the size for the legend title and text 
  annotate(geom="text", x=arch.beg-8, y=0, label="Archery", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.beg-8, y=0, label="Rifle", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.end-8, y=0, label="End Rifle", hjust=0, angle=90, size=3)
p4
##save file as PNG type and export file to working directory
ggsave(file="BullYearFIG.png",width=7,height=4)


########
##Make KM plot for adult elk by sex ######################
reg.out3<-survfit(Surv(enter.bio,exit.bio,event.bio)~Gender,conf.type="log-log",data=adultelk)
summary(reg.out3)

p2 <- autoplot(reg.out3,
               title="",legTitle="",legLabs=c("Female","Male"),
               type="fill",survLineSize=2, palette="Set1",
               alpha=0.2)$plot
p6<-p2 + theme_bw() +
  xlab("Days since Jun 1") + theme(axis.text.x=element_text(size=12, color="black"), axis.title.x=element_text(size=12)) +
  ylab("Adult elk survival probability") + theme(axis.text.y=element_text(size=12, color="black")) + theme(axis.title.y=element_text(size=12)) +
  geom_vline(xintercept=c(arch.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.beg), linetype="dotted",size=1)+geom_vline(xintercept=c(rifl.end), linetype="dotted",size=1)+
  scale_y_continuous(breaks=seq(0,1,by=.2),limits=c(0,1)) +
  scale_x_continuous(breaks=c(seq(0,300,by=100),365)) +
  theme(axis.title=element_text(size=12)) + #this sets the size of the axis lables
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10)) + # this sets the size for the legend title and text 
  annotate(geom="text", x=arch.beg-8, y=0, label="Archery", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.beg-8, y=0, label="Rifle", angle=90, hjust=0, size=3) + 
  annotate(geom="text", x=rifl.end-8, y=0, label="End Rifle", hjust=0, angle=90, size=3)
  p6
##save file as PNG type and export file to working directory
ggsave(file="ElksexFIG_new.png",width=7,height=4)





odbcCloseAll() # Closes connections to databases

