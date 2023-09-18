##### Flux 2023 Figure plots for selected significant main effects

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmerTest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, misty)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/")

Alldata <- read.csv("PSvPAE_Segmentations_LongitudinalRegistration_20230822.csv")

#recoding variables into different data types, for example string into numeric
Alldata <- subset(Alldata, Age < 8.5)
Alldata$ICV <- as.numeric(Alldata$ICV)
Alldata$Age <- as.numeric(Alldata$Age)
Alldata$Subject <- as.factor(Alldata$Subject)
Alldata$PAE <- as.factor(Alldata$PAE)

## Left Posterior Cingulate
lmer_volume_lin_PAE <- lmer(Left.PCgG.posterior.cingulate.gyrus ~ I(Age^2) + Age + PAE + (1|Subject), REML = FALSE, data = Alldata)

#Subsetting Alldata into PAE and PS for Sex analysis and plotting
PAEdata <- subset(Alldata, PAE=="1")
PSdata <- subset(Alldata, PAE=="0")

#Extract values from the coefficients table to use in plotting and change estimates
CoefDataLinPAE <- as.data.frame(fixef(lmer_volume_lin_PAE))
Interc <- as.numeric(CoefDataLinPAE[1,1])
Age2Coef <- as.numeric(CoefDataLinPAE[2,1])
AgeCoef <- as.numeric(CoefDataLinPAE[3,1])
GroupCoef <- as.numeric(CoefDataLinPAE[3,1])


PAEinterc <- Interc+GroupCoef*1

TDinterc <- Interc+GroupCoef*0

#Defining the functions of the two mean trajectories
ROIfnlin_PAE <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + PAEinterc) # outcome is y
}

ROIfnlin_TD <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + TDinterc) # outcome is y
}


#Plot the group difference for PAEvPS with the current ROI
plotROI <- ggplot(Alldata, aes(x=Age, y=Left.PCgG.posterior.cingulate.gyrus, group=Subject)) 
TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=Left.PCgG.posterior.cingulate.gyrus),color='orange',size=.6) + 
  geom_smooth(data=PSdata, aes(x=Age, y=Left.PCgG.posterior.cingulate.gyrus), method=lm, se=FALSE, color="orange", size=.2)
PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=Left.PCgG.posterior.cingulate.gyrus),color='turquoise2',size=.6) + 
  geom_smooth(data=PAEdata, aes(x=Age, y=Left.PCgG.posterior.cingulate.gyrus), method=lm, se=FALSE, color="turquoise2", size=.2) + 
  geom_function(fun = ROIfnlin_PAE, color="turquoise3", size=2) + geom_function(fun = ROIfnlin_TD, color="orange3", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Left Posterior Cingulate ("~mm^3*")")


## Right Caudate
lmer_volume_lin_PAE <- lmer(Right.Caudate ~ I(Age^2) + Age + PAE + (1|Subject), REML = FALSE, data = Alldata)

#Subsetting Alldata into PAE and PS for Sex analysis and plotting
PAEdata <- subset(Alldata, PAE=="1")
PSdata <- subset(Alldata, PAE=="0")

#Extract values from the coefficients table to use in plotting and change estimates
CoefDataLinPAE <- as.data.frame(fixef(lmer_volume_lin_PAE))
Interc <- as.numeric(CoefDataLinPAE[1,1])
Age2Coef <- as.numeric(CoefDataLinPAE[2,1])
AgeCoef <- as.numeric(CoefDataLinPAE[3,1])
GroupCoef <- as.numeric(CoefDataLinPAE[4,1])


PAEinterc <- Interc+GroupCoef*1

TDinterc <- Interc+GroupCoef*0

#Defining the functions of the two mean trajectories
ROIfnlin_PAE <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + PAEinterc) # outcome is y
}

ROIfnlin_TD <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + TDinterc) # outcome is y
}


#Plot the group difference for PAEvPS with the current ROI
plotROI <- ggplot(Alldata, aes(x=Age, y=Right.Caudate, group=Subject)) 
TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=Right.Caudate),color='orange',size=.6) + 
  geom_smooth(data=PSdata, aes(x=Age, y=Right.Caudate), method=lm, se=FALSE, color="orange", size=.2)
PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=Right.Caudate),color='turquoise2',size=.6) + 
  geom_smooth(data=PAEdata, aes(x=Age, y=Right.Caudate), method=lm, se=FALSE, color="turquoise2", size=.2) + 
  geom_function(fun = ROIfnlin_PAE, color="turquoise3", size=2) + geom_function(fun = ROIfnlin_TD, color="orange3", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Caudate ("~mm^3*")")


## Right Middle Frontal Gyrus
lmer_volume_lin_PAE <- lmer(Right.MFG.middle.frontal.gyrus ~ I(Age^2) + Age + PAE + (1|Subject), REML = FALSE, data = Alldata)

#Subsetting Alldata into PAE and PS for Sex analysis and plotting
PAEdata <- subset(Alldata, PAE=="1")
PSdata <- subset(Alldata, PAE=="0")

#Extract values from the coefficients table to use in plotting and change estimates
CoefDataLinPAE <- as.data.frame(fixef(lmer_volume_lin_PAE))
Interc <- as.numeric(CoefDataLinPAE[1,1])
Age2Coef <- as.numeric(CoefDataLinPAE[2,1])
AgeCoef <- as.numeric(CoefDataLinPAE[3,1])
GroupCoef <- as.numeric(CoefDataLinPAE[4,1])


PAEinterc <- Interc+GroupCoef*1

TDinterc <- Interc+GroupCoef*0

#Defining the functions of the two mean trajectories
ROIfnlin_PAE <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + PAEinterc) # outcome is y
}

ROIfnlin_TD <- function(x) {
  as.numeric(Age2Coef*I(x^2) + AgeCoef*x + TDinterc) # outcome is y
}


#Plot the group difference for PAEvPS with the current ROI
plotROI <- ggplot(Alldata, aes(x=Age, y=Right.MFG.middle.frontal.gyrus, group=Subject)) 
TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=Right.MFG.middle.frontal.gyrus),color='orange',size=.6) + 
  geom_smooth(data=PSdata, aes(x=Age, y=Right.MFG.middle.frontal.gyrus), method=lm, se=FALSE, color="orange", size=.2)
PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=Right.MFG.middle.frontal.gyrus),color='turquoise2',size=.6) + 
  geom_smooth(data=PAEdata, aes(x=Age, y=Right.MFG.middle.frontal.gyrus), method=lm, se=FALSE, color="turquoise2", size=.2) + 
  geom_function(fun = ROIfnlin_PAE, color="turquoise3", size=2) + geom_function(fun = ROIfnlin_TD, color="orange3", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Middle Frontal Gyrus ("~mm^3*")")
