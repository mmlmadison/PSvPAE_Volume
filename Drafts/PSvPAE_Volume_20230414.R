###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, lmertest, misty)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/")


Alldata <- read.csv("PSvPAE_Longitudinal_20230414.csv")

#recoding variables into different data types, for example string into numeric
Alldata$ICV <- as.numeric(Alldata$ICV)
Alldata$ICV <- Alldata$ICV/10
#Alldata$Total_Gray <- Alldata$Total_Gray/10
Alldata$Age <- as.numeric(Alldata$Age)
Alldata$Subject <- as.factor(Alldata$Subject)
Alldata$PAE1_TD0 <- as.factor(Alldata$PAE1_TD0)
Alldata$Sex.F1_M0. <- as.factor(Alldata$Sex.F1_M0.)

#Subsetting Alldata into PAE and PS for Sex analyses
PAEdata <- subset(Alldata, PAE1_TD0=="1")
PSdata <- subset(Alldata, PAE1_TD0=="0")

M_data <-subset(PAEdata,Sex.F1_M0.=='0')
F_data <-subset(PAEdata,Sex.F1_M0.=='1')



######################## PS v PAE mixed effects analysis and plots

regions <-as.data.frame(PSdata[,c(126:242)])
regionsPAE <- as.data.frame(PAEdata[,c(126:242)])
regionsPAEvPS <- as.data.frame(Alldata[,c(126:242)])
M_data_regions <- as.data.frame(M_data[,c(126:242)])
F_data_regions <- as.data.frame(F_data[,c(126:242)])

dim(regions)
dim(regionsPAEvPS)
dim(M_data_regions)
dim(F_data_regions)
plot_list <- list()
plot_list_Sex <- list()

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  ROIPAE <- regionsPAE[,i]
  ROIPAEvPS <- regionsPAEvPS[,i]
  ROIM <- M_data_regions[,i]
  ROIF <- F_data_regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  Fplot <-0
  PAEplot <- 0
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  ## Model selection with Preschool sample
  #null model
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  print(colnames(regions)[i])
  
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    
    #Run a linear model to test main and int. effects of Group (PSvPAE in this case)
    lmer_volume_lin_PAE <- lmer(ROIPAEvPS ~ Age*PAE1_TD0 + (1|Subject), REML = FALSE, data = Alldata)
    
    #Run a linear model to test main and int. effects of Sex within the PAE group
    lmer_volume_lin_PAESex <- lmer(ROIPAE ~ Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    
    #Extract values from the coefficients table to use in plotting and change estimates
    CoefDataLinPAE <- as.data.frame(fixef(lmer_volume_lin_PAE))
    Interc <- as.numeric(CoefDataLinPAE[1,1])
    AgeCoef <- as.numeric(CoefDataLinPAE[2,1])
    GroupCoef <- as.numeric(CoefDataLinPAE[3,1])
    AgebyGroupCoef <- as.numeric(CoefDataLinPAE[4,1])
    
    CoefDataLinPAESex <- as.data.frame(fixef(lmer_volume_lin_PAESex))
    Interc_Sex <- as.numeric(CoefDataLinPAESex[1,1])
    AgeCoef_Sex <- as.numeric(CoefDataLinPAESex[2,1])
    GroupCoef_Sex <- as.numeric(CoefDataLinPAESex[3,1])
    AgebyGroupCoef_Sex <- as.numeric(CoefDataLinPAESex[4,1])
    
    PAEslope <- AgeCoef+AgebyGroupCoef*1
    PAEinterc <- Interc+GroupCoef*1
    
    TDslope <- AgeCoef+AgebyGroupCoef*0
    TDinterc <- Interc+GroupCoef*0
    
    Femslope <- AgeCoef_Sex+AgebyGroupCoef_Sex*1
    Feminterc <- Interc_Sex+GroupCoef_Sex*1
    
    Maleslope <- AgeCoef_Sex+AgebyGroupCoef_Sex*0
    Maleinterc <- Interc_Sex+GroupCoef_Sex*0
    
    #Defining the functions of the two mean trajectories
    ROIfnlin_PAE <- function(x) {
      as.numeric(PAEslope*x + PAEinterc) # outcome is y
    }
    
    ROIfnlin_TD <- function(x) {
      as.numeric(TDslope*x + TDinterc) # outcome is y
    }
    
    ROIfnlin_Fem <- function(x) {
      as.numeric(Femslope*x + Feminterc) # outcome is y
    }
    
    ROIfnlin_Male <- function(x) {
      as.numeric(Maleslope*x + Maleinterc) # outcome is y
    }
    
    
    lmerpredlin_PAE
    lmerpredlin_PAESex
    #Plot the group difference for PAEvPS with the current ROI
    plotROI <- ggplot(Alldata, aes(x=Age, y=ROIPAEvPS, group=Subject)) 
    TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=ROI),color='green',size=.6) + 
      geom_line(data=PSdata, aes(x=Age,y=ROI), color="green", size=.2) + geom_function(fun = ROIfnlin_TD, color="green", size=2)
    PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=ROIPAE),color='purple',size=.6) + 
      geom_line(data=PAEdata, aes(x=Age,y=ROIPAE), color="purple", size=.2) + geom_function(fun = ROIfnlin_PAE, color="purple", size=2) +
      theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    #Plot the group difference for Sex in the PAE group with the current ROI
    plotROI_Sex <- ggplot(PAEdata, aes(x=Age, y=ROIPAE, group=Subject)) 
    Mplot <- plotROI_Sex + geom_point(data=M_data, aes(x=Age,y=ROIM),color='blue',size=.6) + 
      geom_line(data=M_data, aes(x=Age,y=ROIM), color="blue", size=.2) + geom_function(fun = ROIfnlin_Male, color="blue", size=2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='red',size=.6) + 
      geom_line(data=F_data, aes(x=Age,y=ROIF), color="red", size=.2) + geom_function(fun = ROIfnlin_Fem, color="red", size=2) +
      theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- AIC_null
    
    Fplot <- ggplot()
    PAEplot <-ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    #Run a 2nd order linear model to test main and int. effects of Group (PSvPAE in this case)
    lmer_volume_quad_PAE <- lmer(ROIPAEvPS ~ I(Age^2)*PAE1_TD0 + Age*PAE1_TD0 + (1|Subject), REML = FALSE, data = Alldata)
    
    #Run a 2nd order linear model to test main and int. effects of Sex within the PAE group
    lmer_volume_quad_PAESex <- lmer(ROIPAE ~ I(Age^2)*Sex.F1_M0. + Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    
    #Extract values from the coefficients table to use in plotting and change estimates
    CoefDataQuadPAE <- as.data.frame(fixef(lmer_volume_quad_PAE))
    Interc <- as.numeric(CoefDataQuadPAE[1,1])
    Age2Coef <- as.numeric(CoefDataQuadPAE[2,1])
    GroupCoef <- as.numeric(CoefDataQuadPAE[3,1])
    AgeCoef <- as.numeric(CoefDataQuadPAE[4,1])
    Age2byGroupCoef <- as.numeric(CoefDataQuadPAE[5,1])
    AgebyGroupCoef <- as.numeric(CoefDataQuadPAE[6,1])
    
    CoefDataQuadPAESex <- as.data.frame(fixef(lmer_volume_quad_PAESex))
    Interc_Sex <- as.numeric(CoefDataQuadPAESex[1,1])
    Age2Coef_Sex <- as.numeric(CoefDataQuadPAESex[2,1])
    GroupCoef_Sex <- as.numeric(CoefDataQuadPAESex[3,1])
    AgeCoef_Sex <- as.numeric(CoefDataQuadPAESex[4,1])
    Age2byGroupCoef_Sex <- as.numeric(CoefDataQuadPAESex[5,1])
    AgebyGroupCoef_Sex <- as.numeric(CoefDataQuadPAESex[6,1])
    
    PAEslope_lin <- AgeCoef+AgebyGroupCoef*1
    PAEslope_quad <- Age2Coef+Age2byGroupCoef*1
    PAEinterc <- Interc+GroupCoef*1
    
    TDslope_lin <- AgeCoef+AgebyGroupCoef*0
    TDslope_quad <- Age2Coef+Age2byGroupCoef*0
    TDinterc <- Interc+GroupCoef*0
    
    Femslope_lin <- AgeCoef_Sex+AgebyGroupCoef_Sex*1
    Femslope_quad <- Age2Coef_Sex+Age2byGroupCoef_Sex*1
    Feminterc <- Interc_Sex+GroupCoef_Sex*1
    
    Maleslope_lin <- AgeCoef_Sex+AgebyGroupCoef_Sex*0
    Maleslope_quad <- Age2Coef_Sex+Age2byGroupCoef_Sex*0
    Maleinterc <- Interc_Sex+GroupCoef_Sex*0
    
    #Defining the functions of the two mean trajectories
    ROIfnquad_PAE <- function(x) {
      as.numeric(PAEslope_quad*x^2 + PAEslope_lin*x + PAEinterc) # outcome is y
    }
    
    ROIfnquad_TD <- function(x) {
      as.numeric(TDslope_quad*x^2 + TDslope_lin*x + TDinterc) # outcome is y
    }
    
    ROIfnquad_Fem <- function(x) {
      as.numeric(Femslope_quad*x^2 + Femslope_lin*x + Feminterc) # outcome is y
    }
    
    ROIfnquad_Male <- function(x) {
      as.numeric(Maleslope_quad*x^2 + Maleslope_lin*x + Maleinterc) # outcome is y
    }
    
    #Plot the group difference for PAEvPS with the current ROI
    plotROI <- ggplot(Alldata, aes(x=Age, y=ROIPAEvPS, group=Subject)) 
    TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=ROI),color='green',size=.6) + 
      geom_line(data=PSdata, aes(x=Age,y=ROI), color="green", size=.2) + geom_function(fun = ROIfnquad_TD, color="green", size=2)
    PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=ROIPAE),color='purple',size=.6) + 
      geom_line(data=PAEdata, aes(x=Age,y=ROIPAE), color="purple", size=.2) + geom_function(fun = ROIfnquad_PAE, color="purple", size=2) +
      theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    #Plot the group difference for Sex in the PAE group with the current ROI
    plotROI_Sex <- ggplot(PAEdata, aes(x=Age, y=ROIPAE, group=Subject)) 
    Mplot <- plotROI_Sex + geom_point(data=M_data, aes(x=Age,y=ROIM),color='blue',size=.6) + 
      geom_line(data=M_data, aes(x=Age,y=ROIM), color="blue", size=.2) + geom_function(fun = ROIfnquad_Male, color="blue", size=2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='red',size=.6) + 
      geom_line(data=F_data, aes(x=Age,y=ROIF), color="red", size=.2) + geom_function(fun = ROIfnquad_Fem, color="red", size=2) +
      theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  

  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  plot_list[[i]] <- ggplotGrob(PAEplot)
  plot_list_Sex[[i]] <- ggplotGrob(Fplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_20230417.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

bigPlot_Sex <- marrangeGrob(grobs=plot_list_Sex, ncol=4, nrow=2)
ggsave("BigPlot_PAE_Sex_20230417.pdf", bigPlot_Sex, width=50, height=25, limitsize = FALSE)



