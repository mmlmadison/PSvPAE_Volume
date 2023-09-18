###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, lmertest, misty, zoo)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/")


Alldata <- read.csv("PSvPAE_Longitudinal_20230420.csv")

#recoding variables into different data types, for example string into numeric
Alldata$ICV <- as.numeric(Alldata$ICV)
#Alldata$ICV <- Alldata$ICV/10
#Alldata$Total_Gray <- Alldata$Total_Gray/10
Alldata$Age <- as.numeric(Alldata$Age)
Alldata$Subject <- as.factor(Alldata$Subject)
Alldata$PAE1_TD0 <- as.factor(Alldata$PAE1_TD0)
Alldata$Sex.F1_M0. <- as.factor(Alldata$Sex.F1_M0.)


lmer_tableinit_PAE <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef.PAE <- as.data.frame(coef(summary(lmer_tableinit_PAE)))
df.master.coef.PAE["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.PAE <- as.data.frame(effectsize::eta_squared(lmer_tableinit_PAE))
df.master.etasq.PAE["Region"] <-c("Right_Inf_Lat_Vent")

lmer_tableinit_Sex <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef.Sex <- as.data.frame(coef(summary(lmer_tableinit_Sex)))
df.master.coef.Sex["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.Sex <- as.data.frame(effectsize::eta_squared(lmer_tableinit_Sex))
df.master.etasq.Sex["Region"] <-c("Right_Inf_Lat_Vent")

######################## PS v PAE mixed effects analysis and plots

regions <-as.data.frame(PSdata[,c(36:36)])
plot_list <- list()
plot_list_Sex <- list()

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  Fplot <-0
  PAEplot <- 0
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  ## Model selection with Preschool sample
  #null model
  lmer_volume_null <- lmer(Right_Calc_calcarine_cortex ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(Right_Calc_calcarine_cortex ~ Age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(Right_Calc_calcarine_cortex ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  print(colnames(regions)[i])
  
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin

  }else {
    bestAIC <- AIC_null
    finalPAE <- lmer_volume_null
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    
    #Run a 2nd order linear model to test main and int. effects of Group (PSvPAE in this case)
    regionsPAEvPS <- as.data.frame(Alldata[,c(36:36)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_quad_PAE <- lmer(ROIPAEvPS ~ I(Age^2)*PAE1_TD0 + Age*PAE1_TD0 + (1|Subject), REML = FALSE, data = Alldata)
    Alldata$lmerpredquad_PAE <- predict(lmer_volume_quad_PAE)
    etasqPAE <- effectsize::eta_squared(lmer_volume_quad_PAE)
    finalPAE <- lmer_volume_quad_PAE
    
    PAEdata <- subset(Alldata, PAE1_TD0=="1")
    PSdata <- subset(Alldata, PAE1_TD0=="0")
    
    #Run a 2nd order linear model to test main and int. effects of Sex within the PAE group
    regionsPAE <- as.data.frame(PAEdata[,c(36:36)])
    ROIPAE <- regionsPAE[,i]
    lmer_volume_quad_PAESex <- lmer(ROIPAE ~ I(Age^2)*Sex.F1_M0. + Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    PAEdata$lmerpredquad_PAESex <- predict(lmer_volume_quad_PAESex)
    etasqSex <- effectsize::eta_squared(lmer_volume_quad_PAESex)
    finalPAESex <- lmer_volume_quad_PAESex
    
    #Extract values from the coefficients table to use in plotting and change estimates
    CoefDataQuadPAE <- as.data.frame(fixef(lmer_volume_quad_PAE))
    Interc <- as.numeric(CoefDataQuadPAE[1,1])
    Age2Coef <- as.numeric(CoefDataQuadPAE[2,1])
    GroupCoef <- as.numeric(CoefDataQuadPAE[3,1])
    AgeCoef <- as.numeric(CoefDataQuadPAE[4,1])
    Age2byGroupCoef <- as.numeric(CoefDataQuadPAE[5,1])
    AgebyGroupCoef <- as.numeric(CoefDataQuadPAE[6,1])
    
    PAEslope_lin <- AgeCoef+AgebyGroupCoef*1
    PAEslope_quad <- Age2Coef+Age2byGroupCoef*1
    PAEinterc <- Interc+GroupCoef*1
    
    TDslope_lin <- AgeCoef+AgebyGroupCoef*0
    TDslope_quad <- Age2Coef+Age2byGroupCoef*0
    TDinterc <- Interc+GroupCoef*0
  
    
    #Defining the functions of the two mean trajectories
    ROIfnquad_PAE <- function(x) {
      as.numeric(PAEslope_quad*x^2 + PAEslope_lin*x + PAEinterc) # outcome is y
    }
    
    ROIfnquad_TD <- function(x) {
      as.numeric(TDslope_quad*x^2 + TDslope_lin*x + TDinterc) # outcome is y
    }
    
    x <- PSdata[,5]
    y <- ROI
    out <- paste(colName,"_LOESS.png", sep="")
    
    LOESS_plot_init <- ggplot(PSdata, aes(x=Age, y=ROI)) 
    
    LOESS_plot_PAE <- LOESS_plot_init + geom_point(data=PAEdata, aes(x=Age, y=ROIPAE, group=Subject), color="purple2", size=.6) + geom_smooth(data=PAEdata, aes(x=Age, y=ROIPAE, group=Subject), method=lm, linewidth=.4, se=FALSE, color="purple3")
    LOESS_plot_Finish <- LOESS_plot_PAE + theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + geom_function(fun = ROIfnquad_PAE, color="purple4", size=2) +
      coord_cartesian(xlim=c(2, 8)) + coord_cartesian(ylim=c(250, 7000)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    LOESS_plot_Finish
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  
  df.etasq.PAE <- as.data.frame(etasqPAE)
  df.etasq.PAE['Region'] <- c(colnames(regions)[i])
  df.master.etasq.PAE <- rbind(df.master.etasq.PAE, df.etasq.PAE)
  
  df.etasq.Sex <- as.data.frame(etasqSex)
  df.etasq.Sex['Region'] <- c(colnames(regions)[i])
  df.master.etasq.Sex <- rbind(df.master.etasq.Sex, df.etasq.Sex)
  
  df.coef.PAE <- as.data.frame(coef(summary(finalPAE)))
  df.coef.PAE['Region'] <- c(colnames(regions)[i])
  df.master.coef.PAE <- rbind(df.master.coef.PAE, df.coef.PAE)
  
  df.coef.Sex <- as.data.frame(coef(summary(finalPAESex)))
  df.coef.Sex['Region'] <- c(colnames(regions)[i])
  df.master.coef.Sex <- rbind(df.master.coef.Sex, df.coef.Sex)
  
  #write.csv(df.master.coef.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_AbsVolume_20230427.csv")
  #write.csv(df.master.etasq.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/EtaSqTable_PSvPAE_AbsVolume_20230427.csv")
  #write.csv(df.master.coef.Sex, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PAESex_AbsVolume_20230427.csv")
  #write.csv(df.master.etasq.Sex, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/EtaSqTable_PAESex_AbsVolume_20230427.csv")

  ggsave(paste0(colName,"_PSvPAE_LOESS_Owerko_20230527.png"), plot=LOESS_plot_Finish, width=11, height=8, limitsize = FALSE)
  print(colName)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}

