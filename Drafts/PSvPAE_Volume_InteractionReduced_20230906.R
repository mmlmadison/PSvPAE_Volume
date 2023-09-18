###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmerTest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, misty)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/")

Alldata <- read.csv("PSvPAE_Segmentations_LongitudinalRegistration_20230822.csv")

#recoding variables into different data types, for example string into numeric
Alldata$ICV <- as.numeric(Alldata$ICV)
#Alldata$ICV <- Alldata$ICV/10
#Alldata$Total_Gray <- Alldata$Total_Gray/10
Alldata$Age <- as.numeric(Alldata$Age)
Alldata$Subject <- as.factor(Alldata$Subject)
Alldata$PAE <- as.factor(Alldata$PAE)
Alldata$Sex.F1_M0. <- as.factor(Alldata$Sex.F1_M0.)

Alldata <- subset(Alldata, Age < 8.5)
PSdata <- subset(Alldata, PAE=="0")
PAEdata <- subset(Alldata, PAE=="1")

lmer_tableinit_PAE <- lmer(Right.Accumbens.Area ~ Age + (1|Subject), REML = FALSE, data = Alldata)
df.master.coef.PAE <- as.data.frame(coef(summary(lmer_tableinit_PAE)))
df.master.coef.PAE["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.PAE <- as.data.frame(effectsize::eta_squared(lmer_tableinit_PAE))
df.master.etasq.PAE["Region"] <-c("Right_Inf_Lat_Vent")

lmer_tableinit_Sex <- lmer(Right.Accumbens.Area ~ Age + (1|Subject), REML = FALSE, data = Alldata)
df.master.coef.Sex <- as.data.frame(coef(summary(lmer_tableinit_Sex)))
df.master.coef.Sex["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.Sex <- as.data.frame(effectsize::eta_squared(lmer_tableinit_Sex))
df.master.etasq.Sex["Region"] <-c("Right_Inf_Lat_Vent")

######################## PS v PAE mixed effects analysis and plots

regions <-as.data.frame(PSdata[,c(10:126)])

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
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
    
    #Run a linear model to test main and int. effects of Group (PSvPAE in this case). Save predicted values.
    regionsPAEvPS <- as.data.frame(Alldata[,c(10:126)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_lin_PAE <- lmer(ROIPAEvPS ~ Age + PAE + (1|Subject), REML = FALSE, data = Alldata)
    #lmer_volume_lin_PAERAND <- lmer(ROIPAEvPS ~ Age*PAE1_TD0 + (1 + Age|Subject), REML = FALSE, data = Alldata)
    #Alldata$lmerpredlin_PAE <- predict(lmer_volume_lin_PAERAND)
    etasqPAE <- effectsize::eta_squared(lmer_volume_lin_PAE)
    finalPAE <- lmer_volume_lin_PAE
    
    #Subsetting Alldata into PAE and PS for Sex analysis and plotting
    PAEdata <- subset(Alldata, PAE=="1")
    #PSdata <- subset(Alldata, PAE=="0")
    
    #Run a linear model to test main and int. effects of Sex within the PAE group. Save predicted values.
    regionsPAE <- as.data.frame(PAEdata[,c(10:126)])
    ROIPAE <- regionsPAE[,i]
    #lmer_volume_lin_PAESex <- lmer(ROIPAE ~ Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    #lmer_volume_lin_PAESexRAND <- lmer(ROIPAE ~ Age + Sex.F1_M0. + (1 + Age|Subject), REML = FALSE, data = PAEdata)
    #PAEdata$lmerpredlin_PAESex <- predict(lmer_volume_lin_PAESexRAND)
    etasqSex <- effectsize::eta_squared(lmer_volume_lin_PAESex)
    finalPAESex <- lmer_volume_lin_PAESex
    
    #Extract values from the coefficients table to use in plotting and change estimates
    
  }else {
    bestAIC <- AIC_null
    finalPAE <- lmer_volume_null
    finalPAESex <- lmer_volume_null
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    
    #Run a 2nd order linear model to test main and int. effects of Group (PSvPAE in this case)
    regionsPAEvPS <- as.data.frame(Alldata[,c(10:126)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_quad_PAE <- lmer(ROIPAEvPS ~ I(Age^2) + Age + PAE + (1|Subject), REML = FALSE, data = Alldata)
    Alldata$lmerpredquad_PAE <- predict(lmer_volume_quad_PAE)
    etasqPAE <- effectsize::eta_squared(lmer_volume_quad_PAE)
    finalPAE <- lmer_volume_quad_PAE
    
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  
  #df.etasq.PAE <- as.data.frame(etasqPAE)
  #df.etasq.PAE['Region'] <- c(colnames(regions)[i])
  #df.master.etasq.PAE <- rbind(df.master.etasq.PAE, df.etasq.PAE)
  
  #df.etasq.Sex <- as.data.frame(etasqSex)
  #df.etasq.Sex['Region'] <- c(colnames(regions)[i])
  #df.master.etasq.Sex <- rbind(df.master.etasq.Sex, df.etasq.Sex)
  
  df.coef.PAE <- as.data.frame(coef(summary(finalPAE)))
  df.coef.PAE['Region'] <- c(colnames(regions)[i])
  df.master.coef.PAE <- rbind(df.master.coef.PAE, df.coef.PAE)
  
  #df.coef.Sex <- as.data.frame(coef(summary(finalPAESex)))
  #df.coef.Sex['Region'] <- c(colnames(regions)[i])
  #df.master.coef.Sex <- rbind(df.master.coef.Sex, df.coef.Sex)
  
  write.csv(df.master.coef.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_InteractionReduced_AbsVolume_20230906.csv")

  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}

#FDR correction
pvaldf <- read.csv("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_InteractionReduced_AbsVolume_FDR_20230906.csv")

qval <- p.adjust(p=pvaldf$pval, method="fdr")
pvaldf$qval <- qval
View(pvaldf)
write.csv(pvaldf, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_InteractionReduced_AbsVolume_FDR_qval_20230906.csv")



########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
####################################### ICV as covariate

lmer_tableinit_PAE <- lmer(Right.Accumbens.Area ~ Age + (1|Subject), REML = FALSE, data = Alldata)
df.master.coef.PAE <- as.data.frame(coef(summary(lmer_tableinit_PAE)))
df.master.coef.PAE["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.PAE <- as.data.frame(effectsize::eta_squared(lmer_tableinit_PAE))
df.master.etasq.PAE["Region"] <-c("Right_Inf_Lat_Vent")

lmer_tableinit_Sex <- lmer(Right.Accumbens.Area ~ Age + (1|Subject), REML = FALSE, data = Alldata)
df.master.coef.Sex <- as.data.frame(coef(summary(lmer_tableinit_Sex)))
df.master.coef.Sex["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.Sex <- as.data.frame(effectsize::eta_squared(lmer_tableinit_Sex))
df.master.etasq.Sex["Region"] <-c("Right_Inf_Lat_Vent")

######################## PS v PAE mixed effects analysis and plots

regions <-as.data.frame(PSdata[,c(10:126)])
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
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + ICV + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + ICV + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  print(colnames(regions)[i])
  
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    
    #Run a linear model to test main and int. effects of Group (PSvPAE in this case). Save predicted values.
    regionsPAEvPS <- as.data.frame(Alldata[,c(10:126)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_lin_PAE <- lmer(ROIPAEvPS ~ Age + PAE + ICV + (1|Subject), REML = FALSE, data = Alldata)
    #lmer_volume_lin_PAERAND <- lmer(ROIPAEvPS ~ Age*PAE1_TD0 + (1 + Age|Subject), REML = FALSE, data = Alldata)
    #Alldata$lmerpredlin_PAE <- predict(lmer_volume_lin_PAERAND)
    etasqPAE <- effectsize::eta_squared(lmer_volume_lin_PAE)
    finalPAE <- lmer_volume_lin_PAE
    
    #Subsetting Alldata into PAE and PS for Sex analysis and plotting
    PAEdata <- subset(Alldata, PAE=="1")
    PSdata <- subset(Alldata, PAE=="0")
    
  }else {
    bestAIC <- AIC_null
    finalPAE <- lmer_volume_null
    finalPAESex <- lmer_volume_null
    Fplot <- ggplot()
    PAEplot <-ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    
    #Run a 2nd order linear model to test main and int. effects of Group (PSvPAE in this case)
    regionsPAEvPS <- as.data.frame(Alldata[,c(10:126)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_quad_PAE <- lmer(ROIPAEvPS ~ I(Age^2) + Age + PAE + ICV + (1|Subject), REML = FALSE, data = Alldata)
    Alldata$lmerpredquad_PAE <- predict(lmer_volume_quad_PAE)
    etasqPAE <- effectsize::eta_squared(lmer_volume_quad_PAE)
    finalPAE <- lmer_volume_quad_PAE
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  df.coef.PAE <- as.data.frame(coef(summary(finalPAE)))
  df.coef.PAE['Region'] <- c(colnames(regions)[i])
  df.master.coef.PAE <- rbind(df.master.coef.PAE, df.coef.PAE)
  
  write.csv(df.master.coef.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_ICVCoV_20230906.csv")
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')

}

#FDR correction
pvaldf <- read.csv("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_ICVCoV_FDR_20230907.csv")

qval <- p.adjust(p=pvaldf$pval, method="fdr")
pvaldf$qval <- qval
View(pvaldf)


