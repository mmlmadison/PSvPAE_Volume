###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, lmertest, misty)

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

regions <-as.data.frame(PSdata[,c(9:125)])
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
  lmer_volume_lin <- lmer(ROI ~ Age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  print(colnames(regions)[i])
  
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    
    #Run a linear model to test main and int. effects of Group (PSvPAE in this case). Save predicted values.
    regionsPAEvPS <- as.data.frame(Alldata[,c(9:125)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_lin_PAE <- lmer(ROIPAEvPS ~ Age*PAE1_TD0 + (1|Subject), REML = FALSE, data = Alldata)
    lmer_volume_lin_PAERAND <- lmer(ROIPAEvPS ~ Age*PAE1_TD0 + (1 + Age|Subject), REML = FALSE, data = Alldata)
    Alldata$lmerpredlin_PAE <- predict(lmer_volume_lin_PAERAND)
    etasqPAE <- effectsize::eta_squared(lmer_volume_lin_PAE)
    finalPAE <- lmer_volume_lin_PAE
    
    #Subsetting Alldata into PAE and PS for Sex analysis and plotting
    PAEdata <- subset(Alldata, PAE1_TD0=="1")
    PSdata <- subset(Alldata, PAE1_TD0=="0")
    
    #Run a linear model to test main and int. effects of Sex within the PAE group. Save predicted values.
    regionsPAE <- as.data.frame(PAEdata[,c(9:125)])
    ROIPAE <- regionsPAE[,i]
    lmer_volume_lin_PAESex <- lmer(ROIPAE ~ Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    lmer_volume_lin_PAESexRAND <- lmer(ROIPAE ~ Age + Sex.F1_M0. + (1 + Age|Subject), REML = FALSE, data = PAEdata)
    PAEdata$lmerpredlin_PAESex <- predict(lmer_volume_lin_PAESexRAND)
    etasqSex <- effectsize::eta_squared(lmer_volume_lin_PAESex)
    finalPAESex <- lmer_volume_lin_PAESex
    
    #Subsetting PAEdata by Sex for plotting
    M_data <-subset(PAEdata,Sex.F1_M0.=='0')
    F_data <-subset(PAEdata,Sex.F1_M0.=='1')
    M_data_regions <- as.data.frame(M_data[,c(9:125)])
    F_data_regions <- as.data.frame(F_data[,c(9:125)])
    ROIM <- M_data_regions[,i]
    ROIF <- F_data_regions[,i]

    
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
    
    #Plot the group difference for PAEvPS with the current ROI
    plotROI <- ggplot(Alldata, aes(x=Age, y=ROIPAEvPS, group=Subject)) 
    TDplot <- plotROI + geom_point(data=PSdata, aes(x=Age,y=ROI),color='green',size=.6) + 
      geom_line(data=PSdata, aes(x=Age,y=lmerpredlin_PAE), color="green", size=.2)
    PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=ROIPAE),color='purple',size=.6) + 
      geom_smooth(data=PAEdata, aes(x=Age, y=ROIPAE), method=lm, se=FALSE, color="purple", size=.2) + 
      geom_function(fun = ROIfnlin_PAE, color="purple3", size=2) + geom_function(fun = ROIfnlin_TD, color="green3", size=2) +
      theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    #Plot the group difference for Sex in the PAE group with the current ROI
    plotROI_Sex <- ggplot(PAEdata, aes(x=Age, y=ROIPAE, group=Subject)) 
    Mplot <- plotROI_Sex + geom_point(data=M_data, aes(x=Age,y=ROIM),color='blue',size=.6) + 
      geom_line(data=M_data, aes(x=Age,y=lmerpredlin_PAESex), color="blue", size=.2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='red',size=.6) + 
      geom_line(data=F_data, aes(x=Age,y=lmerpredlin_PAESex), color="red", size=.2) + 
      geom_function(fun = ROIfnlin_Fem, color="red3", size=2) + geom_function(fun = ROIfnlin_Male, color="blue3", size=2) +
      theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
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
    regionsPAEvPS <- as.data.frame(Alldata[,c(9:125)])
    ROIPAEvPS <- regionsPAEvPS[,i]
    lmer_volume_quad_PAE <- lmer(ROIPAEvPS ~ I(Age^2)*PAE1_TD0 + Age*PAE1_TD0 + (1|Subject), REML = FALSE, data = Alldata)
    Alldata$lmerpredquad_PAE <- predict(lmer_volume_quad_PAE)
    etasqPAE <- effectsize::eta_squared(lmer_volume_quad_PAE)
    finalPAE <- lmer_volume_quad_PAE
    
    PAEdata <- subset(Alldata, PAE1_TD0=="1")
    PSdata <- subset(Alldata, PAE1_TD0=="0")
    
    #Run a 2nd order linear model to test main and int. effects of Sex within the PAE group
    regionsPAE <- as.data.frame(PAEdata[,c(9:125)])
    ROIPAE <- regionsPAE[,i]
    lmer_volume_quad_PAESex <- lmer(ROIPAE ~ I(Age^2)*Sex.F1_M0. + Age*Sex.F1_M0. + (1|Subject), REML = FALSE, data = PAEdata)
    PAEdata$lmerpredquad_PAESex <- predict(lmer_volume_quad_PAESex)
    etasqSex <- effectsize::eta_squared(lmer_volume_quad_PAESex)
    finalPAESex <- lmer_volume_quad_PAESex
    
    
    M_data <-subset(PAEdata,Sex.F1_M0.=='0')
    F_data <-subset(PAEdata,Sex.F1_M0.=='1')
    M_data_regions <- as.data.frame(M_data[,c(9:125)])
    F_data_regions <- as.data.frame(F_data[,c(9:125)])
    ROIM <- M_data_regions[,i]
    ROIF <- F_data_regions[,i]
    
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
      geom_line(data=PSdata, aes(x=Age,y=lmerpredquad_PAE), color="green", size=.2)
    PAEplot <- TDplot + geom_point(data=PAEdata, aes(x=Age,y=ROIPAE),color='purple',size=.6) + 
      geom_smooth(data=PAEdata, aes(x=Age, y=ROIPAE), method=lm, se=FALSE, color="purple", size=.2) + 
      geom_function(fun = ROIfnquad_PAE, color="purple3", size=2) + geom_function(fun = ROIfnquad_TD, color="green3", size=2) +
      theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    #Plot the group difference for Sex in the PAE group with the current ROI
    plotROI_Sex <- ggplot(PAEdata, aes(x=Age, y=ROIPAE, group=Subject)) 
    Mplot <- plotROI_Sex + geom_point(data=M_data, aes(x=Age,y=ROIM),color='blue',size=.6) + 
      geom_line(data=M_data, aes(x=Age,y=lmerpredquad_PAESex), color="blue", size=.2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='red',size=.6) + 
      geom_line(data=F_data, aes(x=Age,y=lmerpredquad_PAESex), color="red", size=.2) + 
      geom_function(fun = ROIfnquad_Fem, color="red3", size=2) + geom_function(fun = ROIfnquad_Male, color="blue3", size=2) +
      theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
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
  
  write.csv(df.master.coef.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PSvPAE_AbsVolume_20230427.csv")
  write.csv(df.master.etasq.PAE, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/EtaSqTable_PSvPAE_AbsVolume_20230427.csv")
  write.csv(df.master.coef.Sex, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/ModelsTable_PAESex_AbsVolume_20230427.csv")
  write.csv(df.master.etasq.Sex, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/EtaSqTable_PAESex_AbsVolume_20230427.csv")

  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  plot_list[[i]] <- ggplotGrob(PAEplot)
  plot_list_Sex[[i]] <- ggplotGrob(Fplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_Abs_20230427.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

bigPlot_Sex <- marrangeGrob(grobs=plot_list_Sex, ncol=4, nrow=2)
ggsave("BigPlot_PAE_Sex_Abs_20230427.pdf", bigPlot_Sex, width=50, height=25, limitsize = FALSE)



