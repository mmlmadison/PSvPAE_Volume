###Looking for sex differences within the PAE group
###Madison Long 2023-03-06

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, sjstats, pkbrtest)
devtools::install_github("cardiomoon/ggiraphExtra")
devtools::install_github("hoxo-m/magicfor")
library(magicfor)

install.packages('lmerTest')
library(lmerTest)
install.packages('pbkrtest')
library(pbkrtest)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/")

PSdata <- read.csv("PSvPAE_Longitudinal_Nov2022.csv")

#recoding variables into differenty data types, for example string into numeric
PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Age <- as.numeric(PSdata$Age)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$Sex.F1_M0. <- as.factor(PSdata$Sex.F1_M0.)
PSdata$PAE1_TD0 <- as.factor(PSdata$PAE1_TD0)


###### PAE vs PS without sex included
regions <-as.data.frame(PSdata[,c(9:12)])
dim(regions)
plot_list <- list()

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

PAE_data <-subset(PSdata,PAE1_TD0=='1')

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel_Age <- 0
  finalmodel_Group <- 0
  bestAIC <- 0 
  
  
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  PAE_data_M <-subset(PAE_data, Sex.F1_M0.=='0')
  PAE_data_F <-subset(PAE_data, Sex.F1_M0.=='1')
  
  PAEM_data_regions <- as.data.frame(PAE_data_M[,c(9:12)])
  ROI_PAEM <- PAEM_data_regions[,i]
  
  PAEF_data_regions <- as.data.frame(PAE_data_F[,c(9:12)])
  ROI_PAEF <- PAEF_data_regions[,i]

  
  
  #null model
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PAE_data)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + (1|Subject), REML = FALSE, data = PAE_data)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PAE_data)
  AIC_quad <- AIC(lmer_volume_quad)
  
  
  print(colnames(regions)[i])
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    lmer_volume_lin_Group <- lmer(ROI ~ Age*PAE1_TD2 + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Group))
    lin_Interc <- as.numeric(CoefDataLin[1,1])
    lin_beta_Age <- as.numeric(CoefDataLin[2,1])
    lin_beta_Group <- as.numeric(CoefDataLin[3,1])
    lin_beta_Sex <- as.numeric(CoefDataLin[4,1])
    lin_beta_AgebyGroup <- as.numeric(CoefDataLin[5,1])
    lin_beta_AgebySex <- as.numeric(CoefDataLin[6,1])
    lin_beta_GroupbySex <- as.numeric(CoefDataLin[7,1])
    lin_beta_AgebyGroupbySex <- as.numeric(CoefDataLin[7,1])
    finalmodel_Age <- lmer_volume_lin
    
    finalmodel_Group <- lmer_volume_lin_Group
    final_etasq <- eta_sq(finalmodel_Group)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(finalmodel_Group)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    lin_IntercPSM <- lin_Interc+lin_beta_Sex*2+lin_beta_Group*2
    lin_LinPSM <- lin_beta_Age+lin_beta_AgebySex*2+lin_beta_AgebyGroup*2+lin_beta_AgebyGroupbySex*2*2
    
    lin_IntercPSF <- lin_Interc+lin_beta_Sex*1+lin_beta_Group*2
    lin_LinPSF <- lin_beta_Age+lin_beta_AgebySex*1+lin_beta_AgebyGroup*2+lin_beta_AgebyGroupbySex*1*2
    
    lin_IntercPAEM <- lin_Interc+lin_beta_Sex*2+lin_beta_Group*1
    lin_LinPAEM <- lin_beta_Age+lin_beta_AgebySex*2+lin_beta_AgebyGroup*1+lin_beta_AgebyGroupbySex*2*1
    
    lin_IntercPAEF <- lin_Interc+lin_beta_Sex*1+lin_beta_Group*1
    lin_LinPAEF <- lin_beta_Age+lin_beta_AgebySex*1+lin_beta_AgebyGroup*1+lin_beta_AgebyGroupbySex*1*1
    
    
    ROIfnlin_PSM <- function(x) {
      as.numeric(lin_LinPSM*x + lin_IntercPSM) # outcome is y
    }
    
    ROIfnlin_PSF <- function(x) {
      as.numeric(lin_LinPSF*x + lin_IntercPSF) # outcome is y
    }
    
    ROIfnlin_PAEM <- function(x) {
      as.numeric(lin_LinPAEM*x + lin_IntercPAEM) # outcome is y
    }
    
    ROIfnlin_PAEF <- function(x) {
      as.numeric(lin_LinPAEF*x + lin_IntercPAEF) # outcome is y
    }
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, group=Subject)) 
    PSMplot <- plotROI + geom_point(data=PS_data_M, aes(x=Age,y=ROI_PSM),color='dark green',size=.4) + geom_line(data=PS_data_M, aes(x=Age,y=ROI_PSM), color="dark green", size=.2)
    PSFplot <- PSMplot + geom_point(data=PS_data_F, aes(x=Age,y=ROI_PSF),color='dark green',size=.4) + geom_line(data=PS_data_F, aes(x=Age,y=ROI_PSF), color="dark green", size=.2, linetype=2) 
    PAEMplot <- PSFplot + geom_point(data=PAE_data_M, aes(x=Age,y=ROI_PAEM),color='plum3',size=.4) + geom_line(data=PAE_data_M, aes(x=Age,y=ROI_PAEM), color="plum3", size=.2) 
    PAEFplot <- PAEMplot + geom_point(data=PAE_data_F, aes(x=Age,y=ROI_PAEF),color='plum3',size=.4) + geom_line(data=PAE_data_F, aes(x=Age,y=ROI_PAEF), color="plum3", size=.2, linetype=2) +
      geom_function(fun = ROIfnlin_PSM, color="dark green", size=2) +
      geom_function(fun = ROIfnlin_PSF, color="dark green", size=2, linetype=2) +
      geom_function(fun = ROIfnlin_PAEM, color="plum3", size=2) +
      geom_function(fun = ROIfnlin_PAEF, color="plum3", size=2, linetype=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    
  }else {
    bestAIC <- AIC_null
    finalmodel_Age <- lmer_volume_null
    finalmodel_Group <- 0
    Fplot <- ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    lmer_volume_quad_Group <- lmer(ROI ~ Age*PAE1_TD2*Sex.F1_M2. + I(Age^2)*PAE1_TD2*Sex.F1_M2.+ (1|Subject), REML = FALSE, data = PSdata)
    CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Group))
    quad_Interc <- as.numeric(CoefDataQuad[1,1])
    quad_beta_Age <- as.numeric(CoefDataQuad[2,1])
    quad_beta_Group <- as.numeric(CoefDataQuad[3,1])
    quad_beta_Sex <- as.numeric(CoefDataQuad[4,1])
    quad_beta_Age2 <- as.numeric(CoefDataQuad[5,1])
    quad_beta_AgebyGroup <- as.numeric(CoefDataQuad[6,1])
    quad_beta_AgebySex <- as.numeric(CoefDataQuad[7,1])
    quad_beta_GroupbySex <- as.numeric(CoefDataQuad[8,1])
    quad_beta_Age2byGroup <- as.numeric(CoefDataQuad[9,1])
    quad_beta_Age2bySex <- as.numeric(CoefDataQuad[10,1])
    quad_beta_AgebyGroupbySex <- as.numeric(CoefDataQuad[11,1])
    quad_beta_Age2byGroupbySex <- as.numeric(CoefDataQuad[12,1])
    
    finalmodel_Age <- lmer_volume_quad
    finalmodel_Group <- lmer_volume_quad_Group
    
    final_etasq <- eta_sq(finalmodel_Group)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(finalmodel_Group)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    quad_IntercPSM <- quad_Interc+quad_beta_Sex*2+quad_beta_Group*2
    quad_LinPSM <- quad_beta_Age+quad_beta_AgebySex*2+quad_beta_AgebyGroup*2+quad_beta_AgebyGroupbySex*2*2
    quad_QuadPSM <- quad_beta_Age2+quad_beta_Age2bySex*2+quad_beta_Age2byGroup*2+quad_beta_Age2byGroupbySex*2*2
    
    quad_IntercPSF <- quad_Interc+quad_beta_Sex*1+quad_beta_Group*2
    quad_LinPSF <- quad_beta_Age+quad_beta_AgebySex*1+quad_beta_AgebyGroup*2+quad_beta_AgebyGroupbySex*1*2
    quad_QuadPSF <- quad_beta_Age2+quad_beta_Age2bySex*1+quad_beta_Age2byGroup*2+quad_beta_Age2byGroupbySex*1*2
    
    quad_IntercPAEM <- quad_Interc+quad_beta_Sex*2+quad_beta_Group*1
    quad_LinPAEM <- quad_beta_Age+quad_beta_AgebySex*2+quad_beta_AgebyGroup*1+quad_beta_AgebyGroupbySex*2*1
    quad_QuadPAEM <- quad_beta_Age2+quad_beta_Age2bySex*2+quad_beta_Age2byGroup*1+quad_beta_Age2byGroupbySex*2*1
    
    quad_IntercPSM <- quad_Interc+quad_beta_Sex*1+quad_beta_Group*1
    quad_LinPSM <- quad_beta_Age+quad_beta_AgebySex*1+quad_beta_AgebyGroup*2+quad_beta_AgebyGroupbySex*1*1
    quad_QuadPSM <- quad_beta_Age2+quad_beta_Age2bySex*1+quad_beta_Age2byGroup*1+quad_beta_Age2byGroupbySex*1*1
    
    
    ROIfnquad_PSM <- function(x) {
      as.numeric(quad_LinPSM*x + quad_QuadPSM*x^2 + quad_IntercPSM) # outcome is y
    }
    
    ROIfnquad_PSF <- function(x) {
      as.numeric(quad_LinPSF*x + quad_QuadPSF*x^2 + quad_IntercPSF) # outcome is y
    }
    
    ROIfnquad_PAEM <- function(x) {
      as.numeric(quad_LinPAEM*x + quad_QuadPAEM*x^2 + quad_IntercPAEM) # outcome is y
    }
    
    ROIfnquad_PAEF <- function(x) {
      as.numeric(quad_LinPAEF*x + quad_QuadPAEF*x^2 + quad_IntercPAEF) # outcome is y
    }
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, group=Subject)) 
    PSMplot <- plotROI + geom_point(data=PS_data_M, aes(x=Age,y=ROI_PSM),color='dark green',size=.4) + geom_line(data=PS_data_M, aes(x=Age,y=ROI_PSM), color="dark green", size=.2)
    PSFplot <- PSMplot + geom_point(data=PS_data_F, aes(x=Age,y=ROI_PSF),color='dark green',size=.4) + geom_line(data=PS_data_F, aes(x=Age,y=ROI_PSF), color="dark green", size=.2, linetype=2) 
    PAEMplot <- PSFplot + geom_point(data=PAE_data_M, aes(x=Age,y=ROI_PAEM),color='plum3',size=.4) + geom_line(data=PAE_data_M, aes(x=Age,y=ROI_PAEM), color="plum3", size=.2) 
    PAEFplot <- PAEMplot + geom_point(data=PAE_data_F, aes(x=Age,y=ROI_PAEF),color='plum3',size=.4) + geom_line(data=PAE_data_F, aes(x=Age,y=ROI_PAEF), color="plum3", size=.2, linetype=2) +
      geom_function(fun = ROIfnquad_PSM, color="dark green", size=2) +
      geom_function(fun = ROIfnquad_PSF, color="dark green", size=2, linetype=2) +
      geom_function(fun = ROIfnquad_PAEM, color="plum3", size=2) +
      geom_function(fun = ROIfnquad_PAEF, color="plum3", size=2, linetype=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
    #print(paste(MaleQuadQuad, MaleLinQuad, MaleIntQuad))
    #print(paste(FemQuadQuad, FemLinQuad, FemIntQuad))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  print("###########")
  #print(summary(finalmodel_Age))
  print("###########")
  print(summary(finalmodel_Group))
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "PSvPAE_Table_finalmodelsummary_trim.csv")
  write.csv(df.master.etasq, "PSvPAE_Table_finalmodeletasq_trim.csv")
  #print(Fplot)
  plot_list[[i]] <- ggplotGrob(PAEFplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_Abs_20221120_.pdf", bigPlot, width=50, height=25, limitsize = FALSE)
