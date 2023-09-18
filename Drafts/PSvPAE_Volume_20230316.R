###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, rstatix, sjstats, pbkrtest, lmertest, misty)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/Archive/")


PAEdata <- as.data.frame(read.csv("PAEdata.csv"))
PSdata <- read.csv("PSvPAE_PSonly_Longitudinal_January2023.csv")
PAEPSdata <- read.csv("PSvPAE_Longitudinal_RetainedRegions_January2023.csv")

#recoding variables into different data types, for example string into numeric
PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$ICV <- PSdata$ICV/10
PSdata$Age <- as.numeric(PSdata$Age)
PSdata$Subject <- as.factor(PSdata$Subject)
#PSdata$Sex.F1_M0. <- as.factor(PSdata$Sex.F1_M0.)

PAEdata$ICV <- as.numeric(PAEdata$ICV)
PAEdata$ICV <- PAEdata$ICV/10
PAEdata$Age <- as.numeric(PAEdata$Age)
PAEdata$Subject <- as.factor(PAEdata$Subject)
PAEdata$Sex.F1_M0. <- as.factor(PAEdata$Sex.F1_M0.)
PAEdata$Total_Gray <- PAEdata$Total_Gray/10

PAEPSdata$ICV <- as.numeric(PAEPSdata$ICV)
PAEPSdata$ICV <- PAEPSdata$ICV/10
PAEPSdata$Total_Gray <- PAEPSdata$Total_Gray/10
PAEPSdata$Age <- as.numeric(PAEPSdata$Age)
PAEPSdata$Subject <- as.factor(PAEPSdata$Subject)
PAEPSdata$PAE1_TD0 <- as.factor(PAEPSdata$PAE1_TD0)


###### Running PS only to retain non-null regions

regions <-as.data.frame(PSdata[,c(9:125)])
dim(regions)
magic_for(silent = TRUE)
lmer_tableinit <- lmer(Right_TTG_transverse_temporal_gyrus ~ Age + I(Age^2) + ICV + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  #plotfilename <- paste0("NoCov_", colnames(regions)[i])
  finalmodel <- 0
  bestAIC <- 0 
  finalplot <-0
  
  
  #null model
  lmer_volume_null <- lmer(ROI ~ ICV + (1|Subject), REML = FALSE, data = PSdata)
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
    finalmodel <- lmer_volume_lin
    
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null

  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  df.coef <- as.data.frame(coef(summary(finalmodel)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/PS_Table_RegionalTrajectories_NullICV.csv")
}

PS_trajectories <- magic_result_as_dataframe()
write.csv(PS_trajectories, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/PS_Trajectories_AIC_NullICV.csv")
magic_free()


############# Running PAE only to retain non-null regions, also running the sex models for everythin non-null
######
sink("PAE_linearOnly_randslopepossible")
regions <-as.data.frame(PAEdata[,c(9:12)])
print(dim(regions))
#magic_for(silent = TRUE)
lmer_tableinit <- lmer(Right_TTG_transverse_temporal_gyrus ~ Age + I(Age^2) + ICV + (1|Subject), REML = FALSE, data = PAEdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")
plot_list <- list()

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  finalplot <-0
  
  
  #null model
  lmer_volume_null <- lmer(ROI ~ ICV + (1|Subject), REML = FALSE, data = PAEdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + ICV + (1|Subject), REML = FALSE, data = PAEdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + ICV + (1 + mc_age|Subject), REML = FALSE, data = PAEdata)
  AIC_linrand <- AIC(lmer_volume_lin_rand)
  
  
  print(colnames(regions)[i])
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    print("Best Model is Linear, fixed slope")
    
    lmer_volume_lin_Sex <- lmer(ROI ~ Age*Sex.F1_M0. + ICV + (1|Subject), REML = FALSE, data = PAEdata)
    
    colNamesAll <- as.data.frame(colnames(regions))
    colName <- colNamesAll[i,1]
    
    F_data <-subset(PAEdata,Sex.F1_M0.=='1')
    M_data <-subset(PAEdata,Sex.F1_M0.=='0')
    
    F_data_regions <- as.data.frame(F_data[,c(9:12)])
    ROI_F <- F_data_regions[,i]
    
    M_data_regions <- as.data.frame(M_data[,c(9:12)])
    ROI_M <- M_data_regions[,i]
    
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Sex))
    lin_Interc <- as.numeric(CoefDataLin[1,1])
    lin_beta_Age <- as.numeric(CoefDataLin[2,1])
    lin_beta_Sex <- as.numeric(CoefDataLin[3,1])
    lin_beta_ICV <- as.numeric(CoefDataLin[4,1])
    lin_beta_AgebySex <- as.numeric(CoefDataLin[5,1])
    
    final_etasq <- eta_sq(lmer_volume_lin_Sex)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(lmer_volume_lin_Sex)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    lin_IntercF <- lin_Interc+lin_beta_Sex*1
    lin_LinF <- lin_beta_Age+lin_beta_AgebySex*1
    
    lin_IntercM <- lin_Interc+lin_beta_Sex*0
    lin_LinM <- lin_beta_Age+lin_beta_AgebySex*0
    
    
    ## x = Age and b = ICV
    ROIfnlin_F <- function(x,b) {
      as.numeric(lin_LinF*x + lin_IntercF + lin_beta_ICV*b) # outcome is y
    }
    
    ROIfnlin_M <- function(x,b) {
      as.numeric(lin_LinM*x + lin_IntercM + lin_beta_ICV*b) # outcome is y
    }
    
    
    plotROI <- ggplot(PAEdata, aes(x=Age, y=ROI, group=Subject)) 
    Fplot <- plotROI + geom_point(data=F_data, aes(x=Age,y=ROI_F),color='red',size=.4) + geom_line(data=F_data, aes(x=Age,y=ROI_F), color="red", size=.2)
    Mplot <- Fplot + geom_point(data=M_data, aes(x=Age,y=ROI_M),color='blue',size=.4) + geom_line(data=M_data, aes(x=Age,y=ROI_M), color="blue", size=.2) 
    Mplot_groupline <- Mplot + geom_function(fun=ROIfnlin_M, args=list(b=127597.2), color="blue", size=2)
    Fplot_groupline <- Mplot_groupline + geom_function(fun=ROIfnlin_F,args=list(b=117505.1), color="red", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))

  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    print("Best Model is Null")
    #Fplot_groupline <- ggplot(PAEdata, aes(x=Age, y=ROI, group=Subject))
    
  }
  if (bestAIC >= AIC_linrand) {
    bestAIC <- AIC_linrand
    finalmodel <- lmer_volume_lin_rand
    print("Best Model is Linear, random slope")
    
    lmer_volume_quad_Sex <- lmer(ROI ~ Age*Sex.F1_M0. + I(Age^2)*Sex.F1_M0. + ICV + (1|Subject), REML = FALSE, data = PAEdata)
    
    CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Sex))
    quad_Interc <- as.numeric(CoefDataQuad[1,1])
    quad_beta_Age <- as.numeric(CoefDataQuad[2,1])
    quad_beta_Sex <- as.numeric(CoefDataQuad[3,1])
    quad_beta_Age2 <- as.numeric(CoefDataQuad[4,1])
    quad_beta_ICV <- as.numeric(CoefDataQuad[5,1])
    quad_beta_AgebySex <- as.numeric(CoefDataQuad[6,1])
    quad_beta_Age2bySex <- as.numeric(CoefDataQuad[7,1])

    colNamesAll <- as.data.frame(colnames(regions))
    colName <- colNamesAll[i,1]
    
    F_data <-subset(PAEdata,Sex.F1_M0.=='1')
    M_data <-subset(PAEdata,Sex.F1_M0.=='0')
    
    F_data_regions <- as.data.frame(F_data[,c(9:12)])
    ROI_F <- F_data_regions[,i]
    
    M_data_regions <- as.data.frame(M_data[,c(9:12)])
    ROI_M <- M_data_regions[,i]
    
    final_etasq <- eta_sq(lmer_volume_quad_Sex)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(lmer_volume_lin_Sex)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    quad_IntercF <- quad_Interc+quad_beta_Sex*1
    quad_LinAgeF <- quad_beta_Age+quad_beta_AgebySex*1
    quad_QuadAgeF <- quad_beta_Age2+quad_beta_Age2bySex*1
    quad_ICVMeanF <- quad_beta_ICV*117505.1
    
    quad_IntercM <- quad_Interc+quad_beta_Sex*0
    quad_LinAgeM <- quad_beta_Age+quad_beta_AgebySex*0
    quad_QuadAgeM <- quad_beta_Age2+quad_beta_Age2bySex*0
    quad_ICVMeanM <- quad_beta_ICV*127597.2
    
    
    ## x = Age and b = ICV
    ROIfnquad_F <- function(x) {
      as.numeric(quad_LinAgeF*x + quad_QuadAgeF*x^2 + quad_IntercF + quad_ICVMeanF) # outcome is y
    }
    
    ROIfnquad_M <- function(x) {
      as.numeric(quad_LinAgeM*x + quad_QuadAgeM*x^2 + quad_IntercM + quad_ICVMeanM) # outcome is y
    }
    
    plotROI <- ggplot(PAEdata, aes(x=Age, y=ROI, group=Subject)) 
    Fplot <- plotROI + geom_point(data=F_data, aes(x=Age,y=ROI_F),color='red',size=.4) + geom_line(data=F_data, aes(x=Age,y=ROI_F), color="red", size=.2)
    Mplot <- Fplot + geom_point(data=M_data, aes(x=Age,y=ROI_M),color='blue',size=.4) + geom_line(data=M_data, aes(x=Age,y=ROI_M), color="blue", size=.2) 
    Mplot_groupline <- Mplot + geom_function(fun=ROIfnquad_M, color="blue", size=2)
    Fplot_groupline <- Mplot_groupline + geom_function(fun=ROIfnquad_F, color="red", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  df.coef <- as.data.frame(coef(summary(finalmodel)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
  
  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  #write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/Table_PAE_Sex_models_20230130.csv")
  #write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/Table_PAE_Sex_etasq_20230130.csv")
  #plot_list[[i]] <- ggplotGrob(Fplot_groupline)
}
sink()
#bigPlot_PAE_Sex <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
#ggsave("BigPlot_PAE_Sex_20230130.pdf", bigPlot_PAE_Sex, width=50, height=25, limitsize = FALSE)

#PAE_trajectories <- magic_result_as_dataframe()
#write.csv(PAE_trajectories, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/PSvPAE_Volume/PAE_Trajectories_AIC_NullICV.csv")
#magic_free()

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################












###### PAE vs PS regions where linear was most complex model ##################
regions <-as.data.frame(PAEPSdata[,c(9:39)])
dim(regions)
plot_list <- list()

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel_Group <- 0
  print(colnames(regions)[i])
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  PS_data <-subset(PAEPSdata,PAE1_TD0=='0')
  PAE_data <-subset(PAEPSdata,PAE1_TD0=='1')
  
  PS_data_regions <- as.data.frame(PS_data[,c(9:39)])
  ROI_PS <- PS_data_regions[,i]
  
  PAE_data_regions <- as.data.frame(PAE_data[,c(9:39)])
  ROI_PAE <- PAE_data_regions[,i]

  lmer_volume_lin_Group <- lmer(ROI ~ Age*PAE1_TD0 + ICV + (1|Subject), REML = FALSE, data = PAEPSdata)
  CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Group))
  lin_Interc <- as.numeric(CoefDataLin[1,1])
  lin_beta_Age <- as.numeric(CoefDataLin[2,1])
  lin_beta_Group <- as.numeric(CoefDataLin[3,1])
  lin_beta_ICV <- as.numeric(CoefDataLin[4,1])
  lin_beta_AgebyGroup <- as.numeric(CoefDataLin[5,1])

  final_etasq <- eta_sq(lmer_volume_lin_Group)
    
  df.etasq <- as.data.frame(final_etasq)
  df.etasq['Region'] <- c(colnames(regions)[i])
  df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
  df.coef <- as.data.frame(coef(summary(lmer_volume_lin_Group)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
    
  lin_IntercPS <- lin_Interc+lin_beta_Group*0
  lin_LinPS <- lin_beta_Age+lin_beta_AgebyGroup*0
    
  lin_IntercPAE <- lin_Interc+lin_beta_Group*1
  lin_LinPAE <- lin_beta_Age+lin_beta_AgebyGroup*1

  
## x = Age and b = ICV
  ROIfnlin_PS <- function(x,b) {
    as.numeric(lin_LinPS*x + lin_IntercPS + lin_beta_ICV*b) # outcome is y
  }

  ROIfnlin_PAE <- function(x,b) {
    as.numeric(lin_LinPAE*x + lin_IntercPAE + lin_beta_ICV*b) # outcome is y
  }

    
  plotROI <- ggplot(PAEPSdata, aes(x=Age, y=ROI, group=Subject)) 
  PSplot <- plotROI + geom_point(data=PS_data, aes(x=Age,y=ROI_PS),color='green3',size=.4) + geom_line(data=PS_data, aes(x=Age,y=ROI_PS), color="green", size=.2)
  PAEplot <- PSplot + geom_point(data=PAE_data, aes(x=Age,y=ROI_PAE),color='plum4',size=.4) + geom_line(data=PAE_data, aes(x=Age,y=ROI_PAE), color="plum4", size=.2) 
  PAEplot_groupline <- PAEplot + geom_function(fun=ROIfnlin_PS, args=list(b=126119.3), color="green3", size=2)
  PSplot_groupline <- PAEplot_groupline + geom_function(fun=ROIfnlin_PAE,args=list(b=122310.9), color="plum4", size=2) +
    theme_classic() +
    theme (plot.title = element_text(size=20, color="black",face="bold")) + 
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
    coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))

  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "PSvPAE_Table_finalmodelsummary_linear_20230130.csv")
  write.csv(df.master.etasq, "PSvPAE_Table_finalmodeletasq_linear_20230130.csv")
  plot_list[[i]] <- ggplotGrob(PSplot_groupline)
}
bigPlot_PSvPAE <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_linear_20230130_.pdf", bigPlot_PSvPAE, width=50, height=25, limitsize = FALSE)


###### PAE vs PS regions where quadratic was most complex model ##################
regions <-as.data.frame(PAEPSdata[,c(40:111)])
dim(regions)
plot_list <- list()

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel_Group <- 0
  print(colnames(regions)[i])
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  PS_data <-subset(PAEPSdata,PAE1_TD0=='0')
  PAE_data <-subset(PAEPSdata,PAE1_TD0=='1')
  
  PS_data_regions <- as.data.frame(PS_data[,c(40:111)])
  ROI_PS <- PS_data_regions[,i]
  
  PAE_data_regions <- as.data.frame(PAE_data[,c(40:111)])
  ROI_PAE <- PAE_data_regions[,i]
  
  lmer_volume_quad_Group <- lmer(ROI ~ Age*PAE1_TD0 + I(Age^2)*PAE1_TD0 + ICV + (1|Subject), REML = FALSE, data = PAEPSdata)
  CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Group))
  quad_Interc <- as.numeric(CoefDataQuad[1,1])
  quad_beta_Age <- as.numeric(CoefDataQuad[2,1])
  quad_beta_Group <- as.numeric(CoefDataQuad[3,1])
  quad_beta_Age2 <- as.numeric(CoefDataQuad[4,1])
  quad_beta_ICV <- as.numeric(CoefDataQuad[5,1])
  quad_beta_AgebyGroup <- as.numeric(CoefDataQuad[6,1])
  quad_beta_Age2byGroup <- as.numeric(CoefDataQuad[7,1])
  
  final_etasq <- eta_sq(lmer_volume_quad_Group)
  
  df.etasq <- as.data.frame(final_etasq)
  df.etasq['Region'] <- c(colnames(regions)[i])
  df.master.etasq <- rbind(df.master.etasq, df.etasq)
  
  df.coef <- as.data.frame(coef(summary(lmer_volume_quad_Group)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
  
  quad_IntercPS <- quad_Interc+quad_beta_Group*0
  quad_LinAgePS <- quad_beta_Age+quad_beta_AgebyGroup*0
  quad_QuadAgePS <- quad_beta_Age2+quad_beta_Age2byGroup*0
  quad_ICVMeanPS <- quad_beta_ICV*126119.3
  
  quad_IntercPAE <- quad_Interc+quad_beta_Group*1
  quad_LinPAE <- quad_beta_Age+quad_beta_AgebyGroup*1
  quad_QuadPAE <- quad_beta_Age2+quad_beta_Age2byGroup*1
  quad_ICVMeanPAE <- quad_beta_ICV*122310.9
  
  
  ## x = Age and b = mean ICV for the group
  ROIfnquad_PS <- function(x) {
    as.numeric(quad_LinPS*x + quad_QuadPS*x^2 + quad_IntercPS + quad_ICVMeanPS) # outcome is y
  }
  
  ROIfnquad_PAE <- function(x) {
    as.numeric(quad_LinPAE*x + quad_QuadPAE*x^2 + quad_IntercPAE + quad_ICVMeanPAE) # outcome is y
  }
  
  
  plotROI <- ggplot(PAEPSdata, aes(x=Age, y=ROI, group=Subject)) 
  PSplot <- plotROI + geom_point(data=PS_data, aes(x=Age,y=ROI_PS),color='green3',size=.4) + geom_line(data=PS_data, aes(x=Age,y=ROI_PS), color="green", size=.2)
  PAEplot <- PSplot + geom_point(data=PAE_data, aes(x=Age,y=ROI_PAE),color='plum4',size=.4) + geom_line(data=PAE_data, aes(x=Age,y=ROI_PAE), color="plum4", size=.2) 
  PAEplot_groupline <- PAEplot + geom_function(fun=ROIfnquad_PS, color="green3", size=2)
  PSplot_groupline <- PAEplot_groupline + geom_function(fun=ROIfnquad_PAE, color="plum4", size=2) +
    theme_classic() +
    theme (plot.title = element_text(size=20, color="black",face="bold")) + 
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
    coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "PSvPAE_Table_finalmodelsummary_quadratic_20230130.csv")
  write.csv(df.master.etasq, "PSvPAE_Table_finalmodeletasq_quadratic_20230130.csv")
  plot_list[[i]] <- ggplotGrob(PSplot_groupline)
}
bigPlot_PSvPAE <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_quadratic_20230130.pdf", bigPlot_PSvPAE, width=50, height=25, limitsize = FALSE)


########
#trying geom_ribbon

PS_data <-subset(PAEPSdata,PAE1_TD0=='0')
PAE_data <-subset(PAEPSdata,PAE1_TD0=='1')

lmer_volume_quad <- lmer(Left_Thalamus_Proper ~ Age + I(Age^2) + ICV + (1|Subject), REML = FALSE, data = PSdata)
confint(lmer_volume_quad, oldNames=FALSE)

CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad))
QuadCoQuad <- as.numeric(CoefDataQuad[3,1])
LinCoQuad <- as.numeric(CoefDataQuad[2,1])
IntercQuad <- as.numeric(CoefDataQuad[1,1])

ROIfnquad <- function(x) {
  as.numeric(QuadCoQuad*I(x^2) + LinCoQuad*x + IntercQuad) # outcome is y
}

PAE_data$Subject <- as.factor(PAE_data$Subject)

mainline <- ggplot(PAE_data, aes(x=Age, y=Left_Thalamus_Proper, group=Subject)) + geom_point(data=PAE_data, aes(x=Age, y=Left_Thalamus_Proper)) + geom_line(data=PAE_data, aes(x=Age, y=Left_Thalamus_Proper)) + theme_classic()
  
PS_Overlay <- mainline + geom_smooth(data=PAEPSdata, aes(x=Age, y=Left_Thalamus_Proper))
PAE_Overlay
