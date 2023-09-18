if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmerTest, gridExtra, grid, ggplot2, sjstats, pbkrtest)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Projects/PSvPAE_Volume/")

df <- read.csv("PSvPAE_Segmentations_LongitudinalRegistrationBeforeNAfter_20230822.csv")
PAEdata_Before <- subset(df, PAE=="1" & AfterLR=="0")
PAEdata_After <- subset(df, PAE=="1" & AfterLR=="1")
PSdata_Before <- subset(df, PAE=="0" & AfterLR=="0")
PSdata_After <- subset(df, PAE=="0" & AfterLR=="1")

PAEdata_After$PAE <- as.factor(PAEdata_After$PAE)
PAEdata_After$Subject <- as.factor(PAEdata_After$Subject)
PAEdata_After$AfterLR <- as.factor(PAEdata_After$AfterLR)
PAEdata_After$SubID <- as.factor(PAEdata_After$SubID)
PAEdata_After$Sex.F1_M0. <- as.factor(PAEdata_After$Sex.F1_M0.)


regions <- as.data.frame(PAEdata_After[,c(9:125)])

plot_list <- list()

for (i in 1:length(regions)) {
  ROIPS <- regions[,i]
  ROIcolNamesAll <- as.data.frame(colnames(regions))
  ROIcolName <- ROIcolNamesAll[i,1]
  
  plot_init <- ggplot(data=PAEdata_After, aes(x=Age, y=ROIPS, group=Subject))
  Pointsplot <- plot_init + geom_point(aes(size=11, color=Sex.F1_M0.))
  Lineplot <- Pointsplot + geom_line(aes(size=5, color=Sex.F1_M0.)) + scale_color_manual(values=c( "blue","red")) 
  Finalplot <- Lineplot + coord_cartesian(xlim=c(2,8)) + ylab(paste(ROIcolName)) +
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + theme(legend.position = "none")
  
  plot_list[[i]] <- ggplotGrob(Finalplot)
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=2, nrow=1)

ggsave("bigPlot_PAE_After_20230822_v4.pdf", bigPlot, width=50, height=25, limitsize = FALSE)





#######################

PAEdata_Before$PAE <- as.factor(PAEdata_Before$PAE)
PAEdata_Before$Subject <- as.factor(PAEdata_Before$Subject)
PAEdata_Before$AfterLR <- as.factor(PAEdata_Before$AfterLR)
PAEdata_Before$SubID <- as.factor(PAEdata_Before$SubID)
PAEdata_Before$Sex.F1_M0. <- as.factor(PAEdata_Before$Sex.F1_M0.)

regions <- as.data.frame(PAEdata_Before[,c(9:125)])

plot_list <- list()

for (i in 1:length(regions)) {
  ROIPS <- regions[,i]
  ROIcolNamesAll <- as.data.frame(colnames(regions))
  ROIcolName <- ROIcolNamesAll[i,1]
  
  plot_init <- ggplot(data=PAEdata_Before, aes(x=Age, y=ROIPS, group=Subject))
  Pointsplot <- plot_init + geom_point(aes(size=11, color=Sex.F1_M0.))
  Lineplot <- Pointsplot + geom_line(aes(size=5, color=Sex.F1_M0.)) + scale_color_manual(values=c( "blue","red"))
  Finalplot <- Lineplot + coord_cartesian(xlim=c(2,8)) + ylab(paste(ROIcolName)) +
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + theme(legend.position = "none")
  
  plot_list[[i]] <- ggplotGrob(Finalplot)
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=2, nrow=1)

ggsave("bigPlot_PAE_B4_20230822_v4.pdf", bigPlot, width=50, height=25, limitsize = FALSE)


##########################

PSdata_After$PAE <- as.factor(PSdata_After$PAE)
PSdata_After$Subject <- as.factor(PSdata_After$Subject)
PSdata_After$AfterLR <- as.factor(PSdata_After$AfterLR)
PSdata_After$SubID <- as.factor(PSdata_After$SubID)
PSdata_After$Sex.F1_M0. <- as.factor(PSdata_After$Sex.F1_M0.)


regions <- as.data.frame(PSdata_After[,c(9:125)])

plot_list <- list()

for (i in 1:length(regions)) {
  ROIPS <- regions[,i]
  ROIcolNamesAll <- as.data.frame(colnames(regions))
  ROIcolName <- ROIcolNamesAll[i,1]
  
  plot_init <- ggplot(data=PSdata_After, aes(x=Age, y=ROIPS, group=Subject))
  Pointsplot <- plot_init + geom_point(aes(size=11, color=Sex.F1_M0.))
  Lineplot <- Pointsplot + geom_line(aes(size=5, color=Sex.F1_M0.)) + scale_color_manual(values=c( "blue","red")) 
  Finalplot <- Lineplot + coord_cartesian(xlim=c(2,8)) + ylab(paste(ROIcolName)) +
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + theme(legend.position = "none")
  
  plot_list[[i]] <- ggplotGrob(Finalplot)
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=2, nrow=1)

ggsave("bigPlot_PS_After_20230822_v4.pdf", bigPlot, width=50, height=25, limitsize = FALSE)



################################
PSdata_Before$PAE <- as.factor(PSdata_Before$PAE)
PSdata_Before$Subject <- as.factor(PSdata_Before$Subject)
PSdata_Before$AfterLR <- as.factor(PSdata_Before$AfterLR)
PSdata_Before$SubID <- as.factor(PSdata_Before$SubID)
PSdata_Before$Sex.F1_M0. <- as.factor(PSdata_Before$Sex.F1_M0.)


regions <- as.data.frame(PSdata_Before[,c(9:125)])

plot_list <- list()

for (i in 1:length(regions)) {
  ROIPS <- regions[,i]
  ROIcolNamesAll <- as.data.frame(colnames(regions))
  ROIcolName <- ROIcolNamesAll[i,1]
  
  plot_init <- ggplot(data=PSdata_Before, aes(x=Age, y=ROIPS, group=Subject))
  Pointsplot <- plot_init + geom_point(aes(size=11, color=Sex.F1_M0.))
  Lineplot <- Pointsplot + geom_line(aes(size=5, color=Sex.F1_M0.)) + scale_color_manual(values=c( "blue","red")) 
  Finalplot <- Lineplot + coord_cartesian(xlim=c(2,8)) + ylab(paste(ROIcolName)) +
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + theme(legend.position = "none")
  
  plot_list[[i]] <- ggplotGrob(Finalplot)
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=2, nrow=1)

ggsave("bigPlot_PS_Before_20230822_v4.pdf", bigPlot, width=50, height=25, limitsize = FALSE)
