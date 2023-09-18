Alldata <- read.csv("PSvPAE_Longitudinal_20230420.csv")

#recoding variables into different data types, for example string into numeric
Alldata$ICV <- as.numeric(Alldata$ICV)
Alldata$Age <- as.numeric(Alldata$Age)
Alldata$Subject <- as.factor(Alldata$Subject)
Alldata$PAE1_TD0 <- as.factor(Alldata$PAE1_TD0)
Alldata$Sex.F1_M0. <- as.factor(Alldata$Sex.F1_M0.)

plot_list <- list()

PAEdata <- subset(Alldata, PAE1_TD0=="1")
PSdata <- subset(Alldata, PAE1_TD0=="0")

regions <- as.data.frame(PSdata[,c(9:125)])
dim(regions)
regionsPAE <- as.data.frame(PAEdata[,c(9:125)])
dim(regionsPAE)

for (n in 1:length(regions)) {
  print(colName)
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[n,1]
  ROI <- as.numeric(regions[,n])
  ROIPAE <- as.numeric(regionsPAE[,n])

  x <- PSdata[,5]
  y <- ROI
  out <- paste(colName,"_LOESS.png", sep="")
  
  QA <- as.data.frame(Quantile.loess(y, x, the.quant = .05,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QB <- as.data.frame(Quantile.loess(y, x, the.quant = .10,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QC <- as.data.frame(Quantile.loess(y, x, the.quant = .25,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QD <- as.data.frame(Quantile.loess(y, x, the.quant = .50,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QE <- as.data.frame(Quantile.loess(y, x, the.quant = .75,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QF <- as.data.frame(Quantile.loess(y, x, the.quant = .90,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))
  QG <- as.data.frame(Quantile.loess(y, x, the.quant = .95,window.size = 30,window.alignment = c("center"), percent.of.overlap.between.two.windows = 0.95))

  LOESS_plot_init <- ggplot(PSdata, aes(x=Age, y=ROI)) + geom_blank(aes(x=Age, y=ROI))
  LOESS_plot_QA <- LOESS_plot_init + geom_line(data=QA, aes(x=x, y=y.loess), linetype=4)
  LOESS_plot_QB <- LOESS_plot_QA + geom_line(data=QB, aes(x=x, y=y.loess), linetype=3)
  LOESS_plot_QC <- LOESS_plot_QB + geom_line(data=QC, aes(x=x, y=y.loess), linetype=2)
  LOESS_plot_QD <- LOESS_plot_QC + geom_line(data=QD, aes(x=x, y=y.loess), linetype=1)
  LOESS_plot_QE <- LOESS_plot_QD + geom_line(data=QE, aes(x=x, y=y.loess), linetype=2)
  LOESS_plot_QF <- LOESS_plot_QE + geom_line(data=QF, aes(x=x, y=y.loess), linetype=3)
  LOESS_plot_QG <- LOESS_plot_QF + geom_line(data=QG, aes(x=x, y=y.loess), linetype=4)
  LOESS_plot_PAE <- LOESS_plot_QG + geom_point(data=PAEdata, aes(x=Age, y=ROIPAE, group=Subject), color="purple3", size=.6) + geom_line(data=PAEdata, aes(x=Age, y=ROIPAE, group=Subject), color="purple3", size=.2)
  
  LOESS_plot_Finish <- LOESS_plot_PAE + theme_classic() +theme (plot.title = element_text(size=20, color="black",face="bold")) + 
    theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
    theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
    coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
  LOESS_plot_Finish

  plot_list[[n]] <- ggplotGrob(LOESS_plot_Finish)
  ggsave(paste0(colName,"_PSvPAE_LOESS_20230428.png"), plot=LOESS_plot_Finish, width=11, height=8, limitsize = FALSE)
  print(colName)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_PSvPAE_LOESS_20230428.png", plot=bigPlot, width=50, height=25, limitsize = FALSE)


  
  