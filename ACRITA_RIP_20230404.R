lmer_volume_lin_rand <- lmer(Right_Hippocampus ~ Age + (1 + Age|Subject), REML = FALSE, data = PAEdata)

lmer_volume_lin_rand_Sex <- lmer(Right_Hippocampus ~ Age + Sex.F1_M0. + (1 + Age|Subject), REML = FALSE, data = PAEdata)

Right_Hippocampus_ICVNorm <- (PAEdata$Right_Hippocampus/PAEdata$ICV)*100

lmer_volume_lin_rand_ICV <- lmer(Right_Hippocampus_ICVNorm ~ Age + (1 + Age|Subject), REML = FALSE, data = PAEdata)

lmer_volume_lin_rand_Sex_ICV <- lmer(Right_Hippocampus_ICVNorm ~ Age + Sex.F1_M0. + (1 + Age|Subject), REML = FALSE, data = PAEdata)

CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_rand))
lin_Interc <- as.numeric(CoefDataLin[1,1])
lin_beta_Age <- as.numeric(CoefDataLin[2,1])


ROIfnlin_PAE <- function(x) {
  as.numeric(lin_beta_Age*x + lin_Interc) # outcome is y
}

#Plot of PAE only abs
plotROI <- ggplot(PAEdata, aes(x=Age, y=Right_Hippocampus, group=Subject)) 
plotROI1 <- plotROI + geom_point(color='plum4',size=.4) + geom_line(aes(y=predict(lmer_volume_lin_rand)), color="plum4", size=.2)
plotROI2 <- plotROI1 + geom_function(fun=ROIfnlin_PAE, color="plum4", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Hippocampus ("~mm^3*")")



CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_rand_ICV))
lin_Interc <- as.numeric(CoefDataLin[1,1])
lin_beta_Age <- as.numeric(CoefDataLin[2,1])


ROIfnlin_PAE <- function(x) {
  as.numeric(lin_beta_Age*x + lin_Interc) # outcome is y
}

#Plot of PAE only norm
plotROI <- ggplot(PAEdata, aes(x=Age, y=Right_Hippocampus_ICVNorm, group=Subject)) 
plotROI1 <- plotROI + geom_point(color='plum4',size=.4) + geom_line(aes(y=predict(lmer_volume_lin_rand_ICV)), color="plum4", size=.2)
plotROI2 <- plotROI1 + geom_function(fun=ROIfnlin_PAE, color="plum4", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Hippocampus (%ICV)")




#############

fnPS <- function(x) {
  as.numeric(-25.7*I(x^2) + 352.11*x + 2625.68) # outcome is y
}

#Plot of PAE with PS group line abs
plotROI <- ggplot(PAEdata, aes(x=Age, y=Right_Hippocampus, group=Subject)) 
plotROI1 <- plotROI + geom_point(color='plum4',size=.4) + geom_line(aes(y=predict(lmer_volume_lin_rand)), color="plum4", size=.2)
plotROI2 <- plotROI1 + geom_function(fun=fnPS, color="green", size=2) +
  geom_function(fun=ROIfnlin_PAE, color="plum4", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Hippocampus ("~mm^3*")")


fnPS <- function(x) {
  as.numeric(-0.00101211718475235*I(x^2) + 0.011185197*x + 0.263932028) # outcome is y
}

plotROI <- ggplot(PAEdata, aes(x=Age, y=Right_Hippocampus_ICVNorm, group=Subject)) 
plotROI1 <- plotROI + geom_point(color='plum4',size=.4) + geom_line(aes(y=predict(lmer_volume_lin_rand_ICV)), color="plum4", size=.2)
plotROI2 <- plotROI1 + geom_function(fun=fnPS, color="green", size=2) +
  geom_function(fun=ROIfnlin_PAE, color="plum4", size=2) +
  theme_classic() + theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Hippocampus (%ICV)")
