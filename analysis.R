rusin<-read.csv('rusin_mouthte.csv')
AbneyT<-read.csv('ToiletTE.csv')
AbneyF<-read.csv('FingerLip_TE.csv')


require(ggplot2)
require(ggpubr)
require(dplyr)

#distribution of TE per organism
windows()
ggplot(rusin)+geom_violin(aes(x=Organism,y=TE,fill=Organism),draw_quantiles = c(0.25,0.5,0.75),alpha=0.7)+
  scale_fill_discrete(name="",labels=c(expression(italic('M. luteus')),"PRD-1",expression(italic('S. rubidea'))))+
  scale_x_discrete(name="",labels=c(expression(italic('M. luteus')),"PRD-1",expression(italic('S. rubidea'))))+
  scale_y_continuous(name="Transfer Efficiency (%)")+
  theme_pubr()+
  theme(legend.position = "none")

#standard deviations

#sd of M. Luteus
sd(rusin$TE[rusin$Organism=="M. Luteus"])

#sd of PRD-1
sd(rusin$TE[rusin$Organism=="PRD-1"])

#sd of S. rubidea
sd(rusin$TE[rusin$Organism=="S. Rubidea"])

  #Summary 
  summary(rusin)
  # Compute the analysis of variance
  R.aov <- aov(TE ~ Organism, data = rusin)
  # Summary of the analysis
  summary(R.aov)
  #Post-Hoc
  tukey.R.aov<-TukeyHSD(R.aov)
  
  tukey.R.aov

#------------------------------------------------------


#distribution of TE per INNOCULUM
windows()
ggplot(AbneyT)+geom_violin(aes(x=Innoculum,y=TE,fill=Innoculum),draw_quantiles = c(0.25,0.5,0.75),alpha=0.7)+
  scale_fill_discrete(name="",labels=c("PBS", "TSB", "FBS", "Tripartite"))+
  scale_x_discrete(name="",labels=c("PBS", "TSB", "FBS", "Tripartite"))+
  scale_y_continuous(name="Transfer Efficiency (%)")+
  theme_pubr()+
  theme(legend.position = "none")

#standard deviations

#sd of PBS
sd(AbneyT$TE[AbneyT$Innoculum=="PBS"])

#sd of TSB
sd(AbneyT$TE[AbneyT$Innoculum=="TSB"])

#sd of FBS
sd(AbneyT$TE[AbneyT$Innoculum=="FBS"])

#sd of Tripartite
sd(AbneyT$TE[AbneyT$Innoculum=="Tripartite"])

  #Summary 
  summary(AbneyT)
  # Compute the analysis of variance
  T.aov <- aov(TE ~ Innoculum, data = AbneyT)
  # Summary of the analysis
  summary(T.aov)
  #Post-Hoc
  tukey.T.aov<-TukeyHSD(T.aov)
  
  tukey.T.aov

#-----------------------------------------------------------------


#distribution of TE per INNOCULUM
windows()
ggplot(AbneyF)+geom_violin(aes(x=Innoculum,y=TE,fill=Innoculum),draw_quantiles = c(0.25,0.5,0.75),alpha=0.7)+
  scale_fill_discrete(name="",labels=c("PBS", "TSB", "FBS", "Tripartite"))+
  scale_x_discrete(name="",labels=c("PBS", "TSB", "FBS", "Tripartite"))+
  scale_y_continuous(name="Transfer Efficiency (%)")+
  theme_pubr()+
  theme(legend.position = "none")

#standard deviations

#sd of PBS
sd(AbneyF$TE[AbneyF$Innoculum=="PBS"])

#sd of TSB
sd(AbneyF$TE[AbneyF$Innoculum=="TSB"])

#sd of FBS
sd(AbneyF$TE[AbneyF$Innoculum=="FBS"])

#sd of Tripartite
sd(AbneyF$TE[AbneyF$Innoculum=="Tripartite"])

  #Summary 
  summary(AbneyF)
  # Compute the analysis of variance
  F.aov <- aov(TE ~ Innoculum, data = AbneyF)
  # Summary of the analysis
  summary(F.aov)
  #Post-Hoc
  tukey.F.aov<-TukeyHSD(F.aov)
  
  tukey.F.aov

