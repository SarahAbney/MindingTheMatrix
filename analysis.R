rusin<-read.csv('rusin_mouthte.csv')
AbneyT<-read.csv('ToiletTE.csv')
AbneyF<-read.csv('FingerLip_TE.csv')
TSB<-read.csv('CombinedTE_TSB.csv')
FBS<-read.csv('CombinedTE_FBS.csv')
PBS<-read.csv('CombinedTE_PBS.csv')
ASTM<-read.csv('CombinedTE_ASTM.csv')
IRN<-read.csv('Abney_noromodel_ALLIR.csv')
IRA<-read.csv('Abney_adenomodelALLIR.csv')

require(ggplot2)
require(ggpubr)
require(ggbreak) #If you use ggbreak in published research, please cite the following paper:S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively utilize plotting space to deal with large datasets and outliers. Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846
require(corrplot)
require(tidyr)
require(truncnorm)
require(triangulr)
require(vtable)
require(dplyr)
require(FSA)
require(effsize)
require(rstatix)
                

#distribution of TE per organism
windows()
   ggplot(rusin)+geom_violin(aes(x=Organism,y=TE,fill=Organism),alpha=0.4,draw_quantiles=c(0.25,0.5,0.75))+
    geom_jitter(aes(x=Organism,y=TE),width=0.01,alpha=0.4)+
    scale_x_discrete(name="",labels=c(expression(italic('M. luteus')),"PRD-1",expression(italic('S. rubidea'))))+
    scale_y_continuous(name="Transfer Efficiency (%)")+
    theme_pubr()
#standard deviations

#sd of M. Luteus
sd(rusin$TE[rusin$Organism=="M. Luteus"])

#sd of PRD-1
sd(rusin$TE[rusin$Organism=="PRD-1"])

#sd of S. rubidea
sd(rusin$TE[rusin$Organism=="S. Rubidea"])

  #Summary 
#  summary(rusin)
  # Compute the analysis of variance
#  R.aov <- aov(TE ~ Organism, data = rusin)
  # Summary of the analysis
#  summary(R.aov)
  #Post-Hoc
#  tukey.R.aov<-TukeyHSD(R.aov)
  
#  tukey.R.aov
  
  kruskal.test(TE ~ Organism, data = rusin) 
  dunnTest(TE ~ Organism, data = rusin, method = "bh")
  
#------------------------------------------------------


#distribution of TE per INNOCULUM
windows()
  ggplot(AbneyT)+geom_violin(aes(x=Inoculum,y=TE,fill=Inoculum),alpha=0.4,draw_quantiles=c(0.25,0.5,0.75))+
    geom_jitter(aes(x=Inoculum,y=TE),width=0.01,alpha=0.4)+
    scale_x_discrete(name="Inoculum")+
    scale_y_continuous(name="Transfer Efficiency (%)")+
    theme_pubr()

#standard deviations

#sd of PBS
sd(AbneyT$TE[AbneyT$Innoculum=="PBS"])

#sd of TSB
sd(AbneyT$TE[AbneyT$Innoculum=="TSB"])

#sd of FBS
sd(AbneyT$TE[AbneyT$Innoculum=="FBS"])

#sd of Tripartite
sd(AbneyT$TE[AbneyT$Innoculum=="ASTM Soil Load"])


#Summary 
  summary(AbneyT)
  # Compute the analysis of variance
    # T.aov <- aov(TE ~ Innoculum, data = AbneyT)
  kruskal.test(TE ~ Inoculum, data = AbneyT)
  kruskal_effsize(TE ~ Inoculum, data = AbneyT)
  # Summary of the analysis
    # summary(T.aov)
  #summary(F2HK.test)
  #Post-Hoc
    # tukey.T.aov<-TukeyHSD(T.aov)
    # tukey.T.aov
  dunnTest(TE ~ Inoculum, data = AbneyT, method = "bh")

 
#-----------------------------------------------------------------


#distribution of TE per INNOCULUM
windows()

  ggplot(AbneyF)+geom_violin(aes(x=Inoculum,y=TE,fill=Inoculum),alpha=0.4,draw_quantiles=c(0.25,0.5,0.75))+
    geom_jitter(aes(x=Inoculum,y=TE),width=0.01,alpha=0.4)+
    scale_x_discrete(name="Inoculum")+
    scale_y_continuous(name="Transfer Efficiency (%)")+
    theme_pubr()

#standard deviations

#sd of PBS
sd(AbneyF$TE[AbneyF$Innoculum=="PBS"])

#sd of TSB
sd(AbneyF$TE[AbneyF$Innoculum=="TSB"])

#sd of FBS
sd(AbneyF$TE[AbneyF$Innoculum=="FBS"])

#sd of Tripartite
sd(AbneyF$TE[AbneyF$Innoculum=="ASTM Soil Load"])

  #Summary 
  summary(AbneyF)
  # Compute the analysis of variance
    #F.aov <- aov(TE ~ Innoculum, data = AbneyF)
  kruskal.test(TE ~ Inoculum, data = AbneyF)
  kruskal_effsize(TE ~ Inoculum , data = AbneyF)
  # Summary of the analysis
  #summary(H2LK.test)
  #Post-Hoc
    #tukey.F.aov<-TukeyHSD(F.aov)
    #tukey.F.aov
  dunnTest(TE ~ Innoculum, data = AbneyF, method = "bh")
  
#-------------------------------------------------------------------
  #STATs to see significance between Toilet and Lip TE as requested in review process
  wilcox_test(TE ~ Innoculum , data = TSB)
    wilcox_effsize(TE ~ Innoculum , data = TSB)
  wilcox_test(TE ~ Innoculum , data = FBS)
    wilcox_effsize(TE ~ Innoculum , data = FBS)
  wilcox_test(TE ~ Innoculum , data = PBS)
    wilcox_effsize(TE ~ Innoculum , data = PBS)
  wilcox_test(TE ~ Innoculum , data = ASTM)
    wilcox_effsize(TE ~ Innoculum , data = ASTM)
    

#--------------------------------------------------------------------#
  #---------------------RISK MODEL Rusin-Lopez---------------#
  #--------------------------------------------------------------------#  
set.seed(5)
  
Fomite.MindingMatrix.LpR<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
#--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency to hand (Lopez) --------------
  te.vir.ts2h<-rtruncnorm(iterations, .186, .55, .3445, .138) #Lopez Averaged (low & High RH) Laminate 2013
  #--------------Transfer efficiency to hand (Rusin) --------------
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.3397,0.1318) #Rusin 2002 assumed distribution
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip 
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_ts<-data.frame(adenovirus.ts = adenovirus.ts,
                       adenovirus.ts.gc = adenovirus.ts.gc,
                       ts2h.c.adenovirus = ts2h.c.adenovirus, 
                       ts2h.adenovirus = ts2h.adenovirus, 
                       ts.adenovirus.ir = ts.adenovirus.ir)
  noro_ts<-data.frame(norovirus.ts = norovirus.ts,
                       norovirus.ts.gc = norovirus.ts.gc,
                       ts2h.c.norovirus = ts2h.c.norovirus, 
                       ts2h.norovirus = ts2h.norovirus, 
                       ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_ts,'LpRusin_adenomodel.csv', row.names = FALSE)
  write.csv(noro_ts, 'LpRusin_noromodel.csv', row.names = FALSE)

  
  #----------Saving Time--------------------
  
  adeno <- read.csv('LpRusin_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'LpR.sum_adeno_model.csv')
  noro <- read.csv('LpRusin_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'LpR.sum_noro_model.csv')
}
Fomite.MindingMatrix.LpR(10000)

#--------------------------------------------------------------------#
#---------------------RISK MODEL PBS ---------------#
#--------------------------------------------------------------------#

Fomite.MindingMatrix.PBS<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
  #--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency to hand (Abney - PBS) --------------
  te.vir.ts2h<-rtruncnorm(iterations, 0, 1, 0.0110, 0.0081) 
  #--------------Transfer efficiency to hand (Abney - PBS) --------------
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.5253,0.0448) 
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip 
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_tsP<-data.frame(adenovirus.ts = adenovirus.ts,
                       adenovirus.ts.gc = adenovirus.ts.gc,
                       ts2h.c.adenovirus = ts2h.c.adenovirus, 
                       ts2h.adenovirus = ts2h.adenovirus, 
                       ts.adenovirus.ir = ts.adenovirus.ir)
  noro_tsP<-data.frame(norovirus.ts = norovirus.ts,
                      norovirus.ts.gc = norovirus.ts.gc,
                      ts2h.c.norovirus = ts2h.c.norovirus, 
                      ts2h.norovirus = ts2h.norovirus, 
                      ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_tsP,'PBSAbney_adenomodel.csv', row.names = FALSE)
  write.csv(noro_tsP, 'PBSAbney_noromodel.csv', row.names = FALSE)
  
  
  #----------Saving Time--------------------
  
  adeno <- read.csv('PBSAbney_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'PBS.sum_adeno_model.csv')
  noro <- read.csv('PBSAbney_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'PBS.sum_noro_model.csv')
}
Fomite.MindingMatrix.PBS(10000)


#--------------------------------------------------------------------#
#---------------------RISK MODEL ASTM ---------------#
#--------------------------------------------------------------------#

Fomite.MindingMatrix.ASTM<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
  #--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency to hand (Abney - ASTM) --------------
  te.vir.ts2h<-rtruncnorm(iterations, 0, 1, 0.0302, 0.0403) 
  #--------------Transfer efficiency to hand (Abney - ASTM) --------------
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.4100,0.1098) 
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip 
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_tsA<-data.frame(adenovirus.ts = adenovirus.ts,
                       adenovirus.ts.gc = adenovirus.ts.gc,
                       ts2h.c.adenovirus = ts2h.c.adenovirus, 
                       ts2h.adenovirus = ts2h.adenovirus, 
                       ts.adenovirus.ir = ts.adenovirus.ir)
  noro_tsA<-data.frame(norovirus.ts = norovirus.ts,
                      norovirus.ts.gc = norovirus.ts.gc,
                      ts2h.c.norovirus = ts2h.c.norovirus, 
                      ts2h.norovirus = ts2h.norovirus, 
                      ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_tsA,'ASTMAbney_adenomodel.csv', row.names = FALSE)
  write.csv(noro_tsA, 'ASTMAbney_noromodel.csv', row.names = FALSE)
  
  
  #----------Saving Time--------------------
  
  adeno <- read.csv('ASTMAbney_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'ASTM.sum_adeno_model.csv')
  noro <- read.csv('ASTMAbney_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'ASTM.sum_noro_model.csv')
}
Fomite.MindingMatrix.ASTM(10000)

#--------------------------------------------------------------------#
#---------------------RISK MODEL TSB ---------------#
#--------------------------------------------------------------------#

Fomite.MindingMatrix.TSB<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
  #--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency to hand (Abney - TSB) --------------
  te.vir.ts2h<-rtruncnorm(iterations, 0, 1, 0.0214, 0.0162) 
  #--------------Transfer efficiency to hand (Abney - TSB) --------------
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.2315,0.2427) 
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip 
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_tsT<-data.frame(adenovirus.ts = adenovirus.ts,
                       adenovirus.ts.gc = adenovirus.ts.gc,
                       ts2h.c.adenovirus = ts2h.c.adenovirus, 
                       ts2h.adenovirus = ts2h.adenovirus, 
                       ts.adenovirus.ir = ts.adenovirus.ir)
  noro_tsT<-data.frame(norovirus.ts = norovirus.ts,
                      norovirus.ts.gc = norovirus.ts.gc,
                      ts2h.c.norovirus = ts2h.c.norovirus, 
                      ts2h.norovirus = ts2h.norovirus, 
                      ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_tsT,'TSBAbney_adenomodel.csv', row.names = FALSE)
  write.csv(noro_tsT, 'TSBAbney_noromodel.csv', row.names = FALSE)
  
  
  #----------Saving Time--------------------
  
  adeno <- read.csv('TSBAbney_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'TSB.sum_adeno_model.csv')
  noro <- read.csv('TSBAbney_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'TSB.sum_noro_model.csv')
}
Fomite.MindingMatrix.TSB(10000)

#--------------------------------------------------------------------#
#---------------------RISK MODEL FBS ---------------#
#--------------------------------------------------------------------#

Fomite.MindingMatrix.FBS<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
  #--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency to hand (Abney - FBS) --------------
  te.vir.ts2h<-rtruncnorm(iterations, 0, 1, 0.0114, 0.0123) 
  #--------------Transfer efficiency to hand (Abney - FBS) --------------
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.3285,0.0983) 
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip  
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_tsF<-data.frame(adenovirus.ts = adenovirus.ts,
                       adenovirus.ts.gc = adenovirus.ts.gc,
                       ts2h.c.adenovirus = ts2h.c.adenovirus, 
                       ts2h.adenovirus = ts2h.adenovirus, 
                       ts.adenovirus.ir = ts.adenovirus.ir)
  noro_tsF<-data.frame(norovirus.ts = norovirus.ts,
                      norovirus.ts.gc = norovirus.ts.gc,
                      ts2h.c.norovirus = ts2h.c.norovirus, 
                      ts2h.norovirus = ts2h.norovirus, 
                      ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_tsF,'FBSAbney_adenomodel.csv', row.names = FALSE)
  write.csv(noro_tsF, 'FBSAbney_noromodel.csv', row.names = FALSE)
  
  
  #----------Saving Time--------------------
  
  adeno <- read.csv('FBSAbney_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'FBS.sum_adeno_model.csv')
  noro <- read.csv('FBSAbney_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'FBS.sum_noro_model.csv')
}
Fomite.MindingMatrix.FBS(10000)

#--------------------------------------------------------------------#
#---------------------RISK MODEL RUSIN ---------------#
#--------------------------------------------------------------------#

Fomite.MindingMatrix.R<-function(iterations){ 
  #Fraction of infectious particles 
  frac.inf<-runif(iterations, 100, 1000)
  #Concentrations on toilet seat
  norovirus.ts.gc<-rtri(iterations, (31217*0.75), (31217*1.25), 31217) #Triangular distributions assumed
  adenovirus.ts.gc<-rtruncnorm(iterations, a=0, 222, 44)
  
  #--------------------Calculating Dose------------------- 
  #organism.ts = fomite concentration on toilet seat (gc/cm2)
  
  norovirus.ts<-norovirus.ts.gc*(1/frac.inf)
  adenovirus.ts<-adenovirus.ts.gc*(1/frac.inf)
  
  #--------------Transfer efficiency (Rusin) --------------
  te.vir.ts2h<-rtruncnorm(iterations, 0, 1, 0.6580, 0.1) #Rusin Phone Reciever (Plastic) 2002 assumed distribution
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.3397,0.1318) #Rusin 2002 assumed distribution
  
  te.vir.h2f<-rtruncnorm(iterations, 0, 1, 0.3285,0.0983) 
  #------ Total Fractional Hand Surface Area ----------
  tfs<-runif(iterations, (0.01), (0.04)) #AuYeung et al 2008 - pinch grip  
  tha<-runif(iterations, 445, 535) #Beamer et. al 2015 & Exposure Factor Handbook
  
  #-----# Calculating Concentration on hand #--------#
  #from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.c.adenovirus<-adenovirus.ts*(tfs)*te.vir.ts2h 
  ts2h.c.norovirus<-norovirus.ts*(tfs)*te.vir.ts2h
  
  #--------------Calculating Dose From Single Touch (non-chain)------------
  #calculating dose from toilet seat touch (Capturing Female to Male use of same toilet)
  ts2h.adenovirus<-((ts2h.c.adenovirus)*te.vir.h2f*tfs*tha) #This is dose need to separate for concentration and then can add the last part for dose 
  ts2h.norovirus<-((ts2h.c.norovirus)*te.vir.h2f*tfs*tha)
  
  #calculating infection risk - TOUCH ADENO
  ts.adenovirus.ir<-1-exp(-6.07E-1*ts2h.adenovirus)
  #TOUCH NORO
  ts.norovirus.ir<-1-((1+(ts2h.norovirus/32.3))^-0.104)
  
  adeno_tsF<-data.frame(adenovirus.ts = adenovirus.ts,
                        adenovirus.ts.gc = adenovirus.ts.gc,
                        ts2h.c.adenovirus = ts2h.c.adenovirus, 
                        ts2h.adenovirus = ts2h.adenovirus, 
                        ts.adenovirus.ir = ts.adenovirus.ir)
  noro_tsF<-data.frame(norovirus.ts = norovirus.ts,
                       norovirus.ts.gc = norovirus.ts.gc,
                       ts2h.c.norovirus = ts2h.c.norovirus, 
                       ts2h.norovirus = ts2h.norovirus, 
                       ts.norovirus.ir = ts.norovirus.ir)
  
  
  write.csv(adeno_tsF,'RAbney_adenomodel.csv', row.names = FALSE)
  write.csv(noro_tsF, 'RAbney_noromodel.csv', row.names = FALSE)
  
  
  #----------Saving Time--------------------
  
  adeno <- read.csv('RAbney_adenomodel.csv', header=T, sep=',')
  adeno.matrix<-st(adeno, digits = 15, add.median = TRUE, out="csv", file= 'R.sum_adeno_model.csv')
  noro <- read.csv('RAbney_noromodel.csv', header=T, sep=',')
  noro.matrix<- st(noro, digits = 15, add.median = TRUE, out="csv", file= 'R.sum_noro_model.csv')
}
Fomite.MindingMatrix.R(10000)


#------------------------------------------------------------------------------------
# STATS for QMRA
#------------------------------------------------------------------------------------
ira<-c(IRA$tsb.adenovirus.ir, IRA$pbs.adenovirus.ir,IRA$fbs.adenovirus.ir,IRA$astm.adenovirus.ir)
matrixtype<-c(rep("tsb",length(IRA$tsb.adenovirus.ir)),
              rep("pbs",length(IRA$pbs.adenovirus.ir)),
              rep("fbs",length(IRA$fbs.adenovirus.ir)),
              rep("astm",length(IRA$astm.adenovirus.ir)))
frame.for.anova<-data.frame(ira,matrixtype)
aovresults<-aov(ira ~ matrixtype,data=frame.for.anova)
summary(aovresults)

#-----Noro-------------------
irn<-c(IRN$tsb.norovirus.ir, IRN$pbs.norovirus.ir,IRN$fbs.norovirus.ir,IRN$astm.norovirus.ir)
matrixtype<-c(rep("tsb",length(IRN$tsb.norovirus.ir)),
              rep("pbs",length(IRN$pbs.norovirus.ir)),
              rep("fbs",length(IRN$fbs.norovirus.ir)),
              rep("astm",length(IRN$astm.norovirus.ir)))
frame.for.anova<-data.frame(irn,matrixtype)
aovresults_noro<-aov(irn ~ matrixtype,data=frame.for.anova)
summary(aovresults_noro)
