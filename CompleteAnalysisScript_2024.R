
## Load package ##

library(metafor)

#########################################################################
############################ PROTECTION #################################
#########################################################################

## Load data ## 

protection <-read.table("HEDGES_protection.txt", header=T)
protectionOR <-read.table("OR_protection.txt", header=T)
protection.s<- read.table("STATISTICS_protection.txt", header=T)


## Calculating effect sizes ##
protectionSMD<-escalc(measure="SMD", m1i= Mean_t, sd1i=SD_t, n1i=N_t,
                      m2i=Mean_c, sd2i=SD_c, n2i=N_c, data=protection)

protectionOR2DN <-escalc(measure="OR2DN", ai=Alive_t, bi=Dead_t, ci=Alive_c,
                         di=Dead_c, data=protectionOR)

## Removing columns we're not using from the table ##
protectionSMD<-protectionSMD[,-c(24,25,26,27,28,29)]
protectionOR2DN<-protectionOR2DN[,-c(24,25,26,27,28,29,30,31)]
protection.s<-protection.s[,-24]

## Merging all effect sizes into the same table ##
protection.data<-merge(protectionSMD,protectionOR2DN,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                                            "Host_strain", "Host_family", "Symbiont_spp", 
                                                            "Symbio_line", "Symbio_type", "Symbio_introgression",
                                                            "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                                                            "Temperature", "Rearing", "NE_infection",
                                                            "Uninfected_symbio_strain", "Experiment",
                                                            "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                       by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                              "Host_strain", "Host_family", "Symbiont_spp", 
                              "Symbio_line", "Symbio_type", "Symbio_introgression",
                              "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                              "Temperature", "Rearing", "NE_infection",
                              "Uninfected_symbio_strain", "Experiment",
                              "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)

protection.data<-merge(protection.data,protection.s,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                                           "Host_strain", "Host_family", "Symbiont_spp", 
                                                           "Symbio_line", "Symbio_type", "Symbio_introgression",
                                                           "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                                                           "Temperature", "Rearing", "NE_infection",
                                                           "Uninfected_symbio_strain", "Experiment",
                                                           "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                       by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                              "Host_strain", "Host_family", "Symbiont_spp", 
                              "Symbio_line", "Symbio_type", "Symbio_introgression",
                              "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                              "Temperature", "Rearing", "NE_infection",
                              "Uninfected_symbio_strain", "Experiment",
                              "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)

protection.data$wi<-1/(protection.data$vi)

## Making a table with effect sizes and moderators we're testing##
protection.data1<- data.frame(Paper_ID=protection.data$Paper_ID, Paper_ES=protection.data$Paper_ES,
                              Symbio_type=protection.data$Symbio_type,Host_family=protection.data$Host_family,
                              Group_measure = protection.data$Group_measure, Symbiont_spp=protection.data$Symbiont_spp, NE_group = protection.data$NE_group,                             
                              yi=protection.data$yi, vi=protection.data$vi, wi=protection.data$wi)

## Removing NAs##
protection.data1<- na.omit(protection.data1)

## Removing moderator categories with less than 10 effect sizes##
protection.data1<-protection.data1[!(protection.data1$Group_measure=="Body_size"),]
protection.data1<-protection.data1[!(protection.data1$Group_measure=="Development_time"),]

## Transforming into factor ##
protection.data1$Paper_ES<-as.factor(protection.data1$Paper_ES)
protection.data1$Paper_ID<-as.factor(protection.data1$Paper_ID)
protection.data1$Host_family<-as.factor(protection.data1$Host_family)
protection.data1$Symbio_type<-as.factor(protection.data1$Symbio_type)
protection.data1$Group_measure<-as.factor(protection.data1$Group_measure)
protection.data1$Symbiont_spp<-as.factor(protection.data1$Symbiont_spp)
protection.data1$NE_group<-as.factor(protection.data1$NE_group)


## Meta-analysis ##


## Overall analysis - study ID and effect size ID as random factors, no moderators included ##
protection.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                           method = "REML", data = protection.data1)
summary(protection.overall)

## Heterogeneity test ##

## Calculating the variance ##
var<-sum(protection.data1$wi*(length(protection.data1$wi)-1))/(sum(protection.data1$wi)^2-sum(protection.data1$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((protection.overall$sigma2[1]+protection.overall$sigma2[2])/
             (protection.overall$sigma2[1]+protection.overall$sigma2[2]+var))*100
I2.total

## Calculating variance between studies ##
Between<-((protection.overall$sigma2[1])/(protection.overall$sigma2[1]+protection.overall$sigma2[2]+var))*100
Between

## Calculating variance within studies ##
Within<-((protection.overall$sigma2[2])/(protection.overall$sigma2[1]+protection.overall$sigma2[2]+var))*100
Within  

## Publication bias - Egger's regression##
egger <- lm(residuals.rma(protection.overall)~protection.data1$wi)
summary(egger)

## Host family as fixed factor, study effect and effect size ID as random factors ##

## Checking the number of effect sizes and studies per moderator category##
table(protection.data1$Host_family, useNA ="always")
pd<-split(protection.data1,protection.data1$Host_family)
table(pd$Aphididae$Paper_ID)
table(pd$Culicidae$Paper_ID)
table(pd$Drosophilidae$Paper_ID)
table(pd$Other$Paper_ID)

## Model ##
protection.gen.fam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, 
                           random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = protection.data1)
summary(protection.gen.fam) 

## Heterogeneity test ##

## Calculating total heterogeneity ##
I2.total.fam<-((protection.gen.fam$sigma2[1]+protection.gen.fam$sigma2[2])/
                 (protection.gen.fam$sigma2[1]+protection.gen.fam$sigma2[2]+var))*100
I2.total.fam

## Calculating variance between studies ##
Between.fam<-((protection.gen.fam$sigma2[1])/(protection.gen.fam$sigma2[1]+protection.gen.fam$sigma2[2]+var))*100
Between.fam

## Calculating variance within studies ##
Within.fam<-((protection.gen.fam$sigma2[2])/(protection.gen.fam$sigma2[1]+protection.gen.fam$sigma2[2]+var))*100
Within.fam  

## Publication bias - Egger's##
egger1 <- lm(residuals.rma(protection.gen.fam)~protection.data1$wi)
summary(egger1)

## Symbiont type as fixed factor, study ID and effect size ID as random factors##

## Checking the number of effect sizes and studies per moderator category ##
table(protection.data1$Symbio_type, useNA ="always")
pd2<-split(protection.data1,protection.data1$Symbio_type)
table(pd2$Natural$Paper_ID)
table(pd2$Not_natural$Paper_ID)

## Model ##
protection.gen.stype<- rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                              random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = protection.data1)
summary(protection.gen.stype)

## Heterogeneity test ##

## Calculating total heterogeneity ##
I2.total.stype<-((protection.gen.stype$sigma2[1]+protection.gen.stype$sigma2[2])/
                   (protection.gen.stype$sigma2[1]+protection.gen.stype$sigma2[2]+var))*100
I2.total.stype

## Calculating variance between studies ##
Between.stype<-((protection.gen.stype$sigma2[1])/(protection.gen.stype$sigma2[1]+protection.gen.stype$sigma2[2]+var))*100
Between.stype

## Calculating variance within studies ##
Within.stype<-((protection.gen.stype$sigma2[2])/(protection.gen.stype$sigma2[1]+protection.gen.stype$sigma2[2]+var))*100
Within.stype  

## Publication bias - Egger's ##
egger2 <- lm(residuals.rma(protection.gen.stype)~protection.data1$wi)
summary(egger2)

##Fitness measure as fixed factor, study ID and study effect size as random factors ##

table(protection.data1$Group_measure, useNA ="always")

pd1<-split(protection.data1,protection.data1$Group_measure)
table(pd1$Fecundity$Paper_ID)
table(pd1$Survival$Paper_ID)
table(pd1$Infection$Paper_ID)
table(pd1$Body_size$Paper_ID)
table(pd1$Development_time$Paper_ID)

## Model ##

protection.gen.fitness<- rma.mv(yi= yi, V =vi, mods = ~Group_measure-1, 
                                random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = protection.data1)
summary(protection.gen.fitness)


## Calculating total heterogeneity ##
I2.total.fitness<-((protection.gen.fitness$sigma2[1]+protection.gen.fitness$sigma2[2])/
                     (protection.gen.fitness$sigma2[1]+protection.gen.fitness$sigma2[2]+var))*100
I2.total.fitness

## Calculating heterogeneity between studies ##
Between.fitness<-((protection.gen.fitness$sigma2[1])/(protection.gen.fitness$sigma2[1]+protection.gen.fitness$sigma2[2]+var))*100
Between.fitness

## Calculating heterogeneity within studies ##
Within.fitness<-((protection.gen.fitness$sigma2[2])/(protection.gen.fitness$sigma2[1]+protection.gen.fitness$sigma2[2]+var))*100
Within.fitness  

## Publication bias -Egger's ##

egger3 <- lm(residuals.rma(protection.gen.fitness)~protection.data1$wi)
summary(egger3)


## Symbiont species as fixed factor, study ID and study effect size as random factors ##

table(protection.data1$Symbiont_spp, useNA ="always")
pds<-split(protection.data1,protection.data1$Symbiont_spp)
table(pds$Wolbachia$Paper_ID)
table(pds$Spiroplasma$Paper_ID)
table(pds$Hamiltonella$Paper_ID)
table(pds$Regiella$Paper_ID)
table(pds$Other$Paper_ID)
table(pds$Hamiltonella_Regiella$Paper_ID)
table(pds$Regiella_Spiroplasma$Paper_ID)

## Model ##

protection.gen.symbio<- rma.mv(yi= yi, V =vi, mods = ~Symbiont_spp-1, 
                               random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = protection.data1)
summary(protection.gen.symbio)

## Calculating total heterogeneity ##
I2.total.symsp<-((protection.gen.symbio$sigma2[1]+protection.gen.symbio$sigma2[2])/
                   (protection.gen.symbio$sigma2[1]+protection.gen.symbio$sigma2[2]+var))*100
I2.total.symsp


## Calculating heterogeneity between studies ##
Between.symsp<-((protection.gen.symbio$sigma2[1])/(protection.gen.symbio$sigma2[1]+protection.gen.symbio$sigma2[2]+var))*100
Between.symsp


## Calculating heterogeneity within studies ##
Within.symsp<-((protection.gen.symbio$sigma2[2])/(protection.gen.symbio$sigma2[1]+protection.gen.symbio$sigma2[2]+var))*100
Within.symsp  

## Publication bias - Eggers's ##
eggersym1 <- lm(residuals.rma(protection.gen.symbio)~protection.data1$wi)
summary(eggersym1)


## Natural enemy group as fixed factor, study ID and study effect size as random factors ##

table(protection.data1$NE_group, useNA ="always")
pdne<-split(protection.data1,protection.data1$NE_group)
table(pdne$Bacteria$Paper_ID)
table(pdne$Fungus$Paper_ID)
table(pdne$Nematode$Paper_ID)
table(pdne$Parasitoid$Paper_ID)
table(pdne$Protozoan$Paper_ID)
table(pdne$Virus$Paper_ID)

## Model ##

protection.gen.ne<- rma.mv(yi= yi, V =vi, mods = ~NE_group-1, 
                           random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = protection.data1)
summary(protection.gen.ne)


## Calculating total heterogeneity ##
I2.total.ne<-((protection.gen.ne$sigma2[1]+protection.gen.ne$sigma2[2])/
                (protection.gen.ne$sigma2[1]+protection.gen.ne$sigma2[2]+var))*100
I2.total.ne

## Calculating heterogeneity between studies ##
Between.ne<-((protection.gen.ne$sigma2[1])/(protection.gen.ne$sigma2[1]+protection.gen.ne$sigma2[2]+var))*100
Between.ne

## Calculating heterogeneity within studies ##
Within.ne<-((protection.gen.ne$sigma2[2])/(protection.gen.ne$sigma2[1]+protection.gen.ne$sigma2[2]+var))*100
Within.ne  

## Publication bias - Egger's ##
eggerne <- lm(residuals.rma(protection.gen.ne)~protection.data1$wi)
summary(eggerne)

#################### Subgroup analyses: Wolbachia only #########################

wolbachia<-protection.data1[protection.data1$Symbiont_spp=="Wolbachia",]

wolbachia$Paper_ES<-as.factor(wolbachia$Paper_ES)
wolbachia$Paper_ID<-as.factor(wolbachia$Paper_ID)
wolbachia$Host_family<-as.factor(wolbachia$Host_family)
wolbachia$Symbio_type<-as.factor(wolbachia$Symbio_type)
wolbachia$Group_measure<-as.factor(wolbachia$Group_measure)
wolbachia$NE_group<-as.factor(wolbachia$NE_group)

wolbachia<-na.omit(wolbachia)

## Removing categories with <10 effect sizes ##

wolbachia<-wolbachia[!(wolbachia$Group_measure=="Fecundity"),]
wolbachia<-wolbachia[!(wolbachia$NE_group=="Fungus"),]
wolbachia<-wolbachia[!(wolbachia$NE_group=="Nematode"),]

## Overall model - no fixed factors ##

pro.w.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                      method = "REML", data = wolbachia)
summary(pro.w.overall)


## Heterogeneity test ##

## Calculating the variance ##
var<-sum(wolbachia$wi*(length(wolbachia$wi)-1))/(sum(wolbachia$wi)^2-sum(wolbachia$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((pro.w.overall$sigma2[1]+pro.w.overall$sigma2[2])/
             (pro.w.overall$sigma2[1]+pro.w.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.w.overall$sigma2[1])/(pro.w.overall$sigma2[1]+pro.w.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.w.overall$sigma2[2])/(pro.w.overall$sigma2[1]+pro.w.overall$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.w1 <- lm(residuals.rma(pro.w.overall)~wolbachia$wi)
summary(egger.w1)

## Host family as fixed factor, study ID and study effect size as random factors ##

table(wolbachia$Host_family, useNA ="always")

pd<-split(wolbachia,wolbachia$Host_family)
table(pd$Culicidae$Paper_ID)
table(pd$Drosophilidae$Paper_ID)
table(pd$Other$Paper_ID)

## Model ##

pro.w.fam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia)
summary(pro.w.fam) 

## Calculating total heterogeneity ##
I2.total<-((pro.w.fam$sigma2[1]+pro.w.fam$sigma2[2])/
             (pro.w.fam$sigma2[1]+pro.w.fam$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.w.fam$sigma2[1])/(pro.w.fam$sigma2[1]+pro.w.fam$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.w.fam$sigma2[2])/(pro.w.fam$sigma2[1]+pro.w.fam$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.w2 <- lm(residuals.rma(pro.w.fam)~wolbachia$wi)
summary(egger.w2)


## Symbiont type as fixed factor, study ID and study effect size as random factors ##

table(wolbachia$Symbio_type, useNA ="always")

pd<-split(wolbachia,wolbachia$Symbio_type)
table(pd$Natural$Paper_ID)
table(pd$Not_natural$Paper_ID)

## Model ##

pro.w.sym<-rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia)
summary(pro.w.sym) 

## Calculating total heterogeneity ##
I2.total<-((pro.w.sym$sigma2[1]+pro.w.sym$sigma2[2])/
             (pro.w.sym$sigma2[1]+pro.w.sym$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.w.sym$sigma2[1])/(pro.w.sym$sigma2[1]+pro.w.sym$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.w.sym$sigma2[2])/(pro.w.sym$sigma2[1]+pro.w.sym$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.w3 <- lm(residuals.rma(pro.w.sym)~wolbachia$wi)
summary(egger.w3)

## Fitness measure as fixed factor and study ID and study effect size as random factors ##

table(wolbachia$Group_measure, useNA ="always")

pd<-split(wolbachia,wolbachia$Group_measure)  
table(pd$Survival$Paper_ID)
table(pd$Infection$Paper_ID)

## Model ##

pro.w.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia)
summary(pro.w.fit) 

## Calculating total heterogeneity ##
I2.total<-((pro.w.fit$sigma2[1]+pro.w.fit$sigma2[2])/
             (pro.w.fit$sigma2[1]+pro.w.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.w.fit$sigma2[1])/(pro.w.fit$sigma2[1]+pro.w.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.w.fit$sigma2[2])/(pro.w.fit$sigma2[1]+pro.w.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.w4 <- lm(residuals.rma(pro.w.fit)~wolbachia$wi)
summary(egger.w4)

##  Natural enemy as fixed factor and study ID and study effect size as random factors ##

table(wolbachia$NE_group, useNA ="always")

pd<-split(wolbachia,wolbachia$NE_group)  
table(pd$Bacteria$Paper_ID)
table(pd$Parasitoid$Paper_ID)
table(pd$Protozoan$Paper_ID)
table(pd$Virus$Paper_ID)

## Model ##

pro.w.ne<-rma.mv(yi = yi, V = vi, mods = ~NE_group-1, 
                 random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia)
summary(pro.w.ne) 

## Calculating total heterogeneity ##
I2.total<-((pro.w.ne$sigma2[1]+pro.w.ne$sigma2[2])/
             (pro.w.ne$sigma2[1]+pro.w.ne$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.w.ne$sigma2[1])/(pro.w.ne$sigma2[1]+pro.w.ne$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.w.ne$sigma2[2])/(pro.w.ne$sigma2[1]+pro.w.ne$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.w5 <- lm(residuals.rma(pro.w.ne)~wolbachia$wi)
summary(egger.w5)


#################### Subgroup analyses: Aphididae only #########################

aphid<-protection.data1[protection.data1$Host_family=="Aphididae",]

aphid$Paper_ES<-as.factor(aphid$Paper_ES)
aphid$Paper_ID<-as.factor(aphid$Paper_ID)
aphid$Host_family<-as.factor(aphid$Host_family)
aphid$Symbio_type<-as.factor(aphid$Symbio_type)
aphid$Group_measure<-as.factor(aphid$Group_measure)
aphid$NE_group<-as.factor(aphid$NE_group)

aphid<-na.omit(aphid)

##Removing categories with <10 effect sizes

table(aphid$Group_measure, useNA ="always")
table(aphid$Symbiont_spp, useNA ="always")
table(aphid$NE_group, useNA ="always")

aphid<-aphid[!(aphid$NE_group=="Bacteria"),]

## Model ##

pro.aph.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                        method = "REML", data = aphid)
summary(pro.aph.overall)



## Calculating variance ##
var<-sum(aphid$wi*(length(aphid$wi)-1))/(sum(aphid$wi)^2-sum(aphid$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((pro.aph.overall$sigma2[1]+pro.aph.overall$sigma2[2])/
             (pro.aph.overall$sigma2[1]+pro.aph.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.aph.overall$sigma2[1])/(pro.aph.overall$sigma2[1]+pro.aph.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.aph.overall$sigma2[2])/(pro.aph.overall$sigma2[1]+pro.aph.overall$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.a1 <- lm(residuals.rma(pro.aph.overall)~aphid$wi)
summary(egger.a1)


## Symbiont type as fixed factor and study ID and study effect size as random factors ##


table(aphid$Symbio_type, useNA ="always")

pd<-split(aphid,aphid$Symbio_type)
table(pd$Natural$Paper_ID)
table(pd$Not_natural$Paper_ID)

## Model ##

pro.aph.sym<-rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                    random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid)
summary(pro.aph.sym) 

## Calculating total heterogeneity ##
I2.total<-((pro.aph.sym$sigma2[1]+pro.aph.sym$sigma2[2])/
             (pro.aph.sym$sigma2[1]+pro.aph.sym$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.aph.sym$sigma2[1])/(pro.aph.sym$sigma2[1]+pro.aph.sym$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.aph.sym$sigma2[2])/(pro.aph.sym$sigma2[1]+pro.aph.sym$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.a2 <- lm(residuals.rma(pro.aph.sym)~aphid$wi)
summary(egger.a2)

## Fitness measure as fixed factor and study ID and study effect size as random factors ##

table(aphid$Group_measure, useNA ="always")

pd<-split(aphid,aphid$Group_measure)  
table(pd$Fecundity$Paper_ID)
table(pd$Infection$Paper_ID)
table(pd$Survival$Paper_ID)

## Model ##

pro.aph.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                    random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid)
summary(pro.aph.fit) 

## Calculating total heterogeneity ##
I2.total<-((pro.aph.fit$sigma2[1]+pro.aph.fit$sigma2[2])/
             (pro.aph.fit$sigma2[1]+pro.aph.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.aph.fit$sigma2[1])/(pro.aph.fit$sigma2[1]+pro.aph.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.aph.fit$sigma2[2])/(pro.aph.fit$sigma2[1]+pro.aph.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.a3 <- lm(residuals.rma(pro.aph.fit)~aphid$wi)
summary(egger.a3)


##  Symbiont species as fixed factor and study ID and study effect size as random factors ##


table(aphid$Symbiont_spp, useNA ="always")

pd<-split(aphid,aphid$Symbiont_spp)  
table(pd$Hamiltonella$Paper_ID)
table(pd$Hamiltonella_Regiella$Paper_ID)
table(pd$Other$Paper_ID)
table(pd$Regiella$Paper_ID)
table(pd$Regiella_Spiroplasma$Paper_ID)
table(pd$Spiroplasma$Paper_ID)

## Model ##

pro.aph.spp<-rma.mv(yi = yi, V = vi, mods = ~Symbiont_spp-1, 
                    random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid)
summary(pro.aph.spp) 


## Calculating total heterogeneity ##
I2.total<-((pro.aph.spp$sigma2[1]+pro.aph.spp$sigma2[2])/
             (pro.aph.spp$sigma2[1]+pro.aph.spp$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.aph.spp$sigma2[1])/(pro.aph.spp$sigma2[1]+pro.aph.spp$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.aph.spp$sigma2[2])/(pro.aph.spp$sigma2[1]+pro.aph.spp$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.a4 <- lm(residuals.rma(pro.aph.spp)~aphid$wi)
summary(egger.a4)


##  Natural enemy as fixed factor and study ID and study effect size as random factors ##

table(aphid$NE_group, useNA ="always")

pd<-split(aphid,aphid$NE_group)  
table(pd$Fungus$Paper_ID)
table(pd$Parasitoid$Paper_ID)

pro.aph.ne<-rma.mv(yi = yi, V = vi, mods = ~NE_group-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid)
summary(pro.aph.ne) 

## Model ##

## Calculating total heterogeneity ##
I2.total<-((pro.aph.ne$sigma2[1]+pro.aph.ne$sigma2[2])/
             (pro.aph.ne$sigma2[1]+pro.aph.ne$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((pro.aph.ne$sigma2[1])/(pro.aph.ne$sigma2[1]+pro.aph.ne$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((pro.aph.ne$sigma2[2])/(pro.aph.ne$sigma2[1]+pro.aph.ne$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.a5 <- lm(residuals.rma(pro.aph.ne)~aphid$wi)
summary(egger.a5)

#####################################################################
########################### COST TO THE HOST ########################
#####################################################################

## load data ##

chost <-read.table("HEDGES_costhost.txt", header=T)
chostOR <-read.table("OR_costhost.txt", header=T)
chost.s<- read.table("STATISTICS_costhost.txt", header=T)

# calculating all fitness measures together, then delete development time ##
chostSMD<-escalc(measure="SMD", m1i= Mean_t, sd1i=SD_t, n1i=N_t,
                 m2i=Mean_c, sd2i=SD_c, n2i=N_c, data=chost)
chostSMD<-chostSMD[!(chostSMD$Group_measure=="Development_time"),]

## calculating development time separately so positive values in the graphs represents decreased development time ##
Devtime<-chost[chost$Group_measure=="Development_time",]

DevtimeSMD<-escalc(measure="SMD", m1i= Mean_c, sd1i=SD_c, n1i=N_c,
                   m2i=Mean_t, sd2i=SD_t, n2i=N_t, data=Devtime)

## calculating odds ratio into SMD ##
chostOR2DN <-escalc(measure="OR2DN", ai=Alive_t, bi=Dead_t, ci=Alive_c,
                    di=Dead_c, data=chostOR)

## deleting columns with data we won't use in the analyses ##
chostSMD<-chostSMD[,-c(21,22,23,24,25,26)]
DevtimeSMD<-DevtimeSMD[,-c(21,22,23,24,25,26)]
chostOR2DN<-chostOR2DN[,-c(21,22,23,24,25,26,27,28)]
chost.s<-chost.s[,-21]

## merging all effect sizes calculated into the same table ##
chost.data<-merge(chostSMD,chostOR2DN,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                             "Host_strain", "Host_family", "Symbiont_spp", 
                                             "Symbio_line", "Symbio_type", "Symbio_introgression",
                                             "Temperature", "Rearing", "Nutrition", "Competition",
                                             "Uninfected_symbio_strain", "Experiment",
                                             "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                  by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                         "Host_strain", "Host_family", "Symbiont_spp", 
                         "Symbio_line", "Symbio_type", "Symbio_introgression",
                         "Temperature", "Rearing", "Nutrition", "Competition",
                         "Uninfected_symbio_strain", "Experiment",
                         "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)

chost.data<-merge(chost.data,chost.s,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                            "Host_strain", "Host_family", "Symbiont_spp", 
                                            "Symbio_line", "Symbio_type", "Symbio_introgression",
                                            "Temperature", "Rearing", "Nutrition", "Competition",
                                            "Uninfected_symbio_strain", "Experiment",
                                            "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                  by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                         "Host_strain", "Host_family", "Symbiont_spp", 
                         "Symbio_line", "Symbio_type", "Symbio_introgression",
                         "Temperature", "Rearing", "Nutrition", "Competition",
                         "Uninfected_symbio_strain", "Experiment",
                         "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)
chost.data<-merge(chost.data,DevtimeSMD,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                               "Host_strain", "Host_family", "Symbiont_spp", 
                                               "Symbio_line", "Symbio_type", "Symbio_introgression",
                                               "Temperature", "Rearing", "Nutrition", "Competition",
                                               "Uninfected_symbio_strain", "Experiment",
                                               "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                  by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                         "Host_strain", "Host_family", "Symbiont_spp", 
                         "Symbio_line", "Symbio_type", "Symbio_introgression",
                         "Temperature", "Rearing", "Nutrition", "Competition",
                         "Uninfected_symbio_strain", "Experiment",
                         "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)

chost.data$wi<-1/(chost.data$vi)

chost.data1<- data.frame(Paper_ID=chost.data$Paper_ID, Paper_ES=chost.data$Paper_ES,
                         Symbio_type=chost.data$Symbio_type,Host_family=chost.data$Host_family,Symbiont_spp=chost.data$Symbiont_spp,
                         Group_measure = chost.data$Group_measure, yi=chost.data$yi,vi=chost.data$vi, wi=chost.data$wi)

## deleting NAs ##
chost.data1<- na.omit(chost.data1)

## transforming into factor ##
chost.data1$Paper_ES<-as.factor(chost.data1$Paper_ES)
chost.data1$Paper_ID<-as.factor(chost.data1$Paper_ID)
chost.data1$Host_family<-as.factor(chost.data1$Host_family)
chost.data1$Symbio_type<-as.factor(chost.data1$Symbio_type)
chost.data1$Group_measure<-as.factor(chost.data1$Group_measure)
chost.data1$Symbiont_spp<-as.factor(chost.data1$Symbiont_spp)

## Overall - no fixed factors, study ID and effect size ID as random factors ##

chost.general<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                      method = "REML", data = chost.data1)
summary(chost.general)

## Calculating variance ##
varch<-sum(chost.data1$wi*(length(chost.data1$wi)-1))/(sum(chost.data1$wi)^2-sum(chost.data1$wi^2))

## Calculating total heterogeneity ##
I2.total.ch<-((chost.general$sigma2[1]+chost.general$sigma2[2])/
                (chost.general$sigma2[1]+chost.general$sigma2[2]+varch))*100
I2.total.ch

## Calculating heterogeneity between studies ##
Between.ch<-((chost.general$sigma2[1])/(chost.general$sigma2[1]+chost.general$sigma2[2]+varch))*100
Between.ch

## Calculating heterogeneity within studies ##
Within.ch<-((chost.general$sigma2[2])/(chost.general$sigma2[1]+chost.general$sigma2[2]+varch))*100
Within.ch  

## Publication bias ##
egger4 <- lm(residuals.rma(chost.general)~chost.data1$wi)
summary(egger4)

## Host family as fixed factor, study ID and study effect size as random factors ##

table(chost.data1$Host_family, useNA ="always")
pd3<-split(chost.data1,chost.data1$Host_family)
table(pd3$Aphididae$Paper_ID)
table(pd3$Culicidae$Paper_ID)
table(pd3$Drosophilidae$Paper_ID)
table(pd3$Other$Paper_ID)

## Model ##

chost.gen.fam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, 
                      random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = chost.data1)
summary(chost.gen.fam) 


## Calculating total heterogeneity ##
I2.total.ch.fam<-((chost.gen.fam$sigma2[1]+chost.gen.fam$sigma2[2])/
                    (chost.gen.fam$sigma2[1]+chost.gen.fam$sigma2[2]+varch))*100
I2.total.ch.fam

## Calculating heterogeneity between tudies ##
Between.ch.fam<-((chost.gen.fam$sigma2[1])/(chost.gen.fam$sigma2[1]+chost.gen.fam$sigma2[2]+varch))*100
Between.ch.fam

## Calculating heterogeneity within studies ##
Within.ch.fam<-((chost.gen.fam$sigma2[2])/(chost.gen.fam$sigma2[1]+chost.gen.fam$sigma2[2]+varch))*100
Within.ch.fam  

## Publication bias ##
egger5 <- lm(residuals.rma(chost.gen.fam)~chost.data1$wi)
summary(egger5)

### Symbiont type as fixed factor, study ID and effect size ID as random factors ##

table(chost.data1$Symbio_type, useNA ="always")
pd4<-split(chost.data1,chost.data1$Symbio_type)
table(pd4$Not_natural$Paper_ID)
table(pd4$Natural$Paper_ID)

## Model ##

chost.gen.stype<- rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                         random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = chost.data1)
summary(chost.gen.stype)


## Calculating total heterogeneity ##
I2.total.ch.stype<-((chost.gen.stype$sigma2[1]+chost.gen.stype$sigma2[2])/
                      (chost.gen.stype$sigma2[1]+chost.gen.stype$sigma2[2]+varch))*100
I2.total.ch.stype

## Calculating heterogeneity between studies ##
Between.ch.stype<-((chost.gen.stype$sigma2[1])/(chost.gen.stype$sigma2[1]+chost.gen.stype$sigma2[2]+varch))*100
Between.ch.stype

## Calculating heterogeneity within studies ##
Within.ch.stype<-((chost.gen.stype$sigma2[2])/(chost.gen.stype$sigma2[1]+chost.gen.stype$sigma2[2]+varch))*100
Within.ch.stype  

## Publication bias ##
egger6 <- lm(residuals.rma(chost.gen.stype)~chost.data1$wi)
summary(egger6)

## Fitness measure as fixed factor, study ID and effect size ID as random factors ##

table(chost.data1$Group_measure, useNA ="always")
pd5<-split(chost.data1,chost.data1$Group_measure)
table(pd5$Fecundity$Paper_ID)
table(pd5$Survival$Paper_ID)
table(pd5$Body_size$Paper_ID)
table(pd5$Development_time$Paper_ID)

## Model ##

chost.gen.fitness<- rma.mv(yi= yi, V =vi, mods = ~Group_measure-1, 
                           random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = chost.data1)
summary(chost.gen.fitness)


## Calculating total heterogeneity ##
I2.total.ch.fitness<-((chost.gen.fitness$sigma2[1]+chost.gen.fitness$sigma2[2])/
                        (chost.gen.fitness$sigma2[1]+chost.gen.fitness$sigma2[2]+varch))*100
I2.total.ch.fitness

## Calculating heterogeneity between studies ##
Between.ch.fitness<-((chost.gen.fitness$sigma2[1])/(chost.gen.fitness$sigma2[1]+chost.gen.fitness$sigma2[2]+varch))*100
Between.ch.fitness

## Calculating heterogeneity within studies ##
Within.ch.fitness<-((chost.gen.fitness$sigma2[2])/(chost.gen.fitness$sigma2[1]+chost.gen.fitness$sigma2[2]+varch))*100
Within.ch.fitness  

## Publication bias ##
egger7<- lm(residuals.rma(chost.gen.fitness)~chost.data1$wi)
summary(egger7)

## Symbiont species as fixed factor, study ID and effect size ID as random factors ##

table(chost.data1$Symbiont_spp, useNA ="always")
pds1<-split(chost.data1,chost.data1$Symbiont_spp)
table(pds1$Cardinium_Wolbachia$Paper_ID)
table(pds1$Hamiltonella$Paper_ID)
table(pds1$Hamiltonella_Rickettsia$Paper_ID)
table(pds1$Hamiltonella_Rickettsiella$Paper_ID)
table(pds1$Regiella$Paper_ID)
table(pds1$Regiella_Spiroplasma$Paper_ID)
table(pds1$Rickettsia$Paper_ID)
table(pds1$Rickettsiella$Paper_ID)
table(pds1$Serratia$Paper_ID)
table(pds1$Spiroplasma$Paper_ID)
table(pds1$Spiroplasma_Wolbachia$Paper_ID)
table(pds1$Wolbachia$Paper_ID)
table(pds1$Other$Paper_ID)

## Model ##

chost.symbiosp<- rma.mv(yi= yi, V =vi, mods = ~Symbiont_spp-1, 
                        random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = chost.data1)
summary(chost.symbiosp)

#heterogeneity
I2.total.ssp<-((chost.symbiosp$sigma2[1]+chost.symbiosp$sigma2[2])/
                 (chost.symbiosp$sigma2[1]+chost.symbiosp$sigma2[2]+varch))*100
I2.total.ssp

## Calculating heterogeneity between studies ##
Between.ssp<-((chost.symbiosp$sigma2[1])/(chost.symbiosp$sigma2[1]+chost.symbiosp$sigma2[2]+varch))*100
Between.ssp

## Calculating heterogeneity within studies ##
Within.ssp<-((chost.symbiosp$sigma2[2])/(chost.symbiosp$sigma2[1]+chost.symbiosp$sigma2[2]+varch))*100
Within.ssp  

## Publication bias ##
eggerssymspch <- lm(residuals.rma(chost.symbiosp)~chost.data1$wi)
summary(eggerssymspch)


#################### Subgroup analyses: Wolbachia only #########################

wolbachia2<-chost.data1[chost.data1$Symbiont_spp=="Wolbachia",]

wolbachia2$Paper_ES<-as.factor(wolbachia2$Paper_ES)
wolbachia2$Paper_ID<-as.factor(wolbachia2$Paper_ID)
wolbachia2$Host_family<-as.factor(wolbachia2$Host_family)
wolbachia2$Symbio_type<-as.factor(wolbachia2$Symbio_type)
wolbachia2$Group_measure<-as.factor(wolbachia2$Group_measure)

wolbachia2<-na.omit(wolbachia2)

## Overall model - no moderators included, study ID and effect size ID as random factors##

cost.w.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                       method = "REML", data = wolbachia2)
summary(cost.w.overall)


## Calculating the variance ##
var<-sum(wolbachia2$wi*(length(wolbachia2$wi)-1))/(sum(wolbachia2$wi)^2-sum(wolbachia2$wi^2))
var


## Calculating total heterogeneity ##
I2.total<-((cost.w.overall$sigma2[1]+cost.w.overall$sigma2[2])/
             (cost.w.overall$sigma2[1]+cost.w.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cost.w.overall$sigma2[1])/(cost.w.overall$sigma2[1]+cost.w.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cost.w.overall$sigma2[2])/(cost.w.overall$sigma2[1]+cost.w.overall$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.chw1 <- lm(residuals.rma(cost.w.overall)~wolbachia2$wi)
summary(egger.chw1)

## Host family as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia2$Host_family, useNA ="always")
pd<-split(wolbachia2,wolbachia2$Host_family)
table(pd$Culicidae$Paper_ID)
table(pd$Drosophilidae$Paper_ID)
table(pd$Other$Paper_ID)

## Model ##

cost.w.fam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia2)
summary(cost.w.fam) 

## Calculating total heterogeneity ##
I2.total<-((cost.w.fam$sigma2[1]+cost.w.fam$sigma2[2])/
             (cost.w.fam$sigma2[1]+cost.w.fam$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cost.w.fam$sigma2[1])/(cost.w.fam$sigma2[1]+cost.w.fam$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cost.w.fam$sigma2[2])/(cost.w.fam$sigma2[1]+cost.w.fam$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.chw2 <- lm(residuals.rma(cost.w.fam)~wolbachia2$wi)
summary(egger.chw2)

## Symbiont type as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia2$Symbio_type, useNA ="always")
pd<-split(wolbachia2,wolbachia2$Symbio_type)
table(pd$Natural$Paper_ID)
table(pd$Not_natural$Paper_ID)

## Model ##

cost.w.sym<-rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia2)
summary(cost.w.sym) 

## Calculating total heterogeneity ##
I2.total<-((cost.w.sym$sigma2[1]+cost.w.sym$sigma2[2])/
             (cost.w.sym$sigma2[1]+cost.w.sym$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cost.w.sym$sigma2[1])/(cost.w.sym$sigma2[1]+cost.w.sym$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cost.w.sym$sigma2[2])/(cost.w.sym$sigma2[1]+cost.w.sym$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.chw3 <- lm(residuals.rma(cost.w.sym)~wolbachia2$wi)
summary(egger.chw3)

## Fitness measure as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia2$Group_measure, useNA ="always")

pd<-split(wolbachia2,wolbachia2$Group_measure)  
table(pd$Body_size$Paper_ID)
table(pd$Development_time$Paper_ID)
table(pd$Survival$Paper_ID)
table(pd$Fecundity$Paper_ID)

## Model ##

cost.w.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia2)
summary(cost.w.fit) 

## Calculating total heterogeneity ##
I2.total<-((cost.w.fit$sigma2[1]+cost.w.fit$sigma2[2])/
             (cost.w.fit$sigma2[1]+cost.w.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cost.w.fit$sigma2[1])/(cost.w.fit$sigma2[1]+cost.w.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cost.w.fit$sigma2[2])/(cost.w.fit$sigma2[1]+cost.w.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.chw4 <- lm(residuals.rma(cost.w.fit)~wolbachia2$wi)
summary(egger.chw4)

#################### Subgroup analyses: Aphididae only #########################

aphid2<-chost.data1[chost.data1$Host_family=="Aphididae",]

aphid2$Paper_ES<-as.factor(aphid2$Paper_ES)
aphid2$Paper_ID<-as.factor(aphid2$Paper_ID)
aphid2$Host_family<-as.factor(aphid2$Host_family)
aphid2$Symbio_type<-as.factor(aphid2$Symbio_type)
aphid2$Group_measure<-as.factor(aphid2$Group_measure)

aphid2<-na.omit(aphid2)

## Overall model - no moderators included,  study ID and effect size ID as random factors ##

ch.aph.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                       method = "REML", data = aphid2)
summary(ch.aph.overall)

## Calculating the variance ##
var<-sum(aphid2$wi*(length(aphid2$wi)-1))/(sum(aphid2$wi)^2-sum(aphid2$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((ch.aph.overall$sigma2[1]+ch.aph.overall$sigma2[2])/
             (ch.aph.overall$sigma2[1]+ch.aph.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ch.aph.overall$sigma2[1])/(ch.aph.overall$sigma2[1]+ch.aph.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ch.aph.overall$sigma2[2])/(ch.aph.overall$sigma2[1]+ch.aph.overall$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cha1 <- lm(residuals.rma(ch.aph.overall)~aphid2$wi)
summary(egger.cha1)

## Symbiont type as fixed factor, study ID and effect size ID as random factors ##

table(aphid2$Symbio_type, useNA ="always")
pd<-split(aphid2,aphid2$Symbio_type)
table(pd$Natural$Paper_ID)
table(pd$Not_natural$Paper_ID)

## Model ##

ch.aph.sym<-rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid2)
summary(ch.aph.sym) 

## Calculating total heterogeneity ##
I2.total<-((ch.aph.sym$sigma2[1]+ch.aph.sym$sigma2[2])/
             (ch.aph.sym$sigma2[1]+ch.aph.sym$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ch.aph.sym$sigma2[1])/(ch.aph.sym$sigma2[1]+ch.aph.sym$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ch.aph.sym$sigma2[2])/(ch.aph.sym$sigma2[1]+ch.aph.sym$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cha2 <- lm(residuals.rma(ch.aph.sym)~aphid2$wi)
summary(egger.cha2)

## Fitness measure as fixed factor, study ID and effect size ID as random factors ##

table(aphid2$Group_measure, useNA ="always")
pd<-split(aphid2,aphid2$Group_measure)  
table(pd$Fecundity$Paper_ID)
table(pd$Body_size$Paper_ID)
table(pd$Development_time$Paper_ID)
table(pd$Survival$Paper_ID)

## Model ##

ch.aph.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid2)
summary(ch.aph.fit) 

## Calculating total heterogeneity ##
I2.total<-((ch.aph.fit$sigma2[1]+ch.aph.fit$sigma2[2])/
             (ch.aph.fit$sigma2[1]+ch.aph.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ch.aph.fit$sigma2[1])/(ch.aph.fit$sigma2[1]+ch.aph.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ch.aph.fit$sigma2[2])/(ch.aph.fit$sigma2[1]+ch.aph.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cha3 <- lm(residuals.rma(ch.aph.fit)~aphid2$wi)
summary(egger.cha3)

##  Symbiont species as fixed factor, study ID and effect size ID as random factors ##

table(aphid2$Symbiont_spp, useNA ="always")

pd<-split(aphid2,aphid2$Symbiont_spp)  
table(pd$Hamiltonella$Paper_ID)
table(pd$Hamiltonella_Rickettsia$Paper_ID)
table(pd$Hamiltonella_Rickettsiella$Paper_ID)
table(pd$Other$Paper_ID)
table(pd$Regiella$Paper_ID)
table(pd$Regiella_Spiroplasma$Paper_ID)
table(pd$Rickettsia$Paper_ID)
table(pd$Rickettsiella$Paper_ID)
table(pd$Serratia$Paper_ID)
table(pd$Spiroplasma$Paper_ID)

## Model ##

ch.aph.spp<-rma.mv(yi = yi, V = vi, mods = ~Symbiont_spp-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid2)
summary(ch.aph.spp) 


## Calculating total heterogeneity ##
I2.total<-((ch.aph.spp$sigma2[1]+ch.aph.spp$sigma2[2])/
             (ch.aph.spp$sigma2[1]+ch.aph.spp$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ch.aph.spp$sigma2[1])/(ch.aph.spp$sigma2[1]+ch.aph.spp$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ch.aph.spp$sigma2[2])/(ch.aph.spp$sigma2[1]+ch.aph.spp$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cha4 <- lm(residuals.rma(ch.aph.spp)~aphid2$wi)
summary(egger.cha4)

##########################################################################
########################## COST TO NATURAL ENEMIES #######################
##########################################################################

## Load data ##

cne <-read.table("HEDGES_naturalenemy.txt", header=T)
cneOR <-read.table("OR_naturalenemy.txt", header=T)
cne.s<- read.table("STATISTICS_naturalenemy.txt", header=T)

## Calculating effect sizes ##

cneSMD<-escalc(measure="SMD", m1i= Mean_t, sd1i=SD_t, n1i=N_t,
               m2i=Mean_c, sd2i=SD_c, n2i=N_c, data=cne)

cneOR2DN <-escalc(measure="OR2DN", ai=Alive_t, bi=Dead_t, ci=Alive_c,
                  di=Dead_c, data=cneOR)

## Removing columns we're not using ##
cneSMD<-cneSMD[,-c(24,25,26,27,28,29)]
cneOR2DN<-cneOR2DN[,-c(24,25,26,27,28,29,30,31)]
cne.s<-cne.s[,-24]

## Merging all effect sizes into the same table ##
cne.data<-merge(cneSMD,cneOR2DN,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                       "Host_strain", "Host_family", "Symbiont_spp", 
                                       "Symbio_line", "Symbio_type", "Symbio_introgression",
                                       "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                                       "Temperature", "Rearing", "NE_infection",
                                       "Uninfected_symbio_strain", "Experiment",
                                       "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                       "Host_strain", "Host_family", "Symbiont_spp", 
                       "Symbio_line", "Symbio_type", "Symbio_introgression",
                       "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                       "Temperature", "Rearing", "NE_infection",
                       "Uninfected_symbio_strain", "Experiment",
                       "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)
cne.data<-merge(cne.data,cne.s,by.x=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                                      "Host_strain", "Host_family", "Symbiont_spp", 
                                      "Symbio_line", "Symbio_type", "Symbio_introgression",
                                      "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                                      "Temperature", "Rearing", "NE_infection",
                                      "Uninfected_symbio_strain", "Experiment",
                                      "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), 
                by.y=c("Paper_ES", "Paper_ID", "Paper", "Host_spp", 
                       "Host_strain", "Host_family", "Symbiont_spp", 
                       "Symbio_line", "Symbio_type", "Symbio_introgression",
                       "Natural_enemy_spp", "NE_group", "NE_line", "NE_type", 
                       "Temperature", "Rearing", "NE_infection",
                       "Uninfected_symbio_strain", "Experiment",
                       "Sex", "Stage", "Measure", "Group_measure", "Obs", "yi", "vi"), all.x=TRUE, all.y=TRUE)

cne.data$wi<-1/(cne.data$vi)

cne.data1<- data.frame(Paper_ID=cne.data$Paper_ID, Paper_ES=cne.data$Paper_ES,
                       Host_family=cne.data$Host_family,Symbio_type=cne.data$Symbio_type,
                       Group_measure = cne.data$Group_measure, Symbiont_spp=cne.data$Symbiont_spp, NE_group = cne.data$NE_group,
                       yi=cne.data$yi, vi=cne.data$vi, wi=cne.data$wi)
## Removing NAs ##
cne.data1<- na.omit(cne.data1)

## Removing moderator levels with less than 10 effect sizes ##
cne.data1<-cne.data1[!(cne.data1$Group_measure=="Development_time"),]
cne.data1<-cne.data1[!(cne.data1$NE_group=="Predator"),]

## Transforming into factors ##
cne.data1$Paper_ES<-as.factor(cne.data1$Paper_ES)
cne.data1$Paper_ID<-as.factor(cne.data1$Paper_ID)
cne.data1$Host_family<-as.factor(cne.data1$Host_family)
cne.data1$Symbio_type<-as.factor(cne.data1$Symbio_type)
cne.data1$Group_measure<-as.factor(cne.data1$Group_measure)
cne.data1$Symbiont_spp<-as.factor(cne.data1$Symbiont_spp)
cne.data1$NE_group<-as.factor(cne.data1$NE_group)


## Overall model - no moderators included, effect size ID and study ID as random factors ##

cne.general<-rma.mv(yi = yi, V = vi, random = list(~1 | Paper_ID, ~1 | Paper_ES), 
                    method = "REML", data = cne.data1)
summary(cne.general)


## Heterogeneity test ##

cne.data1$wi<-1/(cne.data1$vi)

## calculating the variance ##
varcne<-sum(cne.data1$wi*(length(cne.data1$wi)-1))/(sum(cne.data1$wi)^2-sum(cne.data1$wi^2))

## Calculating total heterogeneity ##
I2.total.cne<-((cne.general$sigma2[1]+cne.general$sigma2[2])/
                 (cne.general$sigma2[1]+cne.general$sigma2[2]+varcne))*100
I2.total.cne

## Calculating heterogeneity between studies ##
Between.cne<-((cne.general$sigma2[1])/(cne.general$sigma2[1]+cne.general$sigma2[2]+varcne))*100
Between.cne

## Calculating heterogeneity within studies ##
Within.cne<-((cne.general$sigma2[2])/(cne.general$sigma2[1]+cne.general$sigma2[2]+varcne))*100
Within.cne

## Publication bias ##

egger10 <- lm(residuals.rma(cne.general)~cne.data1$wi)
summary(egger10)

## Host family as fixed factor, study ID and effect size ID as random factor ##

table(cne.data1$Host_family, useNA ="always")
pd6<-split(cne.data1,cne.data1$Host_family)
table(pd6$Aphididae$Paper_ID)
table(pd6$Drosophilidae$Paper_ID)
table(pd6$Other$Paper_ID)
table(pd6$Culicidae$Paper_ID)

## Model ##

cne.hfam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, random = list(~1 | Paper_ID, ~1 | Paper_ES), 
                 method = "REML", data = cne.data1)
summary(cne.hfam)


## Calculating total heterogeneity ##
I2.total<-((cne.hfam$sigma2[1]+cne.hfam$sigma2[2])/
             (cne.hfam$sigma2[1]+cne.hfam$sigma2[2]+varcne))*100
I2.total


## Calculating heterogeneity between studies ##
Between<-((cne.hfam$sigma2[1])/(cne.hfam$sigma2[1]+cne.hfam$sigma2[2]+varcne))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.hfam$sigma2[2])/(cne.hfam$sigma2[1]+cne.hfam$sigma2[2]+varcne))*100
Within  

## Publication bias ##
egger12 <- lm(residuals.rma(cne.hfam)~cne.data1$wi)
summary(egger12)

## Fitness measure as fixed factor, study ID and effect size ID as random factor ##

table(cne.data1$Group_measure, useNA ="always")
pd7<-split(cne.data1,cne.data1$Group_measure)
table(pd7$Survival$Paper_ID)
table(pd7$Body_size$Paper_ID)
table(pd7$NE_load$Paper_ID)

## Model ##

cne.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, random = list(~1 | Paper_ID, ~1 | Paper_ES), 
                method = "REML", data = cne.data1)
summary(cne.fit)

## Calculating total heterogeneity ##
I2.total<-((cne.fit$sigma2[1]+cne.fit$sigma2[2])/
             (cne.fit$sigma2[1]+cne.fit$sigma2[2]+varcne))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cne.fit$sigma2[1])/(cne.fit$sigma2[1]+cne.fit$sigma2[2]+varcne))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.fit$sigma2[2])/(cne.fit$sigma2[1]+cne.fit$sigma2[2]+varcne))*100
Within  

## Publication bias ##
egger13 <- lm(residuals.rma(cne.fit)~cne.data1$wi)
summary(egger13)

## Symbiont type as fixed factor, study ID and effect size ID as random factor ##

table(cne.data1$Symbio_type, useNA ="always")
pd8<-split(cne.data1,cne.data1$Symbio_type)
table(pd8$Natural$Paper_ID)
table(pd8$Not_natural$Paper_ID)


## Model ##

cne.stype<- rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = cne.data1)
summary(cne.stype)


## Calculating total heterogeneity ##
I2.total.cne.stype<-((cne.stype$sigma2[1]+cne.stype$sigma2[2])/
                       (cne.stype$sigma2[1]+cne.stype$sigma2[2]+varcne))*100
I2.total.cne.stype

## Calculating heterogeneity between studies ##
Between.cne.stype<-((cne.stype$sigma2[1])/(cne.stype$sigma2[1]+cne.stype$sigma2[2]+varcne))*100
Between.cne.stype


## Calculating heterogeneity within studies ##
Within.cne.stype<-((cne.stype$sigma2[2])/(cne.stype$sigma2[1]+cne.stype$sigma2[2]+varcne))*100
Within.cne.stype  

## Publication bias ##
egger14 <- lm(residuals.rma(cne.stype)~cne.data1$wi)
summary(egger14)


## Symbiont species as fixed factor, study ID and effect size ID as random factor ##

table(cne.data1$Symbiont_spp, useNA ="always")
pds2<-split(cne.data1,cne.data1$Symbiont_spp)
table(pds2$Hamiltonella$Paper_ID)
table(pds2$Other$Paper_ID)
table(pds2$Spiroplasma$Paper_ID)
table(pds2$Wolbachia$Paper_ID)


## Model ##

cne.symbiosp<- rma.mv(yi= yi, V =vi, mods = ~Symbiont_spp-1, 
                      random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = cne.data1)
summary(cne.symbiosp)


## Calculating total heterogeneity ##
I2.total.cne.symbiosp<-((cne.symbiosp$sigma2[1]+cne.symbiosp$sigma2[2])/
                          (cne.symbiosp$sigma2[1]+cne.symbiosp$sigma2[2]+varcne))*100
I2.total.cne.symbiosp

## Calculating heterogeneity between studies ##
Between.cne.symbiosp<-((cne.symbiosp$sigma2[1])/(cne.symbiosp$sigma2[1]+cne.symbiosp$sigma2[2]+varcne))*100
Between.cne.symbiosp

## Calculating heterogeneity within studies ##
Within.cne.symbiosp<-((cne.symbiosp$sigma2[2])/(cne.symbiosp$sigma2[1]+cne.symbiosp$sigma2[2]+varcne))*100
Within.cne.symbiosp  

## Publication bias ##
eggercnesymbiosp <- lm(residuals.rma(cne.symbiosp)~cne.data1$wi)
summary(eggercnesymbiosp)

## Natural enemy group as fixed factor, study ID and effect size ID as random factor ##

table(cne.data1$NE_group, useNA ="always")
pdsne2<-split(cne.data1,cne.data1$NE_group)
table(pdsne2$Bacteria$Paper_ID)
table(pdsne2$Nematode$Paper_ID)
table(pdsne2$Parasitoid$Paper_ID)
table(pdsne2$Protozoan$Paper_ID)
table(pdsne2$Virus$Paper_ID)


## Model ##

cne.ne<- rma.mv(yi= yi, V =vi, mods = ~NE_group-1, 
                random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = cne.data1)
summary(cne.ne)

## Calculating total heterogeneity ##
I2.total.cne.ne<-((cne.ne$sigma2[1]+cne.ne$sigma2[2])/
                    (cne.ne$sigma2[1]+cne.ne$sigma2[2]+varcne))*100
I2.total.cne.ne

## Calculating heterogeneity between studies ##
Between.cne.ne<-((cne.ne$sigma2[1])/(cne.ne$sigma2[1]+cne.ne$sigma2[2]+varcne))*100
Between.cne.ne

## Calculating heterogeneity within studies ##
Within.cne.ne<-((cne.ne$sigma2[2])/(cne.ne$sigma2[1]+cne.ne$sigma2[2]+varcne))*100
Within.cne.ne  

## Publication bias ##
eggercnene <- lm(residuals.rma(cne.ne)~cne.data1$wi)
summary(eggercnene)

########################## Subgroup analyses: Wolbachia only ##############################

wolbachia3<-cne.data1[cne.data1$Symbiont_spp=="Wolbachia",]

wolbachia3$Paper_ES<-as.factor(wolbachia3$Paper_ES)
wolbachia3$Paper_ID<-as.factor(wolbachia3$Paper_ID)
wolbachia3$Host_family<-as.factor(wolbachia3$Host_family)
wolbachia3$Symbio_type<-as.factor(wolbachia3$Symbio_type)
wolbachia3$Group_measure<-as.factor(wolbachia3$Group_measure)
wolbachia3$NE_group<-as.factor(wolbachia3$NE_group)

wolbachia3<-na.omit(wolbachia3)

##Removing categories with <10 effect sizes
wolbachia3<-wolbachia3[!(wolbachia3$Host_family=="Other"),]

## Overall model - no moderators included, study ID and effect size ID as random factors ##

cne.w.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                      method = "REML", data = wolbachia3)
summary(cne.w.overall)

## Calculating the variance ##
var<-sum(wolbachia3$wi*(length(wolbachia3$wi)-1))/(sum(wolbachia3$wi)^2-sum(wolbachia3$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((cne.w.overall$sigma2[1]+cne.w.overall$sigma2[2])/
             (cne.w.overall$sigma2[1]+cne.w.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cne.w.overall$sigma2[1])/(cne.w.overall$sigma2[1]+cne.w.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.w.overall$sigma2[2])/(cne.w.overall$sigma2[1]+cne.w.overall$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnew1 <- lm(residuals.rma(cne.w.overall)~wolbachia3$wi)
summary(egger.cnew1)

## Host family as fixed factor and study ID and study effect size as random factors ##

table(wolbachia3$Host_family, useNA ="always")
pd<-split(wolbachia3,wolbachia3$Host_family)
table(pd$Culicidae$Paper_ID)
table(pd$Drosophilidae$Paper_ID)

## Model ##

cne.w.fam<-rma.mv(yi = yi, V = vi, mods = ~Host_family-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia3)
summary(cne.w.fam) 


## Calculating total heterogeneity ##
I2.total<-((cne.w.fam$sigma2[1]+cne.w.fam$sigma2[2])/
             (cne.w.fam$sigma2[1]+cne.w.fam$sigma2[2]+var))*100
I2.total

## Calculating total between heterogeneity ##
Between<-((cne.w.fam$sigma2[1])/(cne.w.fam$sigma2[1]+cne.w.fam$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.w.fam$sigma2[2])/(cne.w.fam$sigma2[1]+cne.w.fam$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnew2 <- lm(residuals.rma(cne.w.fam)~wolbachia3$wi)
summary(egger.cnew2)

## Symbiont type as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia3$Symbio_type, useNA ="always")

pd<-split(wolbachia3,wolbachia3$Symbio_type)
table(pd$Natural$Paper_ID)
table(pd$Not_natural$Paper_ID)

## Model ##

cne.w.sym<-rma.mv(yi = yi, V = vi, mods = ~Symbio_type-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia3)
summary(cne.w.sym) 


## Calculating total heterogeneity ##
I2.total<-((cne.w.sym$sigma2[1]+cne.w.sym$sigma2[2])/
             (cne.w.sym$sigma2[1]+cne.w.sym$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cne.w.sym$sigma2[1])/(cne.w.sym$sigma2[1]+cne.w.sym$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.w.sym$sigma2[2])/(cne.w.sym$sigma2[1]+cne.w.sym$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnew3 <- lm(residuals.rma(cne.w.sym)~wolbachia3$wi)
summary(egger.cnew3)

## Fitness measure as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia3$Group_measure, useNA ="always")

pd<-split(wolbachia3,wolbachia3$Group_measure)  
table(pd$Survival$Paper_ID)
table(pd$NE_load$Paper_ID)

## Model ##

cne.w.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                  random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia3)
summary(cne.w.fit) 


## Calculating total heterogeneity ##
I2.total<-((cne.w.fit$sigma2[1]+cne.w.fit$sigma2[2])/
             (cne.w.fit$sigma2[1]+cne.w.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cne.w.fit$sigma2[1])/(cne.w.fit$sigma2[1]+cne.w.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((cne.w.fit$sigma2[2])/(cne.w.fit$sigma2[1]+cne.w.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnew4 <- lm(residuals.rma(cne.w.fit)~wolbachia3$wi)
summary(egger.cnew4)

##  Natural enemy as fixed factor, study ID and effect size ID as random factors ##

table(wolbachia3$NE_group, useNA ="always")

pd<-split(wolbachia3,wolbachia3$NE_group)  
table(pd$Bacteria$Paper_ID)
table(pd$Nematode$Paper_ID)
table(pd$Parasitoid$Paper_ID)
table(pd$Protozoan$Paper_ID)
table(pd$Virus$Paper_ID)

## Model ##

cne.w.ne<-rma.mv(yi = yi, V = vi, mods = ~NE_group-1, 
                 random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = wolbachia3)
summary(cne.w.ne) 


## Calculating total heterogeneity ##
I2.total<-((cne.w.ne$sigma2[1]+cne.w.ne$sigma2[2])/
             (cne.w.ne$sigma2[1]+cne.w.ne$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((cne.w.ne$sigma2[1])/(cne.w.ne$sigma2[1]+cne.w.ne$sigma2[2]+var))*100
Between

## Calculating heterogeneity within##
Within<-((cne.w.ne$sigma2[2])/(cne.w.ne$sigma2[1]+cne.w.ne$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnew5 <- lm(residuals.rma(cne.w.ne)~wolbachia3$wi)
summary(egger.cnew5)

########################### Subgroup analyses: aphids only ################################# 

aphid3<-cne.data1[cne.data1$Host_family=="Aphididae",]

aphid3$Paper_ES<-as.factor(aphid3$Paper_ES)
aphid3$Paper_ID<-as.factor(aphid3$Paper_ID)
aphid3$Symbio_type<-as.factor(aphid3$Symbio_type)
aphid3$Symbiont_spp<-as.factor(aphid3$Symbiont_spp)
aphid3$Group_measure<-as.factor(aphid3$Group_measure)
aphid3$NE_group<-as.factor(aphid3$NE_group)

aphid3<-na.omit(aphid3)

## Overall model - no moderators, study ID and effect size ID as random factors ##

ne.aph.overall<-rma.mv(yi = yi, V = vi, random=list(~1 | Paper_ID, ~1 | Paper_ES),
                       method = "REML", data = aphid3)
summary(ne.aph.overall)

## calculating variance ##
var<-sum(aphid3$wi*(length(aphid3$wi)-1))/(sum(aphid3$wi)^2-sum(aphid3$wi^2))
var

## Calculating total heterogeneity ##
I2.total<-((ne.aph.overall$sigma2[1]+ne.aph.overall$sigma2[2])/
             (ne.aph.overall$sigma2[1]+ne.aph.overall$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ne.aph.overall$sigma2[1])/(ne.aph.overall$sigma2[1]+ne.aph.overall$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ne.aph.overall$sigma2[2])/(ne.aph.overall$sigma2[1]+ne.aph.overall$sigma2[2]+var))*100
Within  

## Publication bias ##

egger.cnea1 <- lm(residuals.rma(ne.aph.overall)~aphid3$wi)
summary(egger.cnea1)

## Symbiont type as fixed factor, study ID and effect size ID as random factors ##

table(aphid3$Symbio_type, useNA ="always") ##all symbionts are natural

## Fitness measure as fixed factor, study ID and effect size ID as random factors ##

table(aphid3$Group_measure, useNA ="always")

pd<-split(aphid3,aphid3$Group_measure)  
table(pd$Body_size$Paper_ID)
table(pd$Survival$Paper_ID)

## Model ##

ne.aph.fit<-rma.mv(yi = yi, V = vi, mods = ~Group_measure-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid3)
summary(ne.aph.fit) 


## Calculating total heterogeneity ##
I2.total<-((ne.aph.fit$sigma2[1]+ne.aph.fit$sigma2[2])/
             (ne.aph.fit$sigma2[1]+ne.aph.fit$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ne.aph.fit$sigma2[1])/(ne.aph.fit$sigma2[1]+ne.aph.fit$sigma2[2]+var))*100
Between

## Calculating heterogeneity within studies ##
Within<-((ne.aph.fit$sigma2[2])/(ne.aph.fit$sigma2[1]+ne.aph.fit$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnea3 <- lm(residuals.rma(ne.aph.fit)~aphid3$wi)
summary(egger.cnea3)

##  Symbiont species as fixed factor, study ID and effect size ID as random factors ##

table(aphid3$Symbiont_spp, useNA ="always")
pd<-split(aphid3,aphid3$Symbiont_spp)  
table(pd$Hamiltonella$Paper_ID)
table(pd$Other$Paper_ID)
table(pd$Spiroplasma$Paper_ID)

## Model ##

ne.aph.spp<-rma.mv(yi = yi, V = vi, mods = ~Symbiont_spp-1, 
                   random = list(~1 | Paper_ID, ~1 | Paper_ES), method = "REML", data = aphid3)
summary(ne.aph.spp)


## Calculating total heterogeneity ##
I2.total<-((ne.aph.spp$sigma2[1]+ne.aph.spp$sigma2[2])/
             (ne.aph.spp$sigma2[1]+ne.aph.spp$sigma2[2]+var))*100
I2.total

## Calculating heterogeneity between studies ##
Between<-((ne.aph.spp$sigma2[1])/(ne.aph.spp$sigma2[1]+ne.aph.spp$sigma2[2]+var))*100
Between

## Calculating heterogeneity within ##
Within<-((ne.aph.spp$sigma2[2])/(ne.aph.spp$sigma2[1]+ne.aph.spp$sigma2[2]+var))*100
Within  

## Publication bias ##
egger.cnea4 <- lm(residuals.rma(ne.aph.spp)~aphid3$wi)
summary(egger.cnea4)
