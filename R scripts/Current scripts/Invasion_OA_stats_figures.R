
# Read in packages --------------------------------------------------------
#### read in useful packages
library(here)
library(mgcv)
library(gratia)
library(tidyr)
library(bbmle) 
library(glmmTMB)
library(doBy)
library(plyr)
library(dplyr)
library(ggplot2) 
library(doBy)
library(grid)
library(lmtest)
library(fitdistrplus)
library(visreg)
library(lme4)
library(coefplot)
library(arm)
library(lmerTest)
library(boot)
library(MASS)
require(scales)
library(car)
library(knitr)
library(tidyverse)
library(kableExtra)
library(multcomp)
library(arm) ## for sim()
library(descr)  ## for LogRegR2
library(reshape2)
library(grid)
library(DHARMa)
library(gap)
library(qrnn)
library(colorspace)
library(cowplot)
library(devtools)
library(installr)


# dataframe management ----------------------------------------------------

invasion.exp.data<-read.csv("C:Biological data/allweeks_cover_counts_without_pres.csv",stringsAsFactors = FALSE, na.strings = c("NA","") )

invasion.exp.data$num.nudi<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.eggs+invasion.exp.data$nudi.hatched

invasion.exp.data.16<-invasion.exp.data %>% filter(Week==16)
invasion.exp.data.8<-invasion.exp.data %>% filter(Week==8)

#ordered and unordered factors
invasion.exp.data.16$oTreatment<-factor(invasion.exp.data.16$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=TRUE)
invasion.exp.data.16$Treatment<-factor(invasion.exp.data.16$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=FALSE)

#I think I need just Invasives if we're going with CO2.Treatment
invasion.exp.data.16$oInvasives<-factor(invasion.exp.data.16$Invasives, levels=c("Absent", "Present"), ordered=TRUE)
invasion.exp.data.16$Invasives<-factor(invasion.exp.data.16$Invasives, levels=c("Absent", "Present"), ordered=FALSE)

invasion.exp.data.16$oCO2.Treatment<-factor(invasion.exp.data.16$CO2.Treatment, levels=c("AIR", "CO2"), ordered=TRUE)
invasion.exp.data.16$CO2.Treatment<-factor(invasion.exp.data.16$CO2.Treatment, levels=c("AIR", "CO2"), ordered=FALSE)




# #some variables need to be read over
# invasion.exp.data.16$caprellid.percent<-invasion.exp.data.16$caprellid
# invasion.exp.data.16$hydroid<-invasion.exp.data.16$hydroid
# invasion.exp.data.16$botryllid<-invasion.exp.data.16$botryllid
# invasion.exp.data.16$folliculina<-invasion.exp.data.16$folliculina
# invasion.exp.data.16$membranipora<-invasion.exp.data.16$membranipora
# invasion.exp.data.16$didemnum<-invasion.exp.data.16$didemnum
# invasion.exp.data.16$total<-invasion.exp.data.16$total
# invasion.exp.data.16$bare<-invasion.exp.data.16$bare
# invasion.exp.data.16$occupied.space<-100-invasion.exp.data.16$bare
# invasion.exp.data.16$occupied.space.001<-0.01*(invasion.exp.data.16$occupied.space)
# invasion.exp.data.16$everything.wet.weight<-invasion.exp.data.16$everything.wet.weight
# invasion.exp.data.16$everything.wet.weight.per.1<-(invasion.exp.data.16$everything.wet.weight)/invasion.exp.data.16$occupied.space
# invasion.exp.data.16$Mussel.wet.weight<-invasion.exp.data.16$Mussel.wet.weight
# invasion.exp.data.16$total_dry_biomass<-invasion.exp.data.16$total_dry_biomass
# invasion.exp.data.16$total_dry_biomass_per1<-invasion.exp.data.16$total_dry_biomass/invasion.exp.data.16$occupied.space
# invasion.exp.data.16$hydroid_dry_biomass<-invasion.exp.data.16$hydroid_dry_biomass
# invasion.exp.data.16$caprellid_dry_biomass<-invasion.exp.data.16$caprellid_dry_biomass
# invasion.exp.data.16$tunicate_dry_biomass<-invasion.exp.data.16$tunicate_dry_biomass
# invasion.exp.data.16$hydtobot<-(invasion.exp.data.16$botryllid)/(invasion.exp.data.16$botryllid+invasion.exp.data.16$hydroid)
# invasion.exp.data.16$rest_dry_biomass<-invasion.exp.data.16$rest_dry_biomass
# 
# #small negative biomass is within error of scale - change to zero
# invasion.exp.data.16$hydroid_dry_biomass[invasion.exp.data.16$hydroid_dry_biomass<0]<-0
# invasion.exp.data.16$tunicate_dry_biomass[invasion.exp.data.16$tunicate_dry_biomass<0]<-0
# invasion.exp.data.16$hydtobot_dry_biomass<-(invasion.exp.data.16$tunicate_dry_biomass)/(invasion.exp.data.16$tunicate_dry_biomass+invasion.exp.data.16$hydroid_dry_biomass)
# 
# 
# invasion.exp.data.16$Mussel.wet.weight.per.1<-(invasion.exp.data.16$Mussel.wet.weight)/(invasion.exp.data.16$mussel+1)
# 
# #making it a proportion instead of % cover
# invasion.exp.data.16$caprellid.percent.001<-(0.01*(invasion.exp.data.16$caprellid.percent))+0.01
 invasion.exp.data.16$hydroid.001<-(0.01*(invasion.exp.data.16$hydroid))+0.01
 invasion.exp.data.16$botryllid.001<-(0.01*(invasion.exp.data.16$botryllid))+0.01
 invasion.exp.data.16$membranipora.001<-(0.01*(invasion.exp.data.16$membranipora))+0.01
 invasion.exp.data.16$didemnum<-invasion.exp.data.16$white.bryo
 invasion.exp.data.16$num.red.bryoporella<-invasion.exp.data.16$red.bryo
 invasion.exp.data.16$folliculina<-invasion.exp.data.16$protozoa
 invasion.exp.data.16$folliculina.001<-(0.01*(invasion.exp.data.16$folliculina))+0.01
 
  invasion.exp.data.16$didemnum.001<-(0.01*(invasion.exp.data.16$didemnum))+0.01

# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
invasion.exp.data.16_zscores<-invasion.exp.data.16
#invasion.exp.data.16_zscores$hydrogen.concentration<-scale(invasion.exp.data.16$hydrogen.concentration, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$av.pH<-scale(invasion.exp.data.16$av.pH, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$min.10.pH<-scale(invasion.exp.data.16$min.10.pH, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$Mesocosm <- as.factor(invasion.exp.data.16$Mesocosm)
invasion.exp.data.16_zscores$av.pH.unscaled <-invasion.exp.data.16_zscores$av.pH * attr(invasion.exp.data.16_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$av.pH, 'scaled:center')
invasion.exp.data.16_zscores$min.10.pH.unscaled <-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# Visualizing histograms of pH to use as continuous vs. discrete variable
invasion.exp.data.16_pres<-invasion.exp.data.16 %>% filter(Invasives=="Present")
invasion.exp.data.16_abs<-invasion.exp.data.16 %>% filter(Invasives=="Absent")




hist(invasion.exp.data.16_abs$min.10.pH, breaks=5)#5
hist(invasion.exp.data.16_pres$min.10.pH, breaks=18)#18 for the whole 16 weeks
#these are actually the exact same values b/c same mesocosms


 ggplot(invasion.exp.data.16, aes(x=min.10.pH))+geom_density()+theme_classic()

# Week 12
#ordered and unordered factors
 invasion.exp.data.8$oTreatment<-factor(invasion.exp.data.8$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=TRUE)
 invasion.exp.data.8$Treatment<-factor(invasion.exp.data.8$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=FALSE)
 
#I think I need just Invasives if we're going with CO2.Treatment
 invasion.exp.data.8$oInvasives<-factor(invasion.exp.data.8$Invasives, levels=c("Absent", "Present"), ordered=TRUE)
 invasion.exp.data.8$Invasives<-factor(invasion.exp.data.8$Invasives, levels=c("Absent", "Present"), ordered=FALSE)
 
 invasion.exp.data.8$oCO2.Treatment<-factor(invasion.exp.data.8$CO2.Treatment, levels=c("AIR", "CO2"), ordered=TRUE)
 invasion.exp.data.8$CO2.Treatment<-factor(invasion.exp.data.8$CO2.Treatment, levels=c("AIR", "CO2"), ordered=FALSE)
 
  
 # #making it a proportion instead of % cover
 # invasion.exp.data.8$caprellid.percent.001<-(0.01*(invasion.exp.data.8$caprellid.percent))+0.01
 invasion.exp.data.8$hydroid.001<-(0.01*(invasion.exp.data.8$hydroid))+0.01
 invasion.exp.data.8$botryllid.001<-(0.01*(invasion.exp.data.8$botryllid))+0.01
 invasion.exp.data.8$membranipora.001<-(0.01*(invasion.exp.data.8$membranipora))+0.01
 invasion.exp.data.8$didemnum<-invasion.exp.data.8$white.bryo
 invasion.exp.data.8$num.red.bryoporella<-invasion.exp.data.8$red.bryo
 invasion.exp.data.8$folliculina<-invasion.exp.data.8$protozoa
 invasion.exp.data.8$folliculina.001<-(0.01*(invasion.exp.data.8$folliculina))+0.01
 invasion.exp.data.8$didemnum.001<-(0.01*(invasion.exp.data.8$didemnum))+0.01

 # need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
 invasion.exp.data.8_zscores<-invasion.exp.data.8
 #invasion.exp.data.8_zscores$hydrogen.concentration<-scale(invasion.exp.data.8$hydrogen.concentration, center=TRUE, scale=TRUE)
 invasion.exp.data.8_zscores$av.pH<-scale(invasion.exp.data.8$av.pH, center=TRUE, scale=TRUE)
 invasion.exp.data.8_zscores$min.10.pH<-scale(invasion.exp.data.8$min.10.pH, center=TRUE, scale=TRUE)
 invasion.exp.data.8_zscores$Mesocosm <- as.factor(invasion.exp.data.8$Mesocosm)
 invasion.exp.data.8_zscores$av.pH.unscaled <-invasion.exp.data.8_zscores$av.pH * attr(invasion.exp.data.8_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$av.pH, 'scaled:center')
 invasion.exp.data.8_zscores$min.10.pH.unscaled <-invasion.exp.data.8_zscores$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')
 
 
# Visualizing histograms of pH to use as continuous vs. discrete variable
 invasion.exp.data.8_pres<-invasion.exp.data.8 %>% filter(Invasives=="Present")
 invasion.exp.data.8_abs<-invasion.exp.data.8 %>% filter(Invasives=="Absent")
 
 hist(invasion.exp.data.8_abs$min.10.pH, breaks=8)#8
 hist(invasion.exp.data.8_pres$min.10.pH, breaks=8)#8
 #these are actually the exact same values b/c same mesocosms

# Notes on contrasts ------------------------------------------------------

#Notes on ordered factors: the factor Invasives is ordered
#options(contrasts = c("contr.sum", "contr.poly"))
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.


# GAM information and resources---------------------------------------------------------------------
#Gavin Simpson's gratia - used to visualize gams

#library(gratia)

#We have a varying coefficient model aka ANCOVA, so use "by"

#We're going to use gam(...,select=TRUE), automatic model selection via null space penalization

##Invasives is an ordered factor - meaning "none" is the base level and the other two relate to that....
# see https://www.fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/



# Plotting settings -------------------------------------------------------

colorset_invasives = c("Present"="#5EC2DA" ,"Absent"="#EB549A")
theme_set(theme_classic(base_size = 12))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))



# GAM beta hydroid / gam.16.beta.hydroid --------------------------------------------------------

#distributions - binomial or beta with various families

gam.16.binomial.hydroid<- gam(formula = cbind(hydroid, 100-hydroid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")
gam.16.beta.hydroid<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.hydroid.1<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.hydroid.2<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.hydroid.3<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.16.beta.hydroid, gam.16.beta.hydroid.1, gam.16.beta.hydroid.2, gam.16.beta.hydroid.3, gam.16.binomial.hydroid)
#cauchit is the best

plot(gam.16.beta.hydroid.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.beta.hydroid.3)
qq_plot(gam.16.beta.hydroid.3, method = 'simulate')
k.check(gam.16.beta.hydroid.3)

gam.16.beta.hydroid.3.unordered<- gam(hydroid.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.16.beta.hydroid.3)
summary(gam.16.beta.hydroid.3.unordered)
# need an unordered and an unordered model to get estimates for overall Invasives effect


### plotting based on gam confidence intervals
#using code from https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

fam.gam.16.hydroid <- family(gam.16.beta.hydroid.3)
fam.gam.16.hydroid 
str(fam.gam.16.hydroid )
ilink.gam.16.hydroid <- fam.gam.16.hydroid$linkinv
ilink.gam.16.hydroid

invasion.exp.data.16_zscores$min.10.pH.unscaled <-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data.16_zscores)

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)
            
mod.hydroid<-gam.16.beta.hydroid.3
ndata.16.hydroid <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.hydroid <- add_column(ndata.16.hydroid, fit = predict(mod.hydroid, newdata = ndata.16.hydroid, type = 'response'))
ndata.16.hydroid <- bind_cols(ndata.16.hydroid, setNames(as_tibble(predict(mod.hydroid, ndata.16.hydroid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform
ndata.16.hydroid <- mutate(ndata.16.hydroid,
                fit_resp  = ilink.gam.16.hydroid(fit_link),
                right_upr = ilink.gam.16.hydroid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.16.hydroid(fit_link - (2 * se_link)))

#make sure pH is unscaled in the plot
ndata.16.hydroid$min.10.pH.unscaled<-ndata.16.hydroid$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


plt.gam.16.hydroid <- ggplot(ndata.16.hydroid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydroid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Obelia")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.hydroid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.gam.16.hydroid 
ggsave("C:Graphs April 2020//hydroid_pred.png")



# GAM beta botryllus / gam.16.beta.botryllid.2 ----------------------------------------------------

#binomial first
gam.16.binomial.botryllid<- gam(formula = cbind(botryllid, 100-botryllid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.16.beta.botryllid<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.1<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.2<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.3<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.16.beta.botryllid, gam.16.beta.botryllid.1, gam.16.beta.botryllid.2, gam.16.beta.botryllid.3, gam.16.binomial.botryllid)


plot(gam.16.beta.botryllid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.beta.botryllid)
qq_plot(gam.16.beta.botryllid, method = 'simulate')
k.check(gam.16.beta.botryllid)
summary(gam.16.beta.botryllid)

gam.16.beta.botryllid.unordered<- gam(botryllid.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)


fam.gam.16.botryllid <- family(gam.16.beta.botryllid)
fam.gam.16.botryllid
ilink.gam.16.botryllid<- fam.gam.16.botryllid$linkinv
ilink.gam.16.botryllid
want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)


mod.botryllid<-gam.16.beta.botryllid
ndata.16.botryllid <- with(invasion.exp.data.16_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                        length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.botryllidel for the new data
ndata.16.botryllid <- add_column(ndata.16.botryllid, fit = predict(mod.botryllid, newdata = ndata.16.botryllid, type = 'response'))


ndata.16.botryllid <- bind_cols(ndata.16.botryllid, setNames(as_tibble(predict(mod.botryllid, ndata.16.botryllid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.botryllid <- mutate(ndata.16.botryllid,
                fit_resp  = ilink.gam.16.botryllid(fit_link),
                right_upr = ilink.gam.16.botryllid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.16.botryllid(fit_link - (2 * se_link)))


invasion.exp.data.16_zscores$min.10.pH.unscaled<-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

ndata.16.botryllid$min.10.pH.unscaled<-ndata.16.botryllid$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.botryllid.16 <- ggplot(ndata.16.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.botryllid.16
ggsave("C:Graphs April 2020//botryllid_pred.16.png")


# GAM beta botryllus / gam.8.beta.botryllid.2 ----------------------------------------------------

#binomial first
gam.8.binomial.botryllid<- gam(formula = cbind(botryllid, 100-botryllid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.8.beta.botryllid<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.1<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.2<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.3<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.8.beta.botryllid, gam.8.beta.botryllid.1, gam.8.beta.botryllid.2, gam.8.beta.botryllid.3, gam.8.binomial.botryllid)
#.3 is best

plot(gam.8.beta.botryllid.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.beta.botryllid.3)
qq_plot(gam.8.beta.botryllid.3, method = 'simulate')
k.check(gam.8.beta.botryllid.3)
summary(gam.8.beta.botryllid.3)

gam.8.beta.botryllid.3.unordered<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=Invasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


fam.gam.8.botryllid <- family(gam.8.beta.botryllid.3)
fam.gam.8.botryllid
ilink.gam.8.botryllid<- fam.gam.8.botryllid$linkinv
ilink.gam.8.botryllid
want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)


mod.botryllid<-gam.8.beta.botryllid.3
ndata.8.botryllid <- with(invasion.exp.data.8_zscores, 
                           data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.botryllidel for the new data
ndata.8.botryllid <- add_column(ndata.8.botryllid, fit = predict(mod.botryllid, newdata = ndata.8.botryllid, type = 'response'))


ndata.8.botryllid <- bind_cols(ndata.8.botryllid, setNames(as_tibble(predict(mod.botryllid, ndata.8.botryllid, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.botryllid <- mutate(ndata.8.botryllid,
                             fit_resp  = ilink.gam.8.botryllid(fit_link),
                             right_upr = ilink.gam.8.botryllid(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.8.botryllid(fit_link - (2 * se_link)))


invasion.exp.data.8_zscores$min.10.pH.unscaled<-invasion.exp.data.8_zscores$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

ndata.8.botryllid$min.10.pH.unscaled<-ndata.8.botryllid$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.botryllid.8 <- ggplot(ndata.8.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.botryllid.8
ggsave("C:Graphs April 2020//botryllid_pred.8.png")



# botryllid + time -- not working with gam b/c too few --------------------------------------------------------
colorset_treatment = c("AIRAbsent"="#5EC2DA","CO2.TreatmentAbsent"="#EBC915", "CO2.TreatmentPresent"="#EB549A", "AIRPresent"="#550133")
#   
theme_set(theme_classic(base_size = 12))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

gam.16.binomial.botryllid<- gam(formula = cbind(botryllid, 100-botryllid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.16.beta.botryllid<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.1<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.2<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.16.beta.botryllid.3<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.16.beta.botryllid, gam.16.beta.botryllid.1, gam.16.beta.botryllid.2, gam.16.beta.botryllid.3, gam.16.binomial.botryllid)
#cauchit is the best

plot(gam.16.beta.botryllid.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.beta.botryllid.3)
qq_plot(gam.16.beta.botryllid.3, method = 'simulate')
k.check(gam.16.beta.botryllid.3)
summary(gam.16.beta.botryllid.3)

gam.16.beta.botryllid.3.unordered<- gam(botryllid.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)


fam.gam.16.botryllid <- family(gam.16.beta.botryllid.3)
fam.gam.16.botryllid
ilink.gam.16.botryllid<- fam.gam.16.botryllid$linkinv
ilink.gam.16.botryllid
want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)


mod.botryllid<-gam.16.beta.botryllid.3
ndata.16.botryllid <- with(invasion.exp.data.16_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                   length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.botryllidel for the new data
ndata.16.botryllid <- add_column(ndata.16.botryllid, fit = predict(mod.botryllid, newdata = ndata.16.botryllid, type = 'response'))


ndata.16.botryllid <- bind_cols(ndata.16.botryllid, setNames(as_tibble(predict(mod.botryllid, ndata.16.botryllid, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.botryllid <- mutate(ndata.16.botryllid,
                          fit_resp  = ilink.gam.16.botryllid(fit_link),
                          right_upr = ilink.gam.16.botryllid(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.16.botryllid(fit_link - (2 * se_link)))


invasion.exp.data.16_zscores$min.10.pH.unscaled<-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

ndata.16.botryllid$min.10.pH.unscaled<-ndata.16.botryllid$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.botryllid <- ggplot(ndata.16.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.botryllid
ggsave("C:Graphs April 2020//botryllid_pred.png")





# GAM beta folliculina / gam.16.beta.folliculina -----------------------------------------------------------

gam.16.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")
gam.16.beta.folliculina<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.folliculina.1<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.folliculina.2<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.folliculina.3<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.16.beta.folliculina, gam.16.beta.folliculina.1, gam.16.beta.folliculina.2,gam.16.binomial.folliculina, gam.16.beta.folliculina.3)
#simplest logit


plot(gam.16.beta.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.beta.folliculina)
qq_plot(gam.16.beta.folliculina, method = 'simulate')
k.check(gam.16.beta.folliculina)
summary(gam.16.beta.folliculina)

gam.16.beta.folliculina.unordered<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.16.folliculina <- family(gam.16.beta.folliculina)
fam.gam.16.folliculina
str(fam.gam.16.folliculina)
ilink.gam.16.folliculina<- fam.gam.16.folliculina$linkinv
ilink.gam.16.folliculina


mod.folliculina<-gam.16.beta.folliculina
ndata.16.folliculina <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.folliculinael for the new data
ndata.16.folliculina <- add_column(ndata.16.folliculina, fit = predict(mod.folliculina, newdata = ndata.16.folliculina, type = 'response'))


ndata.16.folliculina <- bind_cols(ndata.16.folliculina, setNames(as_tibble(predict(mod.folliculina, ndata.16.folliculina, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.folliculina <- mutate(ndata.16.folliculina,
                          fit_resp  = ilink.gam.16.folliculina(fit_link),
                          right_upr = ilink.gam.16.folliculina(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.16.folliculina(fit_link - (2 * se_link)))


ndata.16.folliculina$min.10.pH.unscaled<-ndata.16.folliculina$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.folliculina.16 <- ggplot(ndata.16.folliculina, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = folliculina.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.folliculina.16
ggsave("C:Graphs April 2020//folliculina_pred.16.png")

# GAM beta folliculina / gam.8.beta.folliculina -----------------------------------------------------------

gam.8.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.folliculina<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.folliculina.1<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.folliculina.2<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.folliculina.3<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.8.beta.folliculina, gam.8.beta.folliculina.1, gam.8.beta.folliculina.2,gam.8.binomial.folliculina, gam.8.beta.folliculina.3)
#simplest logit


plot(gam.8.beta.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.beta.folliculina)
qq_plot(gam.8.beta.folliculina, method = 'simulate')
k.check(gam.8.beta.folliculina)
summary(gam.8.beta.folliculina)

gam.8.beta.folliculina.unordered<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.8.folliculina <- family(gam.8.beta.folliculina)
fam.gam.8.folliculina
str(fam.gam.8.folliculina)
ilink.gam.8.folliculina<- fam.gam.8.folliculina$linkinv
ilink.gam.8.folliculina


mod.folliculina<-gam.8.beta.folliculina
ndata.8.folliculina <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.folliculinael for the new data
ndata.8.folliculina <- add_column(ndata.8.folliculina, fit = predict(mod.folliculina, newdata = ndata.8.folliculina, type = 'response'))


ndata.8.folliculina <- bind_cols(ndata.8.folliculina, setNames(as_tibble(predict(mod.folliculina, ndata.8.folliculina, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.folliculina <- mutate(ndata.8.folliculina,
                               fit_resp  = ilink.gam.8.folliculina(fit_link),
                               right_upr = ilink.gam.8.folliculina(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.folliculina(fit_link - (2 * se_link)))


ndata.8.folliculina$min.10.pH.unscaled<-ndata.8.folliculina$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.folliculina.8 <- ggplot(ndata.8.folliculina, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = folliculina.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.folliculina.8
ggsave("C:Graphs April 2020//folliculina_pred.8.png")


# GAM beta membranipora / gam.16.beta.membranipora --------------------------------------------------------

#binomial first
gam.16.binomial.membranipora<- gam(formula = cbind(membranipora, 100-membranipora)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.16.beta.membranipora<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.membranipora.1<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.membranipora.2<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.membranipora.3<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.16.beta.membranipora, gam.16.beta.membranipora.1, gam.16.beta.membranipora.2, gam.16.binomial.membranipora, gam.16.beta.membranipora.3)
#cauchit is best


plot(gam.16.beta.membranipora.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.beta.membranipora.3)
qq_plot(gam.16.beta.membranipora.3, method = 'simulate')
k.check(gam.16.beta.membranipora.3)
summary(gam.16.beta.membranipora.3)
vis.gam(gam.16.beta.membranipora.3)

gam.16.beta.membranipora.3.unordered<- gam(membranipora.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


fam.gam.16.membranipora <- family(gam.16.beta.membranipora.3)
fam.gam.16.membranipora
ilink.gam.16.membranipora<- fam.gam.16.membranipora$linkinv
ilink.gam.16.membranipora


mod.membranipora<-gam.16.beta.membranipora.3
ndata.16.membranipora <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.membraniporael for the new data
ndata.16.membranipora <- add_column(ndata.16.membranipora, fit = predict(mod.membranipora, newdata = ndata.16.membranipora, type = 'response'))


ndata.16.membranipora <- bind_cols(ndata.16.membranipora, setNames(as_tibble(predict(mod.membranipora, ndata.16.membranipora, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.membranipora <- mutate(ndata.16.membranipora,
                          fit_resp  = ilink.gam.16.membranipora(fit_link),
                          right_upr = ilink.gam.16.membranipora(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.16.membranipora(fit_link - (2 * se_link)))

ndata.16.membranipora$min.10.pH.unscaled<-ndata.16.membranipora$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.membranipora.16 <- ggplot(ndata.16.membranipora, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = membranipora.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.membranipora.16
ggsave("C:Graphs April 2020//membranipora_pred.16.png")

# GAM beta membranipora / gam.8.beta.membranipora --------------------------------------------------------

#binomial first
gam.8.binomial.membranipora<- gam(formula = cbind(membranipora, 100-membranipora)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.8.beta.membranipora<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.membranipora.1<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.membranipora.2<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.membranipora.3<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.8.beta.membranipora, gam.8.beta.membranipora.1, gam.8.beta.membranipora.2, gam.8.binomial.membranipora, gam.8.beta.membranipora.3)
#cauchit is best


plot(gam.8.beta.membranipora.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.beta.membranipora.3)
qq_plot(gam.8.beta.membranipora.3, method = 'simulate')
k.check(gam.8.beta.membranipora.3)
summary(gam.8.beta.membranipora.3)
vis.gam(gam.8.beta.membranipora.3)

gam.8.beta.membranipora.3.unordered<- gam(membranipora.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


fam.gam.8.membranipora <- family(gam.8.beta.membranipora.3)
fam.gam.8.membranipora
ilink.gam.8.membranipora<- fam.gam.8.membranipora$linkinv
ilink.gam.8.membranipora


mod.membranipora<-gam.8.beta.membranipora.3
ndata.8.membranipora <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.membraniporael for the new data
ndata.8.membranipora <- add_column(ndata.8.membranipora, fit = predict(mod.membranipora, newdata = ndata.8.membranipora, type = 'response'))


ndata.8.membranipora <- bind_cols(ndata.8.membranipora, setNames(as_tibble(predict(mod.membranipora, ndata.8.membranipora, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.membranipora <- mutate(ndata.8.membranipora,
                                fit_resp  = ilink.gam.8.membranipora(fit_link),
                                right_upr = ilink.gam.8.membranipora(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.8.membranipora(fit_link - (2 * se_link)))

ndata.8.membranipora$min.10.pH.unscaled<-ndata.8.membranipora$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.membranipora.8 <- ggplot(ndata.8.membranipora, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = membranipora.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.membranipora.8
ggsave("C:Graphs April 2020//membranipora_pred.8.png")
 

# GAM negbin barnacles / gam.16.nb.num.barn -----------------------------------------------------------
nbinom.16.barn <- fitdistr(invasion.exp.data.16_zscores$num.barn, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.barn, "nbinom", size = nbinom.16.barn$estimate[[1]], mu = nbinom.16.barn$estimate[[2]])

#negative binomial 
gam.16.nb.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.barn$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.barn.1<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.16.nb.num.barn, gam.16.nb.num.barn.1, gam.16.poisson.num.barn)

plot(gam.16.poisson.num.barn, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.poisson.num.barn)
qq_plot(gam.16.poisson.num.barn, method = 'simulate')
#looks really good
k.check(gam.16.poisson.num.barn)
summary(gam.16.poisson.num.barn)


#a few outside the area
#appraise a bit funnelly

gam.16.poisson.num.barn.unordered<- gam(num.barn ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.barnacles$estimate[[1]]), select=TRUE, method="REML")
fam.gam.16.num.barn <- family(gam.16.poisson.num.barn)
fam.gam.16.num.barn
str(fam.gam.16.num.barn)
ilink.gam.16.num.barn<- fam.gam.16.num.barn$linkinv
ilink.gam.16.num.barn


mod.num.barn<-gam.16.poisson.num.barn
ndata.16.num.barn <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.16.num.barn <- add_column(ndata.16.num.barn, fit = predict(mod.num.barn, newdata = ndata.16.num.barn, type = 'response'))
predict(mod.num.barn, newdata = ndata.16.num.barn, type = 'response')
ndata.16.num.barn <- bind_cols(ndata.16.num.barn, setNames(as_tibble(predict(mod.num.barn, ndata.16.num.barn, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))
## create the interval and backtransform
ndata.16.num.barn <- mutate(ndata.16.num.barn,
                                  fit_resp  = ilink.gam.16.num.barn(fit_link),
                                  right_upr = ilink.gam.16.num.barn(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.16.num.barn(fit_link - (2 * se_link)))
ndata.16.num.barn$min.10.pH.unscaled<-ndata.16.num.barn$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.barn.16 <- ggplot(ndata.16.num.barn, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.barn, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.barn.16
ggsave("C:Graphs April 2020//num.barn_pred.16.png")

# GAM negbin barnacles / gam.8.nb.num.barn -----------------------------------------------------------
nbinom.8.barn <- fitdistr(invasion.exp.data.8_zscores$num.barn, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.barn, "nbinom", size = nbinom.8.barn$estimate[[1]], mu = nbinom.8.barn$estimate[[2]])

#negative binomial 
gam.8.nb.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.barn$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.barn.1<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.barn, gam.8.nb.num.barn.1, gam.8.poisson.num.barn)
#poisson is better

plot(gam.8.poisson.num.barn, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.poisson.num.barn)
qq_plot(gam.8.poisson.num.barn, method = 'simulate')
k.check(gam.8.poisson.num.barn)
summary(gam.8.poisson.num.barn)

gam.8.poisson.num.barn.unordered<- gam(num.barn ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")
fam.gam.8.num.barn <- family(gam.8.poisson.num.barn)
ilink.gam.8.num.barn<- fam.gam.8.num.barn$linkinv
mod.num.barn<-gam.8.poisson.num.barn
ndata.8.num.barn <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                   length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.8.num.barn <- add_column(ndata.8.num.barn, fit = predict(mod.num.barn, newdata = ndata.8.num.barn, type = 'response'))
predict(mod.num.barn, newdata = ndata.8.num.barn, type = 'response')
ndata.8.num.barn <- bind_cols(ndata.8.num.barn, setNames(as_tibble(predict(mod.num.barn, ndata.8.num.barn, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform
ndata.8.num.barn <- mutate(ndata.8.num.barn,
                            fit_resp  = ilink.gam.8.num.barn(fit_link),
                            right_upr = ilink.gam.8.num.barn(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.8.num.barn(fit_link - (2 * se_link)))
ndata.8.num.barn$min.10.pH.unscaled<-ndata.8.num.barn$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')


# plot 
plt.num.barn.8 <- ggplot(ndata.8.num.barn, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.barn, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.barn.8
ggsave("C:Graphs April 2020//num.barn_pred.8.png")


# GAM negbin num.white.bryo / gam.16.nb.num.white.bryo ----------------------------------------------------------
nbinom.16.num.white.bryo <- fitdistr(invasion.exp.data.16_zscores$num.white.bryo, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.white.bryo, "nbinom", size = nbinom.16.num.white.bryo$estimate[[1]], mu = nbinom.16.num.white.bryo$estimate[[2]])
#getting theta

gam.16.nb.num.white.bryo<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.white.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.white.bryo.1<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.white.bryo<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.16.nb.num.white.bryo.1,gam.16.nb.num.white.bryo,gam.16.poisson.num.white.bryo)
#used estimated theta

plot(gam.16.nb.num.white.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.nb.num.white.bryo)
qq_plot(gam.16.nb.num.white.bryo, method = 'simulate')
k.check(gam.16.nb.num.white.bryo)
summary(gam.16.nb.num.white.bryo)

#a few outside the area
#appraise a bit funnelly

gam.16.nb.num.white.bryo.unordered<- gam(num.white.bryo ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.white.bryo$estimate[[1]]), select=TRUE, method="REML")

fam.gam.16.num.white.bryo <- family(gam.16.nb.num.white.bryo)
fam.gam.16.num.white.bryo
str(fam.gam.16.num.white.bryo)
ilink.gam.16.num.white.bryo<- fam.gam.16.num.white.bryo$linkinv
ilink.gam.16.num.white.bryo

mod.num.white.bryo<-gam.16.nb.num.white.bryo
ndata.16.num.white.bryo <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.16.num.white.bryo <- add_column(ndata.16.num.white.bryo, fit = predict(mod.num.white.bryo, newdata = ndata.16.num.white.bryo, type = 'response'))
predict(mod.num.white.bryo, newdata = ndata.16.num.white.bryo, type = 'response')
ndata.16.num.white.bryo <- bind_cols(ndata.16.num.white.bryo, setNames(as_tibble(predict(mod.num.white.bryo, ndata.16.num.white.bryo, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))
## create the interval and backtransform
ndata.16.num.white.bryo <- mutate(ndata.16.num.white.bryo,
                               fit_resp  = ilink.gam.16.num.white.bryo(fit_link),
                               right_upr = ilink.gam.16.num.white.bryo(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.16.num.white.bryo(fit_link - (2 * se_link)))
ndata.16.num.white.bryo$min.10.pH.unscaled<-ndata.16.num.white.bryo$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.num.white.bryo.16 <- ggplot(ndata.16.num.white.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.white.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.white.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.white.bryo.16
ggsave("C:Graphs April 2020//num.white.bryo_pred.16.png")


# GAM negbin num.white.bryo / gam.8.nb.num.white.bryo ----------------------------------------------------------

nbinom.8.num.white.bryo <- fitdistr(invasion.exp.data.8_zscores$num.white.bryo, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.white.bryo, "nbinom", size = nbinom.8.num.white.bryo$estimate[[1]], mu = nbinom.8.num.white.bryo$estimate[[2]])
#getting theta

gam.8.nb.num.white.bryo<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.white.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.white.bryo.1<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.white.bryo<- gam(num.white.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.white.bryo.1,gam.8.nb.num.white.bryo,gam.8.poisson.num.white.bryo)
#poisson is better for week 12

plot(gam.8.poisson.num.white.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.poisson.num.white.bryo)
#not that many points
qq_plot(gam.8.poisson.num.white.bryo, method = 'simulate')
k.check(gam.8.poisson.num.white.bryo)
summary(gam.8.poisson.num.white.bryo)

gam.8.poisson.num.white.bryo.unordered<- gam(num.white.bryo ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

fam.gam.8.num.white.bryo <- family(gam.8.poisson.num.white.bryo)
ilink.gam.8.num.white.bryo<- fam.gam.8.num.white.bryo$linkinv
mod.num.white.bryo<-gam.8.poisson.num.white.bryo
ndata.8.num.white.bryo <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                         length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.8.num.white.bryo <- add_column(ndata.8.num.white.bryo, fit = predict(mod.num.white.bryo, newdata = ndata.8.num.white.bryo, type = 'response'))

predict(mod.num.white.bryo, newdata = ndata.8.num.white.bryo, type = 'response')
ndata.8.num.white.bryo <- bind_cols(ndata.8.num.white.bryo, setNames(as_tibble(predict(mod.num.white.bryo, ndata.8.num.white.bryo, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform
ndata.8.num.white.bryo <- mutate(ndata.8.num.white.bryo,
                                  fit_resp  = ilink.gam.8.num.white.bryo(fit_link),
                                  right_upr = ilink.gam.8.num.white.bryo(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.8.num.white.bryo(fit_link - (2 * se_link)))
ndata.8.num.white.bryo$min.10.pH.unscaled<-ndata.8.num.white.bryo$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.white.bryo.8 <- ggplot(ndata.8.num.white.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.white.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.white.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.white.bryo.8
ggsave("C:Graphs April 2020//num.white.bryo_pred.8.png")


# GAM negbin num.red.bryo / gam.16.nb.num.red.bryo --------------------------------------------------------------
nbinom.16.num.red.bryo <- fitdistr(invasion.exp.data.16_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.red.bryo, "nbinom", size = nbinom.16.num.red.bryo$estimate[[1]], mu = nbinom.16.num.red.bryo$estimate[[2]])
#getting theta

gam.16.nb.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.red.bryo, gam.16.nb.num.red.bryo.1, gam.16.poisson.num.red.bryo)


plot(gam.16.poisson.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.16.poisson.num.red.bryo)
qq_plot(gam.16.poisson.num.red.bryo, method = 'simulate')
k.check(gam.16.poisson.num.red.bryo)
summary(gam.16.poisson.num.red.bryo)

gam.16.poisson.num.red.bryo.unordered<- gam(num.red.bryo ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
summary(gam.16.poisson.num.red.bryo.unordered)

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)

fam.gam.16.num.red.bryo <- family(gam.16.poisson.num.red.bryo)
fam.gam.16.num.red.bryo
str(fam.gam.16.num.red.bryo)
ilink.gam.16.num.red.bryo<- fam.gam.16.num.red.bryo$linkinv
ilink.gam.16.num.red.bryo

mod.num.red.bryo<-gam.16.poisson.num.red.bryo
ndata.16.num.red.bryo <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                   length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.16.num.red.bryo <- add_column(ndata.16.num.red.bryo, fit = predict(mod.num.red.bryo, newdata = ndata.16.num.red.bryo, type = 'response'))

predict(mod.num.red.bryo, newdata = ndata.16.num.red.bryo, type = 'response')
ndata.16.num.red.bryo <- bind_cols(ndata.16.num.red.bryo, setNames(as_tibble(predict(mod.num.red.bryo, ndata.16.num.red.bryo, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))
## create the interval and backtransform
ndata.16.num.red.bryo <- mutate(ndata.16.num.red.bryo,
                           fit_resp  = ilink.gam.16.num.red.bryo(fit_link),
                           right_upr = ilink.gam.16.num.red.bryo(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.16.num.red.bryo(fit_link - (2 * se_link)))
ndata.16.num.red.bryo$min.10.pH.unscaled<-ndata.16.num.red.bryo$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.red.bryo.16 <- ggplot(ndata.16.num.red.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.red.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.red.bryo.16
ggsave("C:Graphs April 2020//num.red.bryo_pred.16.png")


# GAM negbin num.red.bryo / gam.8.nb.num.red.bryo --------------------------------------------------------------
nbinom.8.num.red.bryo <- fitdistr(invasion.exp.data.8_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.red.bryo, "nbinom", size = nbinom.8.num.red.bryo$estimate[[1]], mu = nbinom.8.num.red.bryo$estimate[[2]])
#getting theta

gam.8.nb.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.red.bryo, gam.8.nb.num.red.bryo.1, gam.8.poisson.num.red.bryo)


plot(gam.8.poisson.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.8.poisson.num.red.bryo)
#does not look good
qq_plot(gam.8.poisson.num.red.bryo, method = 'simulate')
k.check(gam.8.poisson.num.red.bryo)
summary(gam.8.poisson.num.red.bryo)

gam.8.poisson.num.red.bryo.unordered<- gam(num.red.bryo ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)
fam.gam.8.num.red.bryo <- family(gam.8.poisson.num.red.bryo)
ilink.gam.8.num.red.bryo<- fam.gam.8.num.red.bryo$linkinv


mod.num.red.bryo<-gam.8.poisson.num.red.bryo
ndata.8.num.red.bryo <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.8.num.red.bryo <- add_column(ndata.8.num.red.bryo, fit = predict(mod.num.red.bryo, newdata = ndata.8.num.red.bryo, type = 'response'))
predict(mod.num.red.bryo, newdata = ndata.8.num.red.bryo, type = 'response')
ndata.8.num.red.bryo <- bind_cols(ndata.8.num.red.bryo, setNames(as_tibble(predict(mod.num.red.bryo, ndata.8.num.red.bryo, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))
## create the interval and backtransform
ndata.8.num.red.bryo <- mutate(ndata.8.num.red.bryo,
                                fit_resp  = ilink.gam.8.num.red.bryo(fit_link),
                                right_upr = ilink.gam.8.num.red.bryo(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.8.num.red.bryo(fit_link - (2 * se_link)))
ndata.8.num.red.bryo$min.10.pH.unscaled<-ndata.8.num.red.bryo$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.red.bryo.8 <- ggplot(ndata.8.num.red.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.red.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,2.5)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.red.bryo.8
ggsave("C:Graphs April 2020//num.red.bryo_pred.8.png")


# GAM poisson num nudi / gam.16.poisson.num.nudi  ------------------------------------------------------------
nbinom.16.num.nudi <- fitdistr(invasion.exp.data.16_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.nudi, "nbinom", size = nbinom.16.num.nudi$estimate[[1]], mu = nbinom.16.num.nudi$estimate[[2]])
#theta

gam.16.nb.num.nudi.1<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.nb.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.16.poisson.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
AICtab(gam.16.nb.num.nudi, gam.16.nb.num.nudi.1, gam.16.poisson.num.nudi)

appraise(gam.16.poisson.num.nudi)
qq_plot(gam.16.poisson.num.nudi, method = 'simulate')
#looks quite good!
plot(gam.16.poisson.num.nudi, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.16.poisson.num.nudi)
summary(gam.16.poisson.num.nudi)
#resids a bit funny but same in neg bin

gam.16.poisson.num.nudi.unordered<- gam(num.nudi ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)
fam.gam.16.num.nudi <- family(gam.16.poisson.num.nudi)
ilink.gam.16.num.nudi<- fam.gam.16.num.nudi$linkinv
mod.num.nudi<-gam.16.poisson.num.nudi
ndata.16.num.nudi <- with(invasion.exp.data.16_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.16.num.nudi <- add_column(ndata.16.num.nudi, fit = predict(mod.num.nudi, newdata = ndata.16.num.nudi, type = 'response'))
predict(mod.num.nudi, newdata = ndata.16.num.nudi, type = 'response')
ndata.16.num.nudi <- bind_cols(ndata.16.num.nudi, setNames(as_tibble(predict(mod.num.nudi, ndata.16.num.nudi, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))
## create the interval and backtransform
ndata.16.num.nudi <- mutate(ndata.16.num.nudi,
                       fit_resp  = ilink.gam.16.num.nudi(fit_link),
                       right_upr = ilink.gam.16.num.nudi(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.16.num.nudi(fit_link - (2 * se_link)))
ndata.16.num.nudi$min.10.pH.unscaled<-ndata.16.num.nudi$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.nudi.16 <- ggplot(ndata.16.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.nudi, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.num.nudi.16
ggsave("C:Graphs April 2020//num.nudi_pred.16.png")


# GAM poisson num nudi / gam.8.poisson.num.nudi  ------------------------------------------------------------
nbinom.8.num.nudi <- fitdistr(invasion.exp.data.8_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.nudi, "nbinom", size = nbinom.8.num.nudi$estimate[[1]], mu = nbinom.8.num.nudi$estimate[[2]])
#theta

gam.8.nb.num.nudi.1<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.nb.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.16.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.8.poisson.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")
AICtab(gam.8.nb.num.nudi, gam.8.nb.num.nudi.1)

appraise(gam.8.nb.num.nudi.1)
qq_plot(gam.8.nb.num.nudi.1, method = 'simulate')
#looks quite good!
plot(gam.8.nb.num.nudi.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.8.nb.num.nudi.1)
summary(gam.8.nb.num.nudi.1)
#resids a bit funny but same in neg bin

gam.8.nb.num.nudi.1.unordered<- gam(num.nudi ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)
fam.gam.8.num.nudi <- family(gam.8.nb.num.nudi.1)
ilink.gam.8.num.nudi<- fam.gam.8.num.nudi$linkinv
mod.num.nudi<-gam.8.nb.num.nudi.1
ndata.8.num.nudi <- with(invasion.exp.data.8_zscores, 
                          data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.8.num.nudi <- add_column(ndata.8.num.nudi, fit = predict(mod.num.nudi, newdata = ndata.8.num.nudi, type = 'response'))
predict(mod.num.nudi, newdata = ndata.8.num.nudi, type = 'response')
ndata.8.num.nudi <- bind_cols(ndata.8.num.nudi, setNames(as_tibble(predict(mod.num.nudi, ndata.8.num.nudi, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))
## create the interval and backtransform
ndata.8.num.nudi <- mutate(ndata.8.num.nudi,
                            fit_resp  = ilink.gam.8.num.nudi(fit_link),
                            right_upr = ilink.gam.8.num.nudi(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.8.num.nudi(fit_link - (2 * se_link)))
ndata.8.num.nudi$min.10.pH.unscaled<-ndata.8.num.nudi$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.nudi.8 <- ggplot(ndata.8.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.nudi, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,4)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.num.nudi.8
ggsave("C:Graphs April 2020//num.nudi_pred.8.png")

# GAM nb() serpulids / gam.16.nb.num.serpulid.1 -----------------------------------------------------------
nbinom.16.num.serpulid <- fitdistr(invasion.exp.data.16_zscores$num.serpulid, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.serpulid, "nbinom", size = nbinom.16.num.serpulid$estimate[[1]], mu = nbinom.16.num.serpulid$estimate[[2]])
#theta

#negative binomial first
gam.16.nb.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.serpulid.1<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.serpulid, gam.16.nb.num.serpulid.1, gam.16.poisson.num.serpulid)
##gam.16.nb.num.serpulid.1 is best

appraise(gam.16.nb.num.serpulid.1)
qq_plot(gam.16.nb.num.serpulid.1, method = 'simulate')
plot(gam.16.nb.num.serpulid.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.16.nb.num.serpulid.1)
summary(gam.16.nb.num.serpulid.1)

gam.16.nb.num.serpulid.1.unordered<- gam(num.serpulid ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)
fam.gam.16.num.serpulid <- family(gam.16.nb.num.serpulid.1)
fam.gam.16.num.serpulid
ilink.gam.16.num.serpulid<- fam.gam.16.num.serpulid$linkinv
ilink.gam.16.num.serpulid
mod.num.serpulid<-gam.16.nb.num.serpulid.1
ndata.16.num.serpulid <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.16.num.serpulid <- add_column(ndata.16.num.serpulid, fit = predict(mod.num.serpulid, newdata = ndata.16.num.serpulid, type = 'response'))

predict(mod.num.serpulid, newdata = ndata.16.num.serpulid, type = 'response')
ndata.16.num.serpulid <- bind_cols(ndata.16.num.serpulid, setNames(as_tibble(predict(mod.num.serpulid, ndata.16.num.serpulid, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))
## create the interval and backtransform
ndata.16.num.serpulid <- mutate(ndata.16.num.serpulid,
                         fit_resp  = ilink.gam.16.num.serpulid(fit_link),
                         right_upr = ilink.gam.16.num.serpulid(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.num.serpulid(fit_link - (2 * se_link)))
ndata.16.num.serpulid$min.10.pH.unscaled<-ndata.16.num.serpulid$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.serpulid.16 <- ggplot(ndata.16.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.serpulid.16
ggsave("C:Graphs April 2020//num.serpulid_pred.16.png")

# GAM nb() serpulids / gam.8.nb.num.serpulid.1 -----------------------------------------------------------
nbinom.8.num.serpulid <- fitdistr(invasion.exp.data.8_zscores$num.serpulid, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.serpulid, "nbinom", size = nbinom.8.num.serpulid$estimate[[1]], mu = nbinom.8.num.serpulid$estimate[[2]])
#theta

#negative binomial first
gam.8.nb.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.serpulid.1<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.serpulid, gam.8.nb.num.serpulid.1, gam.8.poisson.num.serpulid)
##gam.8.poisson is best

appraise(gam.8.poisson.num.serpulid)
qq_plot(gam.8.poisson.num.serpulid, method = 'simulate')
plot(gam.8.poisson.num.serpulid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.8.poisson.num.serpulid)
summary(gam.8.poisson.num.serpulid)

gam.8.poisson.num.serpulid.unordered<- gam(num.serpulid ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")
want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)
fam.gam.8.num.serpulid <- family(gam.8.poisson.num.serpulid)
fam.gam.8.num.serpulid
ilink.gam.8.num.serpulid<- fam.gam.8.num.serpulid$linkinv
ilink.gam.8.num.serpulid
mod.num.serpulid<-gam.8.poisson.num.serpulid
ndata.8.num.serpulid <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.8.num.serpulid <- add_column(ndata.8.num.serpulid, fit = predict(mod.num.serpulid, newdata = ndata.8.num.serpulid, type = 'response'))

predict(mod.num.serpulid, newdata = ndata.8.num.serpulid, type = 'response')
ndata.8.num.serpulid <- bind_cols(ndata.8.num.serpulid, setNames(as_tibble(predict(mod.num.serpulid, ndata.8.num.serpulid, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))
## create the interval and backtransform
ndata.8.num.serpulid <- mutate(ndata.8.num.serpulid,
                                fit_resp  = ilink.gam.8.num.serpulid(fit_link),
                                right_upr = ilink.gam.8.num.serpulid(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.8.num.serpulid(fit_link - (2 * se_link)))
ndata.8.num.serpulid$min.10.pH.unscaled<-ndata.8.num.serpulid$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.serpulid.8 <- ggplot(ndata.8.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,6)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.serpulid.8
ggsave("C:Graphs April 2020//num.serpulid_pred.8.png")



# GAM negbin corella / gam.16.nb.num.corella -------------------------------------------------------------

nbinom.16.num.corella <- fitdistr(invasion.exp.data.16_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.corella, "nbinom", size = nbinom.16.num.corella$estimate[[1]], mu = nbinom.16.num.corella$estimate[[2]])
#extracting theta

gam.16.nb.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.corella.1<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.corella, gam.16.nb.num.corella.1, gam.16.poisson.num.corella)


appraise(gam.16.nb.num.corella)
#looks pretty good - slight pattern
qq_plot(gam.16.nb.num.corella, method = 'simulate')
plot(gam.16.nb.num.corella, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.16.nb.num.corella)
summary(gam.16.nb.num.corella)


gam.16.nb.num.corella.unordered<- gam(num.corella ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.corella$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)
fam.gam.16.num.corella <- family(gam.16.nb.num.corella)
ilink.gam.16.num.corella<- fam.gam.16.num.corella$linkinv


mod.num.corella<-gam.16.nb.num.corella
ndata.16.num.corella <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.16.num.corella <- add_column(ndata.16.num.corella, fit = predict(mod.num.corella, newdata = ndata.16.num.corella, type = 'response'))

predict(mod.num.corella, newdata = ndata.16.num.corella, type = 'response')
ndata.16.num.corella <- bind_cols(ndata.16.num.corella, setNames(as_tibble(predict(mod.num.corella, ndata.16.num.corella, se.fit = TRUE)[1:2]),
                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.num.corella <- mutate(ndata.16.num.corella,
                              fit_resp  = ilink.gam.16.num.corella(fit_link),
                              right_upr = ilink.gam.16.num.corella(fit_link + (2 * se_link)),
                              right_lwr = ilink.gam.16.num.corella(fit_link - (2 * se_link)))


ndata.16.num.corella$min.10.pH.unscaled<-ndata.16.num.corella$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

par(lheight=0.2) 
# plot 
plt.num.corella.16 <- ggplot(ndata.16.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.corella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.corella.16
ggsave("C:Graphs April 2020//num.corella_pred.16.png")

# GAM negbin corella / gam.8.nb.num.corella -------------------------------------------------------------

nbinom.8.num.corella <- fitdistr(invasion.exp.data.8_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.corella, "nbinom", size = nbinom.8.num.corella$estimate[[1]], mu = nbinom.8.num.corella$estimate[[2]])
#extracting theta

gam.8.nb.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.corella.1<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.corella, gam.8.nb.num.corella.1, gam.8.poisson.num.corella)
#poisson

appraise(gam.8.poisson.num.corella)
#looks pretty good - slight pattern
qq_plot(gam.8.poisson.num.corella, method = 'simulate')
plot(gam.8.poisson.num.corella, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.8.poisson.num.corella)
summary(gam.8.poisson.num.corella)

gam.8.poisson.num.corella.unordered<- gam(num.corella ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)
fam.gam.8.num.corella <- family(gam.8.poisson.num.corella)
ilink.gam.8.num.corella<- fam.gam.8.num.corella$linkinv
mod.num.corella<-gam.8.poisson.num.corella
ndata.8.num.corella <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.8.num.corella <- add_column(ndata.8.num.corella, fit = predict(mod.num.corella, newdata = ndata.8.num.corella, type = 'response'))

predict(mod.num.corella, newdata = ndata.8.num.corella, type = 'response')
ndata.8.num.corella <- bind_cols(ndata.8.num.corella, setNames(as_tibble(predict(mod.num.corella, ndata.8.num.corella, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.num.corella <- mutate(ndata.8.num.corella,
                               fit_resp  = ilink.gam.8.num.corella(fit_link),
                               right_upr = ilink.gam.8.num.corella(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.corella(fit_link - (2 * se_link)))


ndata.8.num.corella$min.10.pH.unscaled<-ndata.8.num.corella$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.corella.8 <- ggplot(ndata.8.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.corella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.corella.8
ggsave("C:Graphs April 2020//num.corella_pred.8.png")



# GAM poisson clam / gam.16.poisson.clam ----------------------------------------------------------------

nbinom.8.clam <- fitdistr(invasion.exp.data.16_zscores$clam, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$clam, "nbinom", size = nbinom.8.clam$estimate[[1]], mu = nbinom.8.clam$estimate[[2]])
#theta

#negative binomial first
gam.16.nb.clam<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.8.clam$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.clam.1<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.clam<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.clam, gam.16.nb.clam.1, gam.16.poisson.clam)

###gam poisson is best

#appraise(gam.16.poisson.clam)
#qq_plot(gam.16.poisson.clam, method = 'simulate')
plot(gam.16.poisson.clam, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k.check(gam.16.poisson.clam)
summary(gam.16.poisson.clam)

#residuals a bit patterny as well 
gam.16.poisson.clam.unordered<- gam(clam ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
summary(gam.16.poisson.clam.unordered)

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)

fam.gam.16.clam <- family(gam.16.poisson.clam)
fam.gam.16.clam
str(fam.gam.16.clam)
ilink.gam.16.clam<- fam.gam.16.clam$linkinv
ilink.gam.16.clam

mod.clam<-gam.16.poisson.clam
ndata.16.clam <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                    length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.16.clam <- add_column(ndata.16.clam, fit = predict(mod.clam, newdata = ndata.16.clam, type = 'response'))

predict(mod.clam, newdata = ndata.16.clam, type = 'response')
ndata.16.clam <- bind_cols(ndata.16.clam, setNames(as_tibble(predict(mod.clam, ndata.16.clam, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.clam <- mutate(ndata.16.clam,
                            fit_resp  = ilink.gam.16.clam(fit_link),
                            right_upr = ilink.gam.16.clam(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.16.clam(fit_link - (2 * se_link)))


ndata.16.clam$min.10.pH.unscaled<-ndata.16.clam$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 
plt.clam <- ggplot(ndata.16.clam, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = clam, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.clam
ggsave("C:Graphs April 2020//clam_pred.png")


# Fig 2 plot generation ---------------------------------------------------
library(patchwork)
fig.2<-wrap_plots(plt.gam.16.hydroid,plt.botryllid,plt.folliculina,plt.caprellid.percent,plt.membranipora,plt.didemnum,
          plt.mussel,plt.num.barn,plt.num.white.bryo,plt.num.red.bryo,plt.num.nudi,plt.num.serpulid,
          plt.orange_sponge,plt.num.corella,plt.clam, ncol=5)+
          plot_annotation(tag_levels = 'a')

fig.2

ggplot2::ggsave(plot=fig.2, "C:Data//For submission//For resubmission//RESUB2//First look//Fig2_resized.pdf", width=8.75, height=4.75, units="in")

#ggplot2::ggsave("C:Data//For submission//For resubmission//Fig2.png", width=65, height=35, units="cm")


# Pulling model results to a table ----------------------------------------

hydroid.gam<-summary(gam.16.beta.hydroid.3)
botryllus.gam<-summary(gam.16.beta.botryllid.3)
caprellid.gam<-summary(gam.16.beta.caprellid.percent)
folliculina.gam<-summary(gam.16.beta.folliculina)
membranipora.gam<-summary(gam.16.beta.membranipora)
didemnum.gam<-summary(gam.16.beta.didemnum)
mussel.gam<-summary(gam.16.nb.mussel)
alive.barn.gam<-summary(gam.16.nb.num.barn)
num.white.bryo.gam<-summary(gam.16.nb.num.white.bryo)
num.red.bryo.gam<-summary(gam.16.nb.num.red.bryo)
num.nudi.gam<-summary(gam.16.poisson.num.nudi)
num.serpulid.gam<-summary(gam.16.nb.num.serpulid.1)
clam.gam<-summary(gam.16.poisson.clam)
corella.gam<-summary(gam.16.nb.num.corella)
orange_sponge.gam<-summary(gam.16.nb.orange_sponge)

hydroid.gam.16.unordered<-summary(gam.16.beta.hydroid.3.unordered)
botryllus.gam.16.unordered<-summary(gam.16.beta.botryllid.3.unordered)
caprellid.gam.16.unordered<-summary(gam.16.beta.caprellid.percent.unordered)
folliculina.gam.16.unordered<-summary(gam.16.beta.folliculina.unordered)
membranipora.gam.16.unordered<-summary(gam.16.beta.membranipora.unordered)
didemnum.gam.16.unordered<-summary(gam.16.beta.didemnum.unordered)
mussel.gam.16.unordered<-summary(gam.16.nb.mussel.unordered)
alive.barn.gam.16.unordered<-summary(gam.16.nb.num.barn.unordered)
num.white.bryo.gam.16.unordered<-summary(gam.16.nb.num.white.bryo.unordered)
num.red.bryo.gam.16.unordered<-summary(gam.16.nb.num.red.bryo.unordered)
num.nudi.gam.16.unordered<-summary(gam.16.poisson.num.nudi.unordered)
num.serpulid.gam.16.unordered<-summary(gam.16.nb.num.serpulid.1.unordered)
clam.gam.16.unordered<-summary(gam.16.poisson.clam.unordered)
corella.gam.16.unordered<-summary(gam.16.nb.num.corella.unordered)
orange_sponge.gam.16.unordered<-summary(gam.16.nb.orange_sponge.unordered)

hydroid.gam.16.p.table<-as.data.frame(hydroid.gam.16.unordered$p.table)
hydroid.gam.16.s.table<-as.data.frame(hydroid.gam$s.table)

botryllus.gam.16.p.table<-as.data.frame(botryllus.gam.16.unordered$p.table)
botryllus.gam.16.s.table<-as.data.frame(botryllus.gam$s.table)

caprellid.gam.16.p.table<-as.data.frame(caprellid.gam.16.unordered$p.table)
caprellid.gam.16.s.table<-as.data.frame(caprellid.gam$s.table)

folliculina.gam.16.p.table<-as.data.frame(folliculina.gam.16.unordered$p.table)
folliculina.gam.16.s.table<-as.data.frame(folliculina.gam$s.table)

membranipora.gam.16.p.table<-as.data.frame(membranipora.gam.16.unordered$p.table)
membranipora.gam.16.s.table<-as.data.frame(membranipora.gam$s.table)

didemnum.gam.16.p.table<-as.data.frame(didemnum.gam.16.unordered$p.table)
didemnum.gam.16.s.table<-as.data.frame(didemnum.gam$s.table)

mussel.gam.16.p.table<-as.data.frame(mussel.gam.16.unordered$p.table)
mussel.gam.16.s.table<-as.data.frame(mussel.gam$s.table)

alive.barn.gam.16.p.table<-as.data.frame(alive.barn.gam.16.unordered$p.table)
alive.barn.gam.16.s.table<-as.data.frame(alive.barn.gam$s.table)

num.white.bryo.gam.16.p.table<-as.data.frame(num.white.bryo.gam.16.unordered$p.table)
num.white.bryo.gam.16.s.table<-as.data.frame(num.white.bryo.gam$s.table)

num.red.bryo.gam.16.p.table<-as.data.frame(num.red.bryo.gam.16.unordered$p.table)
num.red.bryo.gam.16.s.table<-as.data.frame(num.red.bryo.gam$s.table)

num.nudi.gam.16.p.table<-as.data.frame(num.nudi.gam.16.unordered$p.table)
num.nudi.gam.16.s.table<-as.data.frame(num.nudi.gam$s.table)

num.serpulid.gam.16.p.table<-as.data.frame(num.serpulid.gam.16.unordered$p.table)
num.serpulid.gam.16.s.table<-as.data.frame(num.serpulid.gam$s.table)

orange_sponge.gam.16.p.table<-as.data.frame(orange_sponge.gam.16.unordered$p.table)
orange_sponge.gam.16.s.table<-as.data.frame(orange_sponge.gam$s.table)

corella.gam.16.p.table<-as.data.frame(corella.gam.16.unordered$p.table)
corella.gam.16.s.table<-as.data.frame(corella.gam$s.table)

clam.gam.16.p.table<-as.data.frame(clam.gam.16.unordered$p.table)
clam.gam.16.s.table<-as.data.frame(clam.gam$s.table)



head(clam.gam.16.p.table)
#### Building the stats table
ptable<-rbind(hydroid.gam.16.p.table, 
               botryllus.gam.16.p.table, 
               caprellid.gam.16.p.table,
               folliculina.gam.16.p.table,
               membranipora.gam.16.p.table,
               didemnum.gam.16.p.table,
               mussel.gam.16.p.table,
               alive.barn.gam.16.p.table,
               num.white.bryo.gam.16.p.table,
               num.red.bryo.gam.16.p.table,
               num.nudi.gam.16.p.table,
               num.serpulid.gam.16.p.table,
               orange_sponge.gam.16.p.table,
               corella.gam.16.p.table,
               clam.gam.16.p.table)


colnames(ptable) <- c("Estimate", "SE", "z", "p")
ptable$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Caprella",7, 9) %>% 
  group_rows("Folliculina", 10, 12) %>% 
  group_rows("Membranipora", 13,15) %>% 
  group_rows("Didemnum", 16, 18) %>% 
  group_rows("Mussels", 19, 21) %>% 
  group_rows("Barnacles", 22, 24) %>% 
  group_rows("num.white.bryo", 25, 27) %>% 
  group_rows("num.red.bryoporella", 28, 30) %>% 
  group_rows("Hermissenda", 31, 33) %>% 
  group_rows("Serpulid", 34, 36) %>% 
  group_rows("Sponge", 37, 39) %>% 
  group_rows("Corella", 40, 42) %>% 
  group_rows("Clams", 43, 45) %>% 
save_kable(file = "C:Data//For submission//ptable.html", self_contained = T)


### s table
stable<-rbind(hydroid.gam.16.s.table, 
              botryllus.gam.16.s.table, 
              caprellid.gam.16.s.table,
              folliculina.gam.16.s.table,
              membranipora.gam.16.s.table,
              didemnum.gam.16.s.table,
              mussel.gam.16.s.table,
              alive.barn.gam.16.s.table,
              num.white.bryo.gam.16.s.table,
              num.red.bryo.gam.16.s.table,
              num.nudi.gam.16.s.table,
              num.serpulid.gam.16.s.table,
              orange_sponge.gam.16.s.table,
              corella.gam.16.s.table,
              clam.gam.16.s.table)


colnames(stable) <- c("Estimated_df", "Reference_df", "Chi_squared", "p_smooth")
stable$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Low quality food", "smooth min.10.pH * High quality food"))

#stable$species<-rep(c("hydroid", 
#                      "botryllus", 
#                       "caprellid",
#                     "folliculina",
#                     "membranipora",
#                    "didemnum",
  #                    "mussel",
  #                   "alive.barn",
  #                   "num.white.bryo",
  #                   "num.red.bryo",
  #                   "num.nudi",
  #                   "num.serpulid",
  #                   "orange_sponge",
  #                   "corella",
  #                   "clam"), each=3)


stable %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Caprella",7, 9) %>% 
  group_rows("Folliculina", 10, 12) %>% 
  group_rows("Membranipora", 13,15) %>% 
  group_rows("Didemnum", 16, 18) %>% 
  group_rows("Mussels", 19, 21) %>% 
  group_rows("Barnacles", 22, 24) %>% 
  group_rows("num.white.bryo", 25, 27) %>% 
  group_rows("num.red.bryoporella", 28, 30) %>% 
  group_rows("Hermissenda", 31, 33) %>% 
  group_rows("Serpulid", 34, 36) %>% 
  group_rows("Sponge", 37, 39) %>% 
  group_rows("Corella", 40, 42) %>% 
  group_rows("Clams", 43, 45) %>% 
  save_kable(file = "C:Data//For submission//stable.html", self_contained = T)

  
pstable<-cbind(ptable, stable)

pstable %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth, Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia, beta", 1, 3) %>%
  group_rows("Botryllus, beta", 4, 6) %>% 
  group_rows("Caprella, beta",7, 9) %>% 
  group_rows("Folliculina, beta", 10, 12) %>% 
  group_rows("Membranipora, beta", 13,15) %>% 
  group_rows("Didemnum, beta", 16, 18) %>% 
  group_rows("Mussels, negative binomial", 19, 21) %>% 
  group_rows("Barnacles, negative binomial", 22, 24) %>% 
  group_rows("num.white.bryo, negative binomial", 25, 27) %>% 
  group_rows("num.red.bryoporella, negative binomial", 28, 30) %>% 
  group_rows("Hermissenda, negative binomial", 31, 33) %>% 
  group_rows("Serpulid, negative binomial", 34, 36) %>% 
  group_rows("Sponge, negative binomial", 37, 39) %>% 
  group_rows("Corella, negative binomial", 40, 42) %>% 
  group_rows("Clams, poisson", 43, 45) %>% 
  save_kable(file = "C:Data//For submission//For resubmission//RESUB2//First look//pstable.html", self_contained = T)


# Richness ----------------------------------------------------------------
gam.16.nb.richness.1<- gam(richness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.richness<- gam(richness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.richness.1, gam.16.poisson.richness)

plot(gam.16.poisson.richness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.poisson.richness)
#okay but qq plot not the best on ends
#qq_plot(gam.16.poisson.richness, method = 'simulate')
#k.check(gam.16.poisson.richness)
summary(gam.16.poisson.richness)
#a few outside the QQ plot on both ends
#low p value for k - but NS and edf is not super close to k-index

gam.16.poisson.richness.unordered<- gam(richness ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
fam.gam.16.richness <- family(gam.16.poisson.richness)
ilink.gam.16.richness<- fam.gam.16.richness$linkinv
ilink.gam.16.richness


mod.richness<-gam.16.poisson.richness
ndata.16.richness <- with(invasion.exp.data.16_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.richness <- add_column(ndata.16.richness, fit = predict(mod.richness, newdata = ndata.16.richness, type = 'response'))

predict(mod.richness, newdata = ndata.16.richness, type = 'response')
ndata.16.richness <- bind_cols(ndata.16.richness, setNames(as_tibble(predict(mod.richness, ndata.16.richness, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.richness <- mutate(ndata.16.richness,
                               fit_resp  = ilink.gam.16.richness(fit_link),
                               right_upr = ilink.gam.16.richness(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.16.richness(fit_link - (2 * se_link)))


ndata.16.richness$min.10.pH.unscaled<-ndata.16.richness$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.richness <- ggplot(ndata.16.richness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = richness, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species richness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.richness,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,20)
plt.richness
ggsave("C:Graphs April 2020//richness_pred.png")







# Evenness ----------------------------------------------------------------

gam.16.lm.evenness<- gam(evenness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.evenness.1<- gam(evenness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

AICtab( gam.16.lm.evenness, gam.16.gamma.evenness.1)

plot(gam.16.lm.evenness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.evenness)
#looks very good
#qq_plot(gam.16.lm.evenness, method = 'simulate')
#k.check(gam.16.lm.evenness)
summary(gam.16.lm.evenness)

gam.16.lm.evenness.unordered<- gam(evenness ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")

fam.gam.16.evenness <- family(gam.16.lm.evenness)
ilink.gam.16.evenness<- fam.gam.16.evenness$linkinv
mod.evenness<-gam.16.lm.evenness
ndata.16.evenness <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.evenness <- add_column(ndata.16.evenness, fit = predict(mod.evenness, newdata = ndata.16.evenness, type = 'response'))

predict(mod.evenness, newdata = ndata.16.evenness, type = 'response')
ndata.16.evenness <- bind_cols(ndata.16.evenness, setNames(as_tibble(predict(mod.evenness, ndata.16.evenness, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.evenness <- mutate(ndata.16.evenness,
                         fit_resp  = ilink.gam.16.evenness(fit_link),
                         right_upr = ilink.gam.16.evenness(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.evenness(fit_link - (2 * se_link)))


ndata.16.evenness$min.10.pH.unscaled<-ndata.16.evenness$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.evenness <- ggplot(ndata.16.evenness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = evenness, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species evenness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.evenness
ggsave("C:Graphs April 2020//evenness_pred.png")


# Occupied space ----------------------------------------------------------

gam.16.lm.occupied.space<- gam(occupied.space ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.occupied.space.1<- gam(occupied.space*0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.binomial.occupied.space<- gam(formula = cbind(occupied.space, 100-occupied.space)~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

gam.16.beta.occupied.space<- gam(occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.occupied.space.1<- gam(occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.occupied.space.2<- gam(occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.occupied.space.3<- gam(occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.16.lm.occupied.space, gam.16.beta.occupied.space.3, gam.16.gamma.occupied.space.1, gam.16.beta.occupied.space, gam.16.beta.occupied.space.1, gam.16.beta.occupied.space.2, gam.16.binomial.occupied.space)
#beta cauchit is best 

plot(gam.16.beta.occupied.space.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.occupied.space.3)
#a bit blocky
#qq_plot(gam.16.beta.occupied.space.3, method = 'simulate')
#k.check(gam.16.beta.occupied.space.3)
#k gettinga bit low but ns
summary(gam.16.beta.occupied.space.3)
gam.16.beta.occupied.space.3.unordered<- gam(occupied.space.001~ s(min.10.pH, k=15)+ Invasives + s(min.10.pH, by=oInvasives, k=15), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.16.beta.occupied.space.3.unordered)

fam.gam.16.occupied.space <- family(gam.16.beta.occupied.space.3)
fam.gam.16.occupied.space
str(fam.gam.16.occupied.space)
ilink.gam.16.occupied.space<- fam.gam.16.occupied.space$linkinv
ilink.gam.16.occupied.space


mod.occupied.space<-gam.16.beta.occupied.space.3
ndata.16.occupied.space <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.occupied.space <- add_column(ndata.16.occupied.space, fit = predict(mod.occupied.space, newdata = ndata.16.occupied.space, type = 'response'))

predict(mod.occupied.space, newdata = ndata.16.occupied.space, type = 'response')
ndata.16.occupied.space <- bind_cols(ndata.16.occupied.space, setNames(as_tibble(predict(mod.occupied.space, ndata.16.occupied.space, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.occupied.space <- mutate(ndata.16.occupied.space,
                         fit_resp  = ilink.gam.16.occupied.space(fit_link),
                         right_upr = ilink.gam.16.occupied.space(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.occupied.space(fit_link - (2 * se_link)))


ndata.16.occupied.space$min.10.pH.unscaled<-ndata.16.occupied.space$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.occupied.space <- ggplot(ndata.16.occupied.space, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = occupied.space.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.occupied.space
ggsave("C:Graphs April 2020//occupied.space_pred.png")



# Everything wet weight ---------------------------------------------------

gam.16.lm.everything.wet.weight<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.everything.wet.weight<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.everything.wet.weight.1, gam.16.lm.log.everything.wet.weight, gam.16.tweedie.everything.wet.weight.1,  gam.16.lm.everything.wet.weight, gam.16.gamma.everything.wet.weight.1)

#gam.16.lm.log.everything.wet.weight is best 


plot(gam.16.lm.log.everything.wet.weight, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.log.everything.wet.weight)
#Looks quite good
#qq_plot(gam.16.lm.log.everything.wet.weight, method = 'simulate')
#k.check(gam.16.lm.log.everything.wet.weight)
summary(gam.16.lm.log.everything.wet.weight)

gam.16.lm.log.everything.wet.weight.unordered<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")

fam.gam.16.everything.wet.weight <- family(gam.16.lm.log.everything.wet.weight)
fam.gam.16.everything.wet.weight
ilink.gam.16.everything.wet.weight<- fam.gam.16.everything.wet.weight$linkinv
ilink.gam.16.everything.wet.weight


mod.everything.wet.weight<-gam.16.lm.log.everything.wet.weight
ndata.16.everything.wet.weight <- with(invasion.exp.data.16_zscores, 
                                    data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                    length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.everything.wet.weight <- add_column(ndata.16.everything.wet.weight, fit = predict(mod.everything.wet.weight, newdata = ndata.16.everything.wet.weight, type = 'response'))

predict(mod.everything.wet.weight, newdata = ndata.16.everything.wet.weight, type = 'response')
ndata.16.everything.wet.weight <- bind_cols(ndata.16.everything.wet.weight, setNames(as_tibble(predict(mod.everything.wet.weight, ndata.16.everything.wet.weight, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.everything.wet.weight <- mutate(ndata.16.everything.wet.weight,
                         fit_resp  = ilink.gam.16.everything.wet.weight(fit_link),
                         right_upr = ilink.gam.16.everything.wet.weight(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.everything.wet.weight(fit_link - (2 * se_link)))


ndata.16.everything.wet.weight$min.10.pH.unscaled<-ndata.16.everything.wet.weight$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight <- ggplot(ndata.16.everything.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(everything.wet.weight), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total wet biomass per mesocosm\n(g, log scale)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.everything.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight
ggsave("C:Graphs April 2020//everything.wet.weight_pred.png")


# Everything wet weight per 1 ---------------------------------------------

gam.16.lm.everything.wet.weight.per.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.everything.wet.weight.per.1<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.everything.wet.weight.per.1.1, gam.16.lm.log.everything.wet.weight.per.1, gam.16.tweedie.everything.wet.weight.per.1.1,gam.16.lm.everything.wet.weight.per.1, gam.16.gamma.everything.wet.weight.per.1.1)

#gam.16.lm.log.everything.wet.weight.per.1 is best by far


plot(gam.16.lm.log.everything.wet.weight.per.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.log.everything.wet.weight.per.1)
#good but mabe a bit funnelly
#qq_plot(gam.16.lm.log.everything.wet.weight.per.1, method = 'simulate')
#k.check(gam.16.lm.log.everything.wet.weight.per.1)
summary(gam.16.lm.log.everything.wet.weight.per.1)

gam.16.lm.log.everything.wet.weight.per.1.unordered<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")

fam.gam.16.everything.wet.weight.per.1 <- family(gam.16.lm.log.everything.wet.weight.per.1)
fam.gam.16.everything.wet.weight.per.1
ilink.gam.16.everything.wet.weight.per.1<- fam.gam.16.everything.wet.weight.per.1$linkinv
ilink.gam.16.everything.wet.weight.per.1


mod.everything.wet.weight.per.1<-gam.16.lm.log.everything.wet.weight.per.1
ndata.16.everything.wet.weight.per.1 <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.everything.wet.weight.per.1 <- add_column(ndata.16.everything.wet.weight.per.1, fit = predict(mod.everything.wet.weight.per.1, newdata = ndata.16.everything.wet.weight.per.1, type = 'response'))

predict(mod.everything.wet.weight.per.1, newdata = ndata.16.everything.wet.weight.per.1, type = 'response')
ndata.16.everything.wet.weight.per.1 <- bind_cols(ndata.16.everything.wet.weight.per.1, setNames(as_tibble(predict(mod.everything.wet.weight.per.1, ndata.16.everything.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.everything.wet.weight.per.1 <- mutate(ndata.16.everything.wet.weight.per.1,
                                      fit_resp  = ilink.gam.16.everything.wet.weight.per.1(fit_link),
                                      right_upr = ilink.gam.16.everything.wet.weight.per.1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.everything.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.16.everything.wet.weight.per.1$min.10.pH.unscaled<-ndata.16.everything.wet.weight.per.1$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight.per.1 <- ggplot(ndata.16.everything.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(everything.wet.weight.per.1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Log Biomass per 1 % cover \n(wet weight)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.everything.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight.per.1
ggsave("C:Graphs April 2020//everything.wet.weight.per.1_pred.png")




# Dry weight --------------------------------------------------------------
gam.16.lm.total_dry_biomass<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.total_dry_biomass<- gam(log(total_dry_biomass) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.total_dry_biomass.1, gam.16.lm.log.total_dry_biomass, gam.16.tweedie.total_dry_biomass.1, gam.16.lm.total_dry_biomass, gam.16.gamma.total_dry_biomass.1)

#gam.16.lm.loglink is best but lm is only 0.1 so going with simplest.... 


plot(gam.16.lm.total_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.total_dry_biomass )
#Looks v good
#qq_plot(gam.16.lm.total_dry_biomass , method = 'simulate')
#k.check(gam.16.lm.total_dry_biomass )
summary(gam.16.lm.total_dry_biomass )

gam.16.lm.total_dry_biomass.unordered<- gam(total_dry_biomass ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
summary(gam.16.lm.total_dry_biomass.unordered)

fam.gam.16.total_dry_biomass <- family(gam.16.lm.total_dry_biomass)
fam.gam.16.total_dry_biomass
str(fam.gam.16.total_dry_biomass)
ilink.gam.16.total_dry_biomass<- fam.gam.16.total_dry_biomass$linkinv
ilink.gam.16.total_dry_biomass


mod.total_dry_biomass<-gam.16.lm.total_dry_biomass
ndata.16.total_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.total_dry_biomass <- add_column(ndata.16.total_dry_biomass, fit = predict(mod.total_dry_biomass, newdata = ndata.16.total_dry_biomass, type = 'response'))

predict(mod.total_dry_biomass, newdata = ndata.16.total_dry_biomass, type = 'response')
ndata.16.total_dry_biomass <- bind_cols(ndata.16.total_dry_biomass, setNames(as_tibble(predict(mod.total_dry_biomass, ndata.16.total_dry_biomass, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.total_dry_biomass <- mutate(ndata.16.total_dry_biomass,
                                      fit_resp  = ilink.gam.16.total_dry_biomass(fit_link),
                                      right_upr = ilink.gam.16.total_dry_biomass(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.total_dry_biomass(fit_link - (2 * se_link)))


ndata.16.total_dry_biomass$min.10.pH.unscaled<-ndata.16.total_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass <- ggplot(ndata.16.total_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (total_dry_biomass), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.total_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//total_dry_biomass_pred.png")


# Dry weight per 1 % cover --------------------------------------------------------------

gam.16.lm.total_dry_biomass_per1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.total_dry_biomass_per1.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.total_dry_biomass_per1<- gam(log(total_dry_biomass_per1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.total_dry_biomass_per1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.total_dry_biomass_per1.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.total_dry_biomass_per1.1, gam.16.lm.log.total_dry_biomass_per1, gam.16.tweedie.total_dry_biomass_per1, gam.16.lm.total_dry_biomass_per1, gam.16.gamma.total_dry_biomass_per1.1)

#tweedie is best by 1.7

plot(gam.16.tweedie.total_dry_biomass_per1 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.tweedie.total_dry_biomass_per1 )
#not available for tweedie
#qq_plot(gam.16.tweedie.total_dry_biomass_per1 , method = 'simulate')
#looks good
#k.check(gam.16.tweedie.total_dry_biomass_per1 )
summary(gam.16.tweedie.total_dry_biomass_per1 )

gam.16.tweedie.total_dry_biomass_per1.unordered<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=Invasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
summary(gam.16.tweedie.total_dry_biomass_per1.unordered)

fam.gam.16.total_dry_biomass_per1 <- family(gam.16.tweedie.total_dry_biomass_per1)
fam.gam.16.total_dry_biomass_per1
str(fam.gam.16.total_dry_biomass_per1)
ilink.gam.16.total_dry_biomass_per1<- fam.gam.16.total_dry_biomass_per1$linkinv
ilink.gam.16.total_dry_biomass_per1


mod.total_dry_biomass_per1<-gam.16.tweedie.total_dry_biomass_per1
ndata.16.total_dry_biomass_per1 <- with(invasion.exp.data.16_zscores, 
                                     data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.total_dry_biomass_per1 <- add_column(ndata.16.total_dry_biomass_per1, fit = predict(mod.total_dry_biomass_per1, newdata = ndata.16.total_dry_biomass_per1, type = 'response'))

predict(mod.total_dry_biomass_per1, newdata = ndata.16.total_dry_biomass_per1, type = 'response')
ndata.16.total_dry_biomass_per1 <- bind_cols(ndata.16.total_dry_biomass_per1, setNames(as_tibble(predict(mod.total_dry_biomass_per1, ndata.16.total_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.total_dry_biomass_per1 <- mutate(ndata.16.total_dry_biomass_per1,
                                  fit_resp  = ilink.gam.16.total_dry_biomass_per1(fit_link),
                                  right_upr = ilink.gam.16.total_dry_biomass_per1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.16.total_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.16.total_dry_biomass_per1$min.10.pH.unscaled<-ndata.16.total_dry_biomass_per1$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass_per1 <- ggplot(ndata.16.total_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (total_dry_biomass_per1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.total_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass_per1
ggplot2::ggsave("C:Graphs April 2020//total_dry_biomass_per1_pred.png")

# hydroid biomass ---------------------------------------------------------

gam.16.lm.hydroid_dry_biomass<- gam(hydroid_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.hydroid_dry_biomass<- gam(log(hydroid_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.16.loglink.hydroid_dry_biomass, gam.16.lm.log.hydroid_dry_biomass, gam.16.tweedie.hydroid_dry_biomass,  gam.16.lm.hydroid_dry_biomass, gam.16.gamma.hydroid_dry_biomass)

#gamma is best 
plot(gam.16.gamma.hydroid_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.gamma.hydroid_dry_biomass )
#looks good, maybe slightly funnelly
#qq_plot(gam.16.gamma.hydroid_dry_biomass , method = 'simulate')
#k.check(gam.16.gamma.hydroid_dry_biomass )
summary(gam.16.gamma.hydroid_dry_biomass )

gam.16.gamma.hydroid_dry_biomass.unordered<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.16.hydroid_dry_biomass <- family(gam.16.gamma.hydroid_dry_biomass)
fam.gam.16.hydroid_dry_biomass
str(fam.gam.16.hydroid_dry_biomass)
ilink.gam.16.hydroid_dry_biomass<- fam.gam.16.hydroid_dry_biomass$linkinv
ilink.gam.16.hydroid_dry_biomass


mod.hydroid_dry_biomass<-gam.16.gamma.hydroid_dry_biomass
ndata.16.hydroid_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.hydroid_dry_biomass <- add_column(ndata.16.hydroid_dry_biomass, fit = predict(mod.hydroid_dry_biomass, newdata = ndata.16.hydroid_dry_biomass, type = 'response'))

predict(mod.hydroid_dry_biomass, newdata = ndata.16.hydroid_dry_biomass, type = 'response')
ndata.16.hydroid_dry_biomass <- bind_cols(ndata.16.hydroid_dry_biomass, setNames(as_tibble(predict(mod.hydroid_dry_biomass, ndata.16.hydroid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.hydroid_dry_biomass <- mutate(ndata.16.hydroid_dry_biomass,
                                  fit_resp  = ilink.gam.16.hydroid_dry_biomass(fit_link),
                                  right_upr = ilink.gam.16.hydroid_dry_biomass(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.16.hydroid_dry_biomass(fit_link - (2 * se_link)))


ndata.16.hydroid_dry_biomass$min.10.pH.unscaled<-ndata.16.hydroid_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.hydroid_dry_biomass <- ggplot(ndata.16.hydroid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (hydroid_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Obelia") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.hydroid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.hydroid_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//hydroid_dry_biomass_pred.png")




# botryllus biomass -------------------------------------------------------
invasion.exp.data.16_zscores$tunicate_dry_biomass[invasion.exp.data.16_zscores$tunicate_dry_biomass<0]<-0

gam.16.lm.tunicate_dry_biomass<- gam(tunicate_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.tunicate_dry_biomass<- gam(log(tunicate_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.tunicate_dry_biomass, gam.16.lm.log.tunicate_dry_biomass, gam.16.tweedie.tunicate_dry_biomass,gam.16.lm.tunicate_dry_biomass, gam.16.gamma.tunicate_dry_biomass)

#gamma is the best

plot(gam.16.gamma.tunicate_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.gamma.tunicate_dry_biomass )
#a bit patterny
#qq_plot(gam.16.gamma.tunicate_dry_biomass , method = 'simulate')
#k.check(gam.16.gamma.tunicate_dry_biomass )
summary(gam.16.gamma.tunicate_dry_biomass )

gam.16.gamma.tunicate_dry_biomass.unordered<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.16.tunicate_dry_biomass <- family(gam.16.gamma.tunicate_dry_biomass)
fam.gam.16.tunicate_dry_biomass
str(fam.gam.16.tunicate_dry_biomass)
ilink.gam.16.tunicate_dry_biomass<- fam.gam.16.tunicate_dry_biomass$linkinv
ilink.gam.16.tunicate_dry_biomass


mod.tunicate_dry_biomass<-gam.16.gamma.tunicate_dry_biomass
ndata.16.tunicate_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.tunicate_dry_biomass <- add_column(ndata.16.tunicate_dry_biomass, fit = predict(mod.tunicate_dry_biomass, newdata = ndata.16.tunicate_dry_biomass, type = 'response'))

predict(mod.tunicate_dry_biomass, newdata = ndata.16.tunicate_dry_biomass, type = 'response')
ndata.16.tunicate_dry_biomass <- bind_cols(ndata.16.tunicate_dry_biomass, setNames(as_tibble(predict(mod.tunicate_dry_biomass, ndata.16.tunicate_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.tunicate_dry_biomass <- mutate(ndata.16.tunicate_dry_biomass,
                                    fit_resp  = ilink.gam.16.tunicate_dry_biomass(fit_link),
                                    right_upr = ilink.gam.16.tunicate_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.16.tunicate_dry_biomass(fit_link - (2 * se_link)))


ndata.16.tunicate_dry_biomass$min.10.pH.unscaled<-ndata.16.tunicate_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.tunicate_dry_biomass <- ggplot(ndata.16.tunicate_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (tunicate_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.tunicate_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.tunicate_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//tunicate_dry_biomass_pred.png")



# caprellid biomass -------------------------------------------------------
gam.16.lm.caprellid_dry_biomass<- gam(caprellid_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.caprellid_dry_biomass<- gam(log(caprellid_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.caprellid_dry_biomass, gam.16.lm.log.caprellid_dry_biomass, gam.16.tweedie.caprellid_dry_biomass, gam.16.lm.caprellid_dry_biomass, gam.16.gamma.caprellid_dry_biomass)
#gamma by 1.2

plot(gam.16.gamma.caprellid_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.gamma.caprellid_dry_biomass )
# look good
#qq_plot(gam.16.gamma.caprellid_dry_biomass , method = 'simulate')
#k.check(gam.16.gamma.caprellid_dry_biomass )
summary(gam.16.gamma.caprellid_dry_biomass )
gam.16.gamma.caprellid_dry_biomass.unordered<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.16.caprellid_dry_biomass <- family(gam.16.gamma.caprellid_dry_biomass)
ilink.gam.16.caprellid_dry_biomass<- fam.gam.16.caprellid_dry_biomass$linkinv
mod.caprellid_dry_biomass<-gam.16.gamma.caprellid_dry_biomass
ndata.16.caprellid_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.caprellid_dry_biomass <- add_column(ndata.16.caprellid_dry_biomass, fit = predict(mod.caprellid_dry_biomass, newdata = ndata.16.caprellid_dry_biomass, type = 'response'))

predict(mod.caprellid_dry_biomass, newdata = ndata.16.caprellid_dry_biomass, type = 'response')
ndata.16.caprellid_dry_biomass <- bind_cols(ndata.16.caprellid_dry_biomass, setNames(as_tibble(predict(mod.caprellid_dry_biomass, ndata.16.caprellid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.caprellid_dry_biomass <- mutate(ndata.16.caprellid_dry_biomass,
                                    fit_resp  = ilink.gam.16.caprellid_dry_biomass(fit_link),
                                    right_upr = ilink.gam.16.caprellid_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.16.caprellid_dry_biomass(fit_link - (2 * se_link)))


ndata.16.caprellid_dry_biomass$min.10.pH.unscaled<-ndata.16.caprellid_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass <- ggplot(ndata.16.caprellid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (caprellid_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.caprellid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//caprellid_dry_biomass_pred.png")



# caprellid biomass per invididual -------------------------------------------------------

invasion.exp.data.16_zscores$caprellid_dry_biomass_per1<-invasion.exp.data.16_zscores$caprellid_dry_biomass/(food.caprellid.data_zscores$total.caprellids)

gam.16.lm.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.caprellid_dry_biomass_per1<- gam(log(caprellid_dry_biomass_per1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.caprellid_dry_biomass_per1, gam.16.lm.log.caprellid_dry_biomass_per1, gam.16.tweedie.caprellid_dry_biomass_per1,  gam.16.lm.caprellid_dry_biomass_per1, gam.16.gamma.caprellid_dry_biomass_per1)

#gamma is the best by 2.2
plot(gam.16.gamma.caprellid_dry_biomass_per1 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.gamma.caprellid_dry_biomass_per1 )
#a bit funnelly and qq right tail
#qq_plot(gam.16.gamma.caprellid_dry_biomass_per1 , method = 'simulate')
#k.check(gam.16.gamma.caprellid_dry_biomass_per1 )
summary(gam.16.gamma.caprellid_dry_biomass_per1 )
gam.16.gamma.caprellid_dry_biomass_per1.unordered<- gam(caprellid_dry_biomass_per1+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.16.caprellid_dry_biomass_per1 <- family(gam.16.gamma.caprellid_dry_biomass_per1)
fam.gam.16.caprellid_dry_biomass_per1
str(fam.gam.16.caprellid_dry_biomass_per1)
ilink.gam.16.caprellid_dry_biomass_per1<- fam.gam.16.caprellid_dry_biomass_per1$linkinv
ilink.gam.16.caprellid_dry_biomass_per1


mod.caprellid_dry_biomass_per1<-gam.16.gamma.caprellid_dry_biomass_per1
ndata.16.caprellid_dry_biomass_per1 <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.caprellid_dry_biomass_per1 <- add_column(ndata.16.caprellid_dry_biomass_per1, fit = predict(mod.caprellid_dry_biomass_per1, newdata = ndata.16.caprellid_dry_biomass_per1, type = 'response'))

predict(mod.caprellid_dry_biomass_per1, newdata = ndata.16.caprellid_dry_biomass_per1, type = 'response')
ndata.16.caprellid_dry_biomass_per1 <- bind_cols(ndata.16.caprellid_dry_biomass_per1, setNames(as_tibble(predict(mod.caprellid_dry_biomass_per1, ndata.16.caprellid_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.caprellid_dry_biomass_per1 <- mutate(ndata.16.caprellid_dry_biomass_per1,
                                      fit_resp  = ilink.gam.16.caprellid_dry_biomass_per1(fit_link),
                                      right_upr = ilink.gam.16.caprellid_dry_biomass_per1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.caprellid_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.16.caprellid_dry_biomass_per1$min.10.pH.unscaled<-ndata.16.caprellid_dry_biomass_per1$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass_per1 <- ggplot(ndata.16.caprellid_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (caprellid_dry_biomass_per1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.caprellid_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,0.012)
plt.caprellid_dry_biomass_per1
ggplot2::ggsave("C:Graphs April 2020//caprellid_dry_biomass_per1_pred.png")


# rest biomass ------------------------------------------------------------

#kcheck was significant so increase k from 10 to 11
gam.16.lm.rest_dry_biomass<- gam(rest_dry_biomass ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.rest_dry_biomass<- gam(log(rest_dry_biomass+0.1) ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.rest_dry_biomass, gam.16.lm.log.rest_dry_biomass, gam.16.tweedie.rest_dry_biomass, gam.16.lm.rest_dry_biomass, gam.16.gamma.rest_dry_biomass)

#tweedie the best 


plot(gam.16.tweedie.rest_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#qq_plot(gam.16.tweedie.rest_dry_biomass , method = 'simulate')
#k.check(gam.16.tweedie.rest_dry_biomass )
summary(gam.16.tweedie.rest_dry_biomass)

gam.16.tweedie.rest_dry_biomass.unordered<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
summary(gam.16.tweedie.rest_dry_biomass.unordered)

fam.gam.16.rest_dry_biomass <- family(gam.16.tweedie.rest_dry_biomass)
fam.gam.16.rest_dry_biomass
str(fam.gam.16.rest_dry_biomass)
ilink.gam.16.rest_dry_biomass<- fam.gam.16.rest_dry_biomass$linkinv
ilink.gam.16.rest_dry_biomass


mod.rest_dry_biomass<-gam.16.tweedie.rest_dry_biomass
ndata.16.rest_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.16.rest_dry_biomass <- add_column(ndata.16.rest_dry_biomass, fit = predict(mod.rest_dry_biomass, newdata = ndata.16.rest_dry_biomass, type = 'response'))

predict(mod.rest_dry_biomass, newdata = ndata.16.rest_dry_biomass, type = 'response')
ndata.16.rest_dry_biomass <- bind_cols(ndata.16.rest_dry_biomass, setNames(as_tibble(predict(mod.rest_dry_biomass, ndata.16.rest_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.rest_dry_biomass <- mutate(ndata.16.rest_dry_biomass,
                                    fit_resp  = ilink.gam.16.rest_dry_biomass(fit_link),
                                    right_upr = ilink.gam.16.rest_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.16.rest_dry_biomass(fit_link - (2 * se_link)))


ndata.16.rest_dry_biomass$min.10.pH.unscaled<-ndata.16.rest_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.rest_dry_biomass <- ggplot(ndata.16.rest_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (rest_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression("Remaining dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.rest_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.rest_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//rest_dry_biomass_pred.png")




# Mussel wet weight -------------------------------------------------------

gam.16.lm.Mussel.wet.weight<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.Mussel.wet.weight.1<- gam(Mussel.wet.weight+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.Mussel.wet.weight<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.Mussel.wet.weight.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.Mussel.wet.weight.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.16.loglink.Mussel.wet.weight.1, gam.16.lm.log.Mussel.wet.weight, gam.16.tweedie.Mussel.wet.weight.1,gam.16.lm.Mussel.wet.weight, gam.16.gamma.Mussel.wet.weight.1)

#gam.16.lm.log.Mussel.wet.weight is best 

plot(gam.16.lm.log.Mussel.wet.weight, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.log.Mussel.wet.weight)
#look pretty good, esp qq
#qq_plot(gam.16.lm.log.Mussel.wet.weight, method = 'simulate')
#k.check(gam.16.lm.log.Mussel.wet.weight)
summary(gam.16.lm.log.Mussel.wet.weight)

#residuals many in a straight line but otherwise good 

gam.16.lm.log.Mussel.wet.weight.unordered<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
fam.gam.16.Mussel.wet.weight <- family(gam.16.lm.log.Mussel.wet.weight)
fam.gam.16.Mussel.wet.weight
ilink.gam.16.Mussel.wet.weight<- fam.gam.16.Mussel.wet.weight$linkinv
ilink.gam.16.Mussel.wet.weight


mod.Mussel.wet.weight<-gam.16.lm.log.Mussel.wet.weight
ndata.16.Mussel.wet.weight <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.Mussel.wet.weight <- add_column(ndata.16.Mussel.wet.weight, fit = predict(mod.Mussel.wet.weight, newdata = ndata.16.Mussel.wet.weight, type = 'response'))

predict(mod.Mussel.wet.weight, newdata = ndata.16.Mussel.wet.weight, type = 'response')
ndata.16.Mussel.wet.weight <- bind_cols(ndata.16.Mussel.wet.weight, setNames(as_tibble(predict(mod.Mussel.wet.weight, ndata.16.Mussel.wet.weight, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.Mussel.wet.weight <- mutate(ndata.16.Mussel.wet.weight,
                                      fit_resp  = ilink.gam.16.Mussel.wet.weight(fit_link),
                                      right_upr = ilink.gam.16.Mussel.wet.weight(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.Mussel.wet.weight(fit_link - (2 * se_link)))


ndata.16.Mussel.wet.weight$min.10.pH.unscaled<-ndata.16.Mussel.wet.weight$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight <- ggplot(ndata.16.Mussel.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(Mussel.wet.weight+0.1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Mytilus") ~"wet weight (g, log scale)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.Mussel.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight
ggsave("C:Graphs April 2020//Mussel.wet.weight_pred.png")






# Mussel wet weight per individual mussel ---------------------------------
#big outlier .... 5 mussels weigh 105.5 g?? so 17 g per mussel? 

invasion.exp.data.16_zscores$Mussel.wet.weight[invasion.exp.data.16_zscores$Mesocosm==8]<-14.0
invasion.exp.data.16_zscores$Mussel.wet.weight.per.1<-(invasion.exp.data.16_zscores$Mussel.wet.weight)/(invasion.exp.data.16_zscores$mussel+1)

gam.16.lm.Mussel.wet.weight.per.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.Mussel.wet.weight.per.1<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.lm.log.Mussel.wet.weight.per.1<- gam(log(Mussel.wet.weight.per.1+0.01) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.tweedie.Mussel.wet.weight.per.1.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = tw, select=TRUE, method="REML")
gam.16.loglink.Mussel.wet.weight.per.1.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.16.loglink.Mussel.wet.weight.per.1.1, gam.16.lm.log.Mussel.wet.weight.per.1, gam.16.tweedie.Mussel.wet.weight.per.1.1,  gam.16.lm.Mussel.wet.weight.per.1, gam.16.gamma.Mussel.wet.weight.per.1)

#gamma is best 
plot(gam.16.gamma.Mussel.wet.weight.per.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.gamma.Mussel.wet.weight.per.1)
# good but a bit funnelly
#qq_plot(gam.16.gamma.Mussel.wet.weight.per.1, method = 'simulate')
#k.check(gam.16.gamma.Mussel.wet.weight.per.1)
summary(gam.16.gamma.Mussel.wet.weight.per.1)

gam.16.gamma.Mussel.wet.weight.per.1.unordered<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.16.Mussel.wet.weight.per.1 <- family(gam.16.gamma.Mussel.wet.weight.per.1)
fam.gam.16.Mussel.wet.weight.per.1
str(fam.gam.16.Mussel.wet.weight.per.1)
ilink.gam.16.Mussel.wet.weight.per.1<- fam.gam.16.Mussel.wet.weight.per.1$linkinv
ilink.gam.16.Mussel.wet.weight.per.1


mod.Mussel.wet.weight.per.1<-gam.16.gamma.Mussel.wet.weight.per.1
ndata.16.Mussel.wet.weight.per.1 <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.Mussel.wet.weight.per.1 <- add_column(ndata.16.Mussel.wet.weight.per.1, fit = predict(mod.Mussel.wet.weight.per.1, newdata = ndata.16.Mussel.wet.weight.per.1, type = 'response'))

predict(mod.Mussel.wet.weight.per.1, newdata = ndata.16.Mussel.wet.weight.per.1, type = 'response')
ndata.16.Mussel.wet.weight.per.1 <- bind_cols(ndata.16.Mussel.wet.weight.per.1, setNames(as_tibble(predict(mod.Mussel.wet.weight.per.1, ndata.16.Mussel.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.Mussel.wet.weight.per.1 <- mutate(ndata.16.Mussel.wet.weight.per.1,
                                  fit_resp  = ilink.gam.16.Mussel.wet.weight.per.1(fit_link),
                                  right_upr = ilink.gam.16.Mussel.wet.weight.per.1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.16.Mussel.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.16.Mussel.wet.weight.per.1$min.10.pH.unscaled<-ndata.16.Mussel.wet.weight.per.1$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight.per.1 <- ggplot(ndata.16.Mussel.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = Mussel.wet.weight.per.1, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Biomass per mussel \n(wet weight)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.Mussel.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight.per.1
ggsave("C:Graphs April 2020//Mussel.wet.weight.per.1_pred.png")




# Hydtobot ----------------------------------------------------------------
invasion.exp.data.16_zscores$hydtobot[invasion.exp.data.16_zscores$hydtobot==1]<-0.99
invasion.exp.data.16_zscores$hydtobot[invasion.exp.data.16_zscores$hydtobot==0]<-0.01

gam.16.lm.hydtobot<- gam(log(hydtobot+0.01)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.log.lm.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.16.tweedie.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family =tw(), select=TRUE, method="REML")

gam.16.beta.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.hydtobot.1<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.hydtobot.2<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")


AICtab( gam.16.tweedie.hydtobot ,gam.16.beta.hydtobot, gam.16.lm.hydtobot, gam.16.log.lm.hydtobot, gam.16.beta.hydtobot.1, gam.16.beta.hydtobot.2)
### beta logit is the best 

plot(gam.16.beta.hydtobot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.hydtobot)
#qq a bit wiggly
#qq_plot(gam.16.beta.hydtobot, method = 'simulate')
#k.check(gam.16.beta.hydtobot)
gam.16.check(gam.16.tweedie.hydtobot)
summary(gam.16.beta.hydtobot)
#not the best qq plot
#resids are in rows

gam.16.beta.hydtobot.unordered<- gam(hydtobot~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
fam.gam.16.hydtobot <- family(gam.16.beta.hydtobot)
fam.gam.16.hydtobot 

ilink.gam.16.hydtobot <- fam.gam.16.hydtobot$linkinv
ilink.gam.16.hydtobot


invasion.exp.data.16_zscores$min.10.pH.unscaled <-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data.16_zscores)

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)

mod.hydtobot<-gam.16.beta.hydtobot
ndata.16.hydtobot <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                               length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.hydtobot <- add_column(ndata.16.hydtobot, fit = predict(mod.hydtobot, newdata = ndata.16.hydtobot, type = 'response'))


ndata.16.hydtobot <- bind_cols(ndata.16.hydtobot, setNames(as_tibble(predict(mod.hydtobot, ndata.16.hydtobot, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.hydtobot <- mutate(ndata.16.hydtobot,
                        fit_resp  = ilink.gam.16.hydtobot(fit_link),
                        right_upr = ilink.gam.16.hydtobot(fit_link + (2 * se_link)),
                        right_lwr = ilink.gam.16.hydtobot(fit_link - (2 * se_link)))




ndata.16.hydtobot$min.10.pH.unscaled<-ndata.16.hydtobot$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.16.hydtobot <- ggplot(ndata.16.hydtobot, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydtobot, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "cover ratio"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.hydtobot,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.2, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=6), legend.title = element_text(size=7))+ 
  geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.gam.16.hydtobot 


# Hydtobot by weight ------------------------------------------------------
invasion.exp.data.16_zscores$hydtobot_dry_biomass[invasion.exp.data.16_zscores$hydtobot_dry_biomass==1]<-0.99
invasion.exp.data.16_zscores$hydtobot_dry_biomass[invasion.exp.data.16_zscores$hydtobot_dry_biomass==0]<-0.01

gam.16.lm.hydtobot_dry_biomass<- gam(log(hydtobot_dry_biomass+0.01)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.log.lm.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.16.tweedie.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family =tw(), select=TRUE, method="REML")
gam.16.beta.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.hydtobot_dry_biomass.1<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.hydtobot_dry_biomass.2<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")

AICtab(gam.16.tweedie.hydtobot_dry_biomass ,gam.16.beta.hydtobot_dry_biomass, gam.16.lm.hydtobot_dry_biomass, gam.16.log.lm.hydtobot_dry_biomass, gam.16.beta.hydtobot_dry_biomass.1, gam.16.beta.hydtobot_dry_biomass.2, lm.hydtobot_dry_biomass, lm.log.hydtobot_dry_biomass)
### logit is the best 

plot(gam.16.beta.hydtobot_dry_biomass, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.hydtobot_dry_biomass)
#look ok, qq is bloacky
#qq_plot(gam.16.beta.hydtobot_dry_biomass, method = 'simulate')
#k.check(gam.16.beta.hydtobot_dry_biomass)
summary(gam.16.beta.hydtobot_dry_biomass)
gam.16.beta.hydtobot_dry_biomass.unordered<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.16.hydtobot_dry_biomass <- family(gam.16.beta.hydtobot_dry_biomass)
fam.gam.16.hydtobot_dry_biomass 
str(fam.gam.16.hydtobot_dry_biomass )
ilink.gam.16.hydtobot_dry_biomass <- fam.gam.16.hydtobot_dry_biomass$linkinv
ilink.gam.16.hydtobot_dry_biomass

invasion.exp.data.16_zscores$min.10.pH.unscaled <-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data.16_zscores)

want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)

mod.hydtobot_dry_biomass<-gam.16.beta.hydtobot_dry_biomass
ndata.16.hydtobot_dry_biomass <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.hydtobot_dry_biomass <- add_column(ndata.16.hydtobot_dry_biomass, fit = predict(mod.hydtobot_dry_biomass, newdata = ndata.16.hydtobot_dry_biomass, type = 'response'))


ndata.16.hydtobot_dry_biomass <- bind_cols(ndata.16.hydtobot_dry_biomass, setNames(as_tibble(predict(mod.hydtobot_dry_biomass, ndata.16.hydtobot_dry_biomass, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.hydtobot_dry_biomass <- mutate(ndata.16.hydtobot_dry_biomass,
                         fit_resp  = ilink.gam.16.hydtobot_dry_biomass(fit_link),
                         right_upr = ilink.gam.16.hydtobot_dry_biomass(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.hydtobot_dry_biomass(fit_link - (2 * se_link)))




ndata.16.hydtobot_dry_biomass$min.10.pH.unscaled<-ndata.16.hydtobot_dry_biomass$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.16.hydtobot_dry_biomass <- ggplot(ndata.16.hydtobot_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydtobot_dry_biomass, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "biomass ratio"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.hydtobot_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)#+coord_cartesian(ylim = c(0, 1)) 
plt.gam.16.hydtobot_dry_biomass 


#mesocosm #8 didn't have either tunicates or hydroids weight - 0% hydroids, 2% tunicates


# CAP1 --------------------------------------------------------------------
#negative values so can't do gamma
gam.16.lm.CAP1<- gam(CAP1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.CAP1.1<- gam(CAP1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( gam.16.loglink.CAP1.1, gam.16.lm.CAP1)
#gam.16.lm.CAP1

plot(gam.16.lm.CAP1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.CAP1)
#look good
#qq_plot(gam.16.lm.CAP1, method = 'simulate')
#k.check(gam.16.lm.CAP1)
summary(gam.16.lm.CAP1)
gam.16.lm.CAP1.unordered<- gam(CAP1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
summary(gam.16.lm.CAP1.unordered)

fam.gam.16.CAP1 <- family(gam.16.lm.CAP1)
fam.gam.16.CAP1
str(fam.gam.16.CAP1)
ilink.gam.16.CAP1<- fam.gam.16.CAP1$linkinv
ilink.gam.16.CAP1


mod.CAP1<-gam.16.lm.CAP1
ndata.16.CAP1 <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.CAP1 <- add_column(ndata.16.CAP1, fit = predict(mod.CAP1, newdata = ndata.16.CAP1, type = 'response'))

predict(mod.CAP1, newdata = ndata.16.CAP1, type = 'response')
ndata.16.CAP1 <- bind_cols(ndata.16.CAP1, setNames(as_tibble(predict(mod.CAP1, ndata.16.CAP1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.CAP1 <- mutate(ndata.16.CAP1,
                                      fit_resp  = ilink.gam.16.CAP1(fit_link),
                                      right_upr = ilink.gam.16.CAP1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.CAP1(fit_link - (2 * se_link)))


ndata.16.CAP1$min.10.pH.unscaled<-ndata.16.CAP1$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.CAP1 <- ggplot(ndata.16.CAP1, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(CAP1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(36% of constrained variation)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.1, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.CAP1
ggsave("C:Graphs April 2020//CAP1_pred.png")




# Distances ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.16.lm.distances<- gam(distcentroid ~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.distances.1<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.16.gamma.distances<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.16.loglink.distances.1,  gam.16.lm.distances, gam.16.gamma.distances)
#gam.16.lm.distances although both are equal


plot(gam.16.lm.distances, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.distances)
#looks good
#qq_plot(gam.16.lm.distances, method = 'simulate')
#k.check(gam.16.lm.distances)
#good now
summary(gam.16.lm.distances)

gam.16.lm.distances.unordered<- gam(distcentroid~ s(min.10.pH, k=11)+ Invasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
summary(gam.16.lm.distances.unordered)

fam.gam.16.distances <- family(gam.16.lm.distances)
fam.gam.16.distances
str(fam.gam.16.distances)
ilink.gam.16.distances<- fam.gam.16.distances$linkinv
ilink.gam.16.distances


mod.distances<-gam.16.lm.distances
ndata.16.distances <- with(invasion.exp.data.16_zscores, tibble(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                             length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.distances <- add_column(ndata.16.distances, fit = predict(mod.distances, newdata = ndata.16.distances, type = 'response'))

predict(mod.distances, newdata = ndata.16.distances, type = 'response')
ndata.16.distances <- bind_cols(ndata.16.distances, setNames(as_tibble(predict(mod.distances, ndata.16.distances, se.fit = TRUE)[1:2]),
                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.distances <- mutate(ndata.16.distances,
                     fit_resp  = ilink.gam.16.distances(fit_link),
                     right_upr = ilink.gam.16.distances(fit_link + (2 * se_link)),
                     right_lwr = ilink.gam.16.distances(fit_link - (2 * se_link)))


ndata.16.distances$min.10.pH.unscaled<-ndata.16.distances$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.distances <- ggplot(ndata.16.distances, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(distcentroid), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.distances
ggsave("C:Graphs April 2020//distances_pred.png")

# Community plotting ------------------------------------------------------
library(cowplot)

#### revised community fig
fig.4.community<-wrap_plots( plt.occupied.space,plt.total_dry_biomass,
                            plt.richness, plt.evenness,
                            plt.CAP1, plt.distances, ncol=2)+
                            plot_annotation(tag_levels = 'a')

fig.4.community

ggplot2::ggsave(plot=fig.4.community, "C:Data//For submission//For resubmission//RESUB2//First look//Fig4.community.pdf", width=3, height=5, units="in")



#hyd tobot figs
fig.s3.hydtobot<-wrap_plots(plt.gam.16.hydtobot, plt.gam.16.hydtobot_dry_biomass, ncol=2) + plot_annotation(tag_levels = 'a')

fig.s3.hydtobot

ggplot2::ggsave("C:Data//For submission//For resubmission//RESUB2//First look//Fig.S3.hydotobot.png", width=6, height=3, units="in", dpi=600)



# Community level tables --------------------------------------------------

richness.gam<- summary(gam.16.poisson.richness)
richness.gam.16.unordered<- summary(gam.16.poisson.richness.unordered)

evenness.gam<-summary(gam.16.lm.evenness)
evenness.gam.16.unordered<-summary(gam.16.lm.evenness.unordered)


occupied.space.gam<- summary(gam.16.beta.occupied.space.3)
occupied.space.gam.16.unordered<- summary(gam.16.beta.occupied.space.3.unordered)

distances.gam<- summary(gam.16.lm.distances)
distances.gam.16.unordered<- summary(gam.16.lm.distances.unordered)


CAP1.gam <- summary(gam.16.lm.CAP1)
CAP1.gam.16.unordered <- summary(gam.16.lm.CAP1.unordered)

#dry biomass
total_dry_biomass.gam <- summary(gam.16.lm.total_dry_biomass)
total_dry_biomass.gam.16.unordered <- summary(gam.16.lm.total_dry_biomass.unordered)

hydroid_dry_biomass.gam.16.unordered <- summary(gam.16.gamma.hydroid_dry_biomass.unordered) 
hydroid_dry_biomass.gam<- summary(gam.16.gamma.hydroid_dry_biomass) 

caprellid_dry_biomass.gam.16.unordered <- summary(gam.16.gamma.caprellid_dry_biomass.unordered)
caprellid_dry_biomass.gam <- summary(gam.16.gamma.caprellid_dry_biomass)

tunicate_dry_biomass.gam <- summary(gam.16.gamma.tunicate_dry_biomass)
tunicate_dry_biomass.gam.16.unordered <- summary(gam.16.gamma.tunicate_dry_biomass.unordered)

rest_dry_biomass.gam.16.unordered <- summary(gam.16.tweedie.rest_dry_biomass.unordered)
rest_dry_biomass.gam <- summary(gam.16.tweedie.rest_dry_biomass)

hydtobot_dry_biomass.gam.16.unordered<-summary(gam.16.beta.hydtobot_dry_biomass.unordered)
hydtobot_dry_biomass.gam<-summary(gam.16.beta.hydtobot_dry_biomass)

#wet biomass
everything.wet.weight.gam <-summary(gam.16.lm.log.everything.wet.weight)
everything.wet.weight.gam.16.unordered <- summary(gam.16.lm.log.everything.wet.weight.unordered)

everything.wet.weight.per.1.gam <-summary(gam.16.lm.log.everything.wet.weight.per.1)
everything.wet.weight.per.1.gam.16.unordered <- summary(gam.16.lm.log.everything.wet.weight.per.1.unordered)


Mussel.wet.weight.gam <- summary(gam.16.lm.log.Mussel.wet.weight)
Mussel.wet.weight.gam.16.unordered <-summary(gam.16.lm.log.Mussel.wet.weight.unordered)

Mussel.wet.weight.per.1.gam <- summary(gam.16.gamma.Mussel.wet.weight.per.1)
Mussel.wet.weight.per.1.gam.16.unordered <-summary(gam.16.gamma.Mussel.wet.weight.per.1.unordered)


#competition metric 
hydtobot.gam <- summary(gam.16.beta.hydtobot)
hydtobot.gam.16.unordered <- summary(gam.16.beta.hydtobot.unordered)

#ptable building
richness.gam.16.p.table<-as.data.frame(richness.gam.16.unordered$p.table)
richness.gam.16.s.table<-as.data.frame(richness.gam$s.table)

evenness.gam.16.p.table<-as.data.frame(evenness.gam.16.unordered$p.table)
evenness.gam.16.s.table<-as.data.frame(evenness.gam$s.table)

occupied.space.gam.16.p.table<-as.data.frame(occupied.space.gam.16.unordered$p.table)
occupied.space.gam.16.s.table<-as.data.frame(occupied.space.gam$s.table)

distances.gam.16.p.table<-as.data.frame(distances.gam.16.unordered$p.table)
distances.gam.16.s.table<-as.data.frame(distances.gam$s.table)

CAP1.gam.16.p.table<-as.data.frame(CAP1.gam.16.unordered$p.table)
CAP1.gam.16.s.table<-as.data.frame(CAP1.gam$s.table)

total_dry_biomass.gam.16.p.table<-as.data.frame(total_dry_biomass.gam.16.unordered$p.table)
total_dry_biomass.gam.16.s.table<-as.data.frame(total_dry_biomass.gam$s.table)

hydroid_dry_biomass.gam.16.p.table<-as.data.frame(hydroid_dry_biomass.gam.16.unordered$p.table)
hydroid_dry_biomass.gam.16.s.table<-as.data.frame(hydroid_dry_biomass.gam$s.table)

caprellid_dry_biomass.gam.16.p.table<-as.data.frame(caprellid_dry_biomass.gam.16.unordered$p.table)
caprellid_dry_biomass.gam.16.s.table<-as.data.frame(caprellid_dry_biomass.gam$s.table)

tunicate_dry_biomass.gam.16.p.table<-as.data.frame(tunicate_dry_biomass.gam.16.unordered$p.table)
tunicate_dry_biomass.gam.16.s.table<-as.data.frame(tunicate_dry_biomass.gam$s.table)

rest_dry_biomass.gam.16.p.table<-as.data.frame(rest_dry_biomass.gam.16.unordered$p.table)
rest_dry_biomass.gam.16.s.table<-as.data.frame(rest_dry_biomass.gam$s.table)

everything.wet.weight.gam.16.p.table<-as.data.frame(everything.wet.weight.gam.16.unordered$p.table)
everything.wet.weight.gam.16.s.table<-as.data.frame(everything.wet.weight.gam$s.table)

everything.wet.weight.per.1.gam.16.p.table<-as.data.frame(everything.wet.weight.per.1.gam.16.unordered$p.table)
everything.wet.weight.per.1.gam.16.s.table<-as.data.frame(everything.wet.weight.per.1.gam$s.table)

Mussel.wet.weight.gam.16.p.table<-as.data.frame(Mussel.wet.weight.gam.16.unordered$p.table)
Mussel.wet.weight.gam.16.s.table<-as.data.frame(Mussel.wet.weight.gam$s.table)

Mussel.wet.weight.per.1.gam.16.p.table<-as.data.frame(Mussel.wet.weight.per.1.gam.16.unordered$p.table)
Mussel.wet.weight.per.1.gam.16.s.table<-as.data.frame(Mussel.wet.weight.per.1.gam$s.table)

hydtobot.gam.16.p.table<-as.data.frame(hydtobot.gam.16.unordered$p.table)
hydtobot.gam.16.s.table<-as.data.frame(hydtobot.gam$s.table)

hydtobot_dry_biomass.gam.16.p.table<-as.data.frame(hydtobot_dry_biomass.gam.16.unordered$p.table)
hydtobot_dry_biomass.gam.16.s.table<-as.data.frame(hydtobot_dry_biomass.gam$s.table)

hydtobot_dry_biomass.gam.16.p.table
hydtobot_dry_biomass.gam.16.s.table

#richness.gam.16.p.table and  hydtobot.gam.16.p.table, is with z value 
colnames(richness.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot_dry_biomass.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(occupied.space.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.community.t<-rbind(richness.gam.16.p.table,
              evenness.gam.16.p.table,
              occupied.space.gam.16.p.table,
              total_dry_biomass.gam.16.p.table,
              CAP1.gam.16.p.table,
              distances.gam.16.p.table,
              hydtobot.gam.16.p.table, 
              hydtobot_dry_biomass.gam.16.p.table,
              everything.wet.weight.gam.16.p.table
              )


colnames(ptable.community.t) <- c("Estimate", "SE", "t", "p")
ptable.community.t$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.community.t %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (z)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:Data//For submission//ptable.community.t.html", self_contained = T)


#again hydtobot and richness
#richness.gam.16.p.table and  hydtobot.gam.16.p.table, is with Chisq
colnames(richness.gam.16.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(hydtobot.gam.16.s.table) <- c("edf", "Ref.df",  "F", "p-value")
colnames(hydtobot_dry_biomass.gam.16.s.table) <- c("edf", "Ref.df",  "F", "p-value")

colnames(occupied.space.gam.16.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.community.f<-rbind(richness.gam.16.s.table,
                          evenness.gam.16.s.table,
                          occupied.space.gam.16.s.table,
                          total_dry_biomass.gam.16.s.table,
                          CAP1.gam.16.s.table,
                          distances.gam.16.s.table,
                          hydtobot.gam.16.s.table, 
                          hydtobot_dry_biomass.gam.16.s.table,
                          everything.wet.weight.gam.16.s.table
)


colnames(stable.community.f) <- c("Estimated_df", "Reference_df", "F", "p_smooth")
stable.community.f$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Low quality food", "smooth min.10.pH * High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

stable.community.f %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi-square)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi-square)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:Data//For submission//stable.community.f.html", self_contained = T)

pstable.community<-cbind(ptable.community.t, stable.community.f)


pstable.community %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.051, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth, Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi-square, z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi-square, z)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square, z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square, z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  save_kable(file = "C:Data//For submission//For resubmission//RESUB2//First look//pstable.community.html ", self_contained = T)
