
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
invasion.exp.data<-invasion.exp.data %>% filter(Week==16)

head(invasion.exp.data)

str(invasion.exp.data)

#ordered and unordered factors
invasion.exp.data$oTreatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2.TreatmentAbsent", "CO2.TreatmentPresent", "AIRPresent"), ordered=TRUE)
invasion.exp.data$Treatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2.TreatmentAbsent", "CO2.TreatmentPresent", "AIRPresent"), ordered=FALSE)

#I think I need just Invasives if we're going with CO2.Treatment
invasion.exp.data$oInvasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=TRUE)
invasion.exp.data$Invasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=FALSE)




# #some variables need to be read over
# invasion.exp.data$caprellid.percent<-invasion.exp.data$caprellid
# invasion.exp.data$hydroid<-invasion.exp.data$hydroid
# invasion.exp.data$botryllid<-invasion.exp.data$botryllid
# invasion.exp.data$folliculina<-invasion.exp.data$folliculina
# invasion.exp.data$membranipora<-invasion.exp.data$membranipora
# invasion.exp.data$didemnum<-invasion.exp.data$didemnum
# invasion.exp.data$total<-invasion.exp.data$total
# invasion.exp.data$bare<-invasion.exp.data$bare
# invasion.exp.data$occupied.space<-100-invasion.exp.data$bare
# invasion.exp.data$occupied.space.001<-0.01*(invasion.exp.data$occupied.space)
# invasion.exp.data$everything.wet.weight<-invasion.exp.data$everything.wet.weight
# invasion.exp.data$everything.wet.weight.per.1<-(invasion.exp.data$everything.wet.weight)/invasion.exp.data$occupied.space
# invasion.exp.data$Mussel.wet.weight<-invasion.exp.data$Mussel.wet.weight
# invasion.exp.data$total_dry_biomass<-invasion.exp.data$total_dry_biomass
# invasion.exp.data$total_dry_biomass_per1<-invasion.exp.data$total_dry_biomass/invasion.exp.data$occupied.space
# invasion.exp.data$hydroid_dry_biomass<-invasion.exp.data$hydroid_dry_biomass
# invasion.exp.data$caprellid_dry_biomass<-invasion.exp.data$caprellid_dry_biomass
# invasion.exp.data$tunicate_dry_biomass<-invasion.exp.data$tunicate_dry_biomass
# invasion.exp.data$hydtobot<-(invasion.exp.data$botryllid)/(invasion.exp.data$botryllid+invasion.exp.data$hydroid)
# invasion.exp.data$rest_dry_biomass<-invasion.exp.data$rest_dry_biomass
# 
# #small negative biomass is within error of scale - change to zero
# invasion.exp.data$hydroid_dry_biomass[invasion.exp.data$hydroid_dry_biomass<0]<-0
# invasion.exp.data$tunicate_dry_biomass[invasion.exp.data$tunicate_dry_biomass<0]<-0
# invasion.exp.data$hydtobot_dry_biomass<-(invasion.exp.data$tunicate_dry_biomass)/(invasion.exp.data$tunicate_dry_biomass+invasion.exp.data$hydroid_dry_biomass)
# 
# 
# invasion.exp.data$Mussel.wet.weight.per.1<-(invasion.exp.data$Mussel.wet.weight)/(invasion.exp.data$mussel+1)
# 
# #making it a proportion instead of % cover
# invasion.exp.data$caprellid.percent.001<-(0.01*(invasion.exp.data$caprellid.percent))+0.01
 invasion.exp.data$hydroid.001<-(0.01*(invasion.exp.data$hydroid))+0.01
 invasion.exp.data$botryllid.001<-(0.01*(invasion.exp.data$botryllid))+0.01
 invasion.exp.data$membranipora.001<-(0.01*(invasion.exp.data$membranipora))+0.01
 invasion.exp.data$didemnum<-invasion.exp.data$white.bryo
 invasion.exp.data$num.red.bryoporella<-invasion.exp.data$red.bryo
 invasion.exp.data$folliculina<-invasion.exp.data$protozoa
 invasion.exp.data$folliculina.001<-(0.01*(invasion.exp.data$folliculina))+0.01
 
  invasion.exp.data$didemnum.001<-(0.01*(invasion.exp.data$didemnum))+0.01

# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
invasion.exp.data_zscores<-invasion.exp.data
#invasion.exp.data_zscores$hydrogen.concentration<-scale(invasion.exp.data$hydrogen.concentration, center=TRUE, scale=TRUE)
invasion.exp.data_zscores$av.pH<-scale(invasion.exp.data$av.pH, center=TRUE, scale=TRUE)
invasion.exp.data_zscores$min.10.pH<-scale(invasion.exp.data$min.10.pH, center=TRUE, scale=TRUE)
invasion.exp.data_zscores$Mesocosm <- as.factor(invasion.exp.data$Mesocosm)
invasion.exp.data_zscores$av.pH.unscaled <-invasion.exp.data_zscores$av.pH * attr(invasion.exp.data_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$av.pH, 'scaled:center')
invasion.exp.data_zscores$min.10.pH.unscaled <-invasion.exp.data_zscores$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# Visualizing histograms of pH to use as continuous vs. discrete variable
invasion.exp.data_pres<-invasion.exp.data %>% filter(Invasives=="Present")
invasion.exp.data_abs<-invasion.exp.data %>% filter(Invasives=="Absent")




hist(invasion.exp.data_abs$min.10.pH, breaks=5)#5
hist(invasion.exp.data_pres$min.10.pH, breaks=18)#18 for the whole 16 weeks
#these are actually the exact same values b/c same mesocosms


 ggplot(invasion.exp.data, aes(x=min.10.pH))+geom_density()+theme_classic()


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
theme_set(theme_classic(base_size = 6))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))



# GAM beta hydroid / gam.beta.hydroid --------------------------------------------------------

#distributions - binomial or beta with various families

gam.binomial.hydroid<- gam(formula = cbind(hydroid, 100-hydroid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")


gam.beta.hydroid<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydroid.1<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydroid.2<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.hydroid.3<- gam(hydroid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.beta.hydroid, gam.beta.hydroid.1, gam.beta.hydroid.2, gam.beta.hydroid.3, gam.binomial.hydroid)
#cauchit is the best

plot(gam.beta.hydroid.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.hydroid.3)
qq_plot(gam.beta.hydroid.3, method = 'simulate')
k_check(gam.beta.hydroid.3)

gam.beta.hydroid.3.unordered<- gam(hydroid.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.hydroid.3)
summary(gam.beta.hydroid.3.unordered)
# need an unordered and an unordered model to get estimates for overall Invasives effect


### plotting based on gam confidence intervals
#using code from https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

fam.gam.hydroid <- family(gam.beta.hydroid.3)
fam.gam.hydroid 
str(fam.gam.hydroid )
ilink.gam.hydroid <- fam.gam.hydroid$linkinv
ilink.gam.hydroid

invasion.exp.data_zscores$min.10.pH.unscaled <-invasion.exp.data_zscores$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data_zscores)

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)
            
mod.hydroid<-gam.beta.hydroid.3
ndata.hydroid <- with(invasion.exp.data_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydroid <- add_column(ndata.hydroid, fit = predict(mod.hydroid, newdata = ndata.hydroid, type = 'response'))
ndata.hydroid <- bind_cols(ndata.hydroid, setNames(as_tibble(predict(mod.hydroid, ndata.hydroid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform
ndata.hydroid <- mutate(ndata.hydroid,
                fit_resp  = ilink.gam.hydroid(fit_link),
                right_upr = ilink.gam.hydroid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.hydroid(fit_link - (2 * se_link)))

#make sure pH is unscaled in the plot
ndata.hydroid$min.10.pH.unscaled<-ndata.hydroid$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


plt.gam.hydroid <- ggplot(ndata.hydroid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydroid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Obelia")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.hydroid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.gam.hydroid 
ggsave("C:Graphs April 2020//hydroid_pred.png")



# GAM beta botryllus / gam.beta.botryllid.2 ----------------------------------------------------

#binomial first
gam.binomial.botryllid<- gam(formula = cbind(botryllid, 100-botryllid)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.beta.botryllid<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.beta.botryllid.1<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.beta.botryllid.2<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.beta.botryllid.3<- gam(botryllid.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.beta.botryllid, gam.beta.botryllid.1, gam.beta.botryllid.2, gam.beta.botryllid.3, gam.binomial.botryllid)
#cauchit is the best

plot(gam.beta.botryllid.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.botryllid.3)
qq_plot(gam.beta.botryllid.3, method = 'simulate')
k_check(gam.beta.botryllid.3)
summary(gam.beta.botryllid.3)

gam.beta.botryllid.3.unordered<- gam(botryllid.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)


fam.gam.botryllid <- family(gam.beta.botryllid.3)
fam.gam.botryllid
ilink.gam.botryllid<- fam.gam.botryllid$linkinv
ilink.gam.botryllid
want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)


mod.botryllid<-gam.beta.botryllid.3
ndata.botryllid <- with(invasion.exp.data_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                        length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.botryllidel for the new data
ndata.botryllid <- add_column(ndata.botryllid, fit = predict(mod.botryllid, newdata = ndata.botryllid, type = 'response'))


ndata.botryllid <- bind_cols(ndata.botryllid, setNames(as_tibble(predict(mod.botryllid, ndata.botryllid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.botryllid <- mutate(ndata.botryllid,
                fit_resp  = ilink.gam.botryllid(fit_link),
                right_upr = ilink.gam.botryllid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.botryllid(fit_link - (2 * se_link)))


invasion.exp.data_zscores$min.10.pH.unscaled<-invasion.exp.data_zscores$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

ndata.botryllid$min.10.pH.unscaled<-ndata.botryllid$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.botryllid <- ggplot(ndata.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.botryllid
ggsave("C:Graphs April 2020//botryllid_pred.png")



# GAM negbin caprellid / gam.nb.caprellid -----------------------------------------------------------

#variety of potential distributions
gam.nb.caprellid<- gam(total.caprellids ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, family = negbin(nbinom12.caprellid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.caprellid.1<- gam(total.caprellids ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.caprellid<- gam(total.caprellids ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, family = poisson, select=TRUE, method="REML")
gam.lm.caprellid<- gam(total.caprellids ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")
gam.log.lm.caprellid<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")


AICtab(gam.lm.caprellid, gam.log.lm.caprellid, gam.nb.caprellid.1,gam.nb.caprellid, gam.poisson.caprellid)
#gam.log.lm.caprellid by far the best

#appraise(gam.log.lm.caprellid)
#qq_plot(gam.log.lm.caprellid, method = 'simulate')
plot(gam.log.lm.caprellid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.log.lm.caprellid)
summary(gam.log.lm.caprellid.unordered)

gam.log.lm.caprellid.unordered<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = food.caprellid.data_zscores, select=TRUE, method="REML")

fam.gam.caprellid <- family(gam.log.lm.caprellid)
fam.gam.caprellid
ilink.gam.caprellid<- fam.gam.caprellid$linkinv
ilink.gam.caprellid

mod.caprellid<-gam.log.lm.caprellid
ndata.caprellid <- with(food.caprellid.data_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                        length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.caprellid <- add_column(ndata.caprellid, fit = predict(mod.caprellid, newdata = ndata.caprellid, type = 'response'))

predict(mod.caprellid, newdata = ndata.caprellid, type = 'response')
ndata.caprellid <- bind_cols(ndata.caprellid, setNames(as_tibble(predict(mod.caprellid, ndata.caprellid, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid <- mutate(ndata.caprellid,
                          fit_resp  = ilink.gam.caprellid(fit_link),
                          right_upr = ilink.gam.caprellid(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.caprellid(fit_link - (2 * se_link)))


food.caprellid.data_zscores$min.10.pH.unscaled<-food.caprellid.data_zscores$min.10.pH * attr(food.caprellid.data_zscores$min.10.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$min.10.pH, 'scaled:center')
ndata.caprellid$min.10.pH.unscaled<-ndata.caprellid$min.10.pH * attr(food.caprellid.data_zscores$min.10.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$min.10.pH, 'scaled:center')




# plot 
plt.caprellid <- ggplot(ndata.caprellid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(total.caprellids+1), shape=CO2.Treatment, colour=oInvasives), data = food.caprellid.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Caprella")~ "abundance"), textstyle("(Log # of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid

ggsave("C:Graphs April 2020//caprellid_pred.png")


# GAM caprellid percent -----------------------------------------------------------

gam.binomial.caprellid.percent<- gam(formula = cbind(caprellid.percent, 100-caprellid.percent)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.caprellid.percent<- gam(caprellid.percent.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.1<- gam(caprellid.percent.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.2<- gam(caprellid.percent.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.caprellid.percent.3<- gam(caprellid.percent.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.caprellid.percent, gam.beta.caprellid.percent.1, gam.beta.caprellid.percent.3, gam.beta.caprellid.percent.2,gam.binomial.caprellid.percent)
#12,12.1, 12.2 are tied ... go with simplest logit


plot(gam.beta.caprellid.percent, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

#appraise(gam.beta.caprellid.percent)
#qq_plot(gam.beta.caprellid.percent, method = 'simulate')
#does not look great
#k_check(gam.beta.caprellid.percent)
summary(gam.beta.caprellid.percent.unordered)

gam.beta.caprellid.percent.unordered<- gam(caprellid.percent.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.caprellid.percent <- family(gam.beta.caprellid.percent)
fam.gam.caprellid.percent
str(fam.gam.caprellid.percent)
ilink.gam.caprellid.percent<- fam.gam.caprellid.percent$linkinv
ilink.gam.caprellid.percent


mod.caprellid.percent<-gam.beta.caprellid.percent
ndata.caprellid.percent <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.caprellid.percentel for the new data
ndata.caprellid.percent <- add_column(ndata.caprellid.percent, fit = predict(mod.caprellid.percent, newdata = ndata.caprellid.percent, type = 'response'))


ndata.caprellid.percent <- bind_cols(ndata.caprellid.percent, setNames(as_tibble(predict(mod.caprellid.percent, ndata.caprellid.percent, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid.percent <- mutate(ndata.caprellid.percent,
                          fit_resp  = ilink.gam.caprellid.percent(fit_link),
                          right_upr = ilink.gam.caprellid.percent(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.caprellid.percent(fit_link - (2 * se_link)))


ndata.caprellid.percent$min.10.pH.unscaled<-ndata.caprellid.percent$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.caprellid.percent <- ggplot(ndata.caprellid.percent, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = caprellid.percent.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Caprella")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid.percent,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid.percent
ggsave("C:Graphs April 2020//caprellid.percent_pred.png")



# GAM beta folliculina / gam.beta.folliculina -----------------------------------------------------------

gam.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")
gam.beta.folliculina<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.folliculina.1<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.folliculina.2<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.folliculina.3<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.folliculina, gam.beta.folliculina.1, gam.beta.folliculina.2,gam.binomial.folliculina, gam.beta.folliculina.3)
#12, 12.1, 12.2 best - go with simplest logit


plot(gam.beta.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.folliculina)
qq_plot(gam.beta.folliculina, method = 'simulate')
k_check(gam.beta.folliculina)
summary(gam.beta.folliculina)

gam.beta.folliculina.unordered<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.folliculina <- family(gam.beta.folliculina)
fam.gam.folliculina
str(fam.gam.folliculina)
ilink.gam.folliculina<- fam.gam.folliculina$linkinv
ilink.gam.folliculina


mod.folliculina<-gam.beta.folliculina
ndata.folliculina <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.folliculinael for the new data
ndata.folliculina <- add_column(ndata.folliculina, fit = predict(mod.folliculina, newdata = ndata.folliculina, type = 'response'))


ndata.folliculina <- bind_cols(ndata.folliculina, setNames(as_tibble(predict(mod.folliculina, ndata.folliculina, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.folliculina <- mutate(ndata.folliculina,
                          fit_resp  = ilink.gam.folliculina(fit_link),
                          right_upr = ilink.gam.folliculina(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.folliculina(fit_link - (2 * se_link)))


ndata.folliculina$min.10.pH.unscaled<-ndata.folliculina$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.folliculina <- ggplot(ndata.folliculina, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = folliculina.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.folliculina
ggsave("C:Graphs April 2020//folliculina_pred.png")



# GAM beta membranipora / gam.beta.membranipora --------------------------------------------------------

#binomial first
gam.binomial.membranipora<- gam(formula = cbind(membranipora, 100-membranipora)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.membranipora<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.membranipora.1<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.membranipora.2<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.membranipora.3<- gam(membranipora.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.beta.membranipora, gam.beta.membranipora.1, gam.beta.membranipora.2, gam.binomial.membranipora, gam.beta.membranipora.3)
#12, 12.1, 12.2, 12.3 all equal, go with logit


plot(gam.beta.membranipora, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.membranipora)
qq_plot(gam.beta.membranipora, method = 'simulate')
k_check(gam.beta.membranipora)
summary(gam.beta.membranipora)
vis.gam(gam.beta.membranipora)

gam.beta.membranipora.unordered<- gam(membranipora.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.membranipora <- family(gam.beta.membranipora)
fam.gam.membranipora
ilink.gam.membranipora<- fam.gam.membranipora$linkinv
ilink.gam.membranipora


mod.membranipora<-gam.beta.membranipora
ndata.membranipora <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.membraniporael for the new data
ndata.membranipora <- add_column(ndata.membranipora, fit = predict(mod.membranipora, newdata = ndata.membranipora, type = 'response'))


ndata.membranipora <- bind_cols(ndata.membranipora, setNames(as_tibble(predict(mod.membranipora, ndata.membranipora, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.membranipora <- mutate(ndata.membranipora,
                          fit_resp  = ilink.gam.membranipora(fit_link),
                          right_upr = ilink.gam.membranipora(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.membranipora(fit_link - (2 * se_link)))

ndata.membranipora$min.10.pH.unscaled<-ndata.membranipora$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.membranipora <- ggplot(ndata.membranipora, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = membranipora.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.membranipora
ggsave("C:Graphs April 2020//membranipora_pred.png")


# GAM beta didemnum / gam.beta.didemnum ------------------------------------------------------------

#binomial first
gam.binomial.didemnum<- gam(formula = cbind(didemnum, 100-didemnum)~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.didemnum<- gam(didemnum.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.didemnum.1<- gam(didemnum.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.didemnum.2<- gam(didemnum.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.didemnum.3<- gam(didemnum.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.didemnum, gam.beta.didemnum.1, gam.beta.didemnum.2,  gam.binomial.didemnum, gam.beta.didemnum.3)
#logit, all the betas are equal go with logit

plot(gam.beta.didemnum, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.didemnum)
#qq_plot(gam.beta.didemnum, method = 'simulate')
#not very good
#k_check(gam.beta.didemnum)
summary(gam.beta.didemnum)
vis.gam(gam.beta.didemnum)

#appraise doesn't fit that well .... but I think there's just not enough data

gam.beta.didemnum.unordered<- gam(didemnum.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.didemnum <- family(gam.beta.didemnum)
fam.gam.didemnum
str(fam.gam.didemnum)
ilink.gam.didemnum<- fam.gam.didemnum$linkinv
ilink.gam.didemnum


mod.didemnum<-gam.beta.didemnum
ndata.didemnum <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.didemnumel for the new data
ndata.didemnum <- add_column(ndata.didemnum, fit = predict(mod.didemnum, newdata = ndata.didemnum, type = 'response'))


ndata.didemnum <- bind_cols(ndata.didemnum, setNames(as_tibble(predict(mod.didemnum, ndata.didemnum, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.didemnum <- mutate(ndata.didemnum,
                          fit_resp  = ilink.gam.didemnum(fit_link),
                          right_upr = ilink.gam.didemnum(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.didemnum(fit_link - (2 * se_link)))

ndata.didemnum$min.10.pH.unscaled<-ndata.didemnum$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.didemnum <- ggplot(ndata.didemnum, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = didemnum.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Didemnum")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.didemnum,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.didemnum
ggsave("C:Graphs April 2020//didemnum_pred.png")



# GAM negbin mussel incomplete / gam.nb.mussel ---------------------------------------------------

poisson<-fitdistr(invasion.exp.data_zscores$mussel, "Poisson")
qqp(invasion.exp.data_zscores$mussel, "pois", lambda=poisson$estimate[[1]])
#estimating lambda

nbinom12.mussel <- fitdistr(invasion.exp.data_zscores$mussel, "Negative Binomial")
qqp(invasion.exp.data_zscores$mussel, "nbinom", size = nbinom12.mussel$estimate[[1]], mu = nbinom12.mussel$estimate[[2]])
#estimating theta

gam.nb.mussel<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")

gam.nb.mussel.1<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(link="log"), select=TRUE, method="REML")
gam.nb.mussel.2<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(link="sqrt"), select=TRUE, method="REML")
gam.nb.mussel.3<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")

gam.poisson.mussel<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson(link="log"), select=TRUE, method="REML")
gam.poisson.mussel.1<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson(link="identity"), select=TRUE, method="REML")
#gam.poisson.mussel.2<- gam(mussel ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson(link="sqrt"), select=TRUE, method="REML")
#won't run

AICtab(gam.nb.mussel, glm.nb.mussel.hydrogen, gam.nb.mussel.1, gam.nb.mussel.2, gam.nb.mussel.3, gam.poisson.mussel.1, gam.poisson.mussel)
#gam with theta estaimted from data is best
#gam.nb.mussel 


plot(gam.nb.mussel, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.mussel)
#qq_plot(gam.nb.mussel, method = 'simulate')
#looks good!
#k_check(gam.nb.mussel)
summary(gam.nb.mussel)


gam.nb.mussel.unordered<- gam(mussel ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")
summary(gam.nb.mussel.unordered)


fam.gam.mussel <- family(gam.nb.mussel)
fam.gam.mussel
ilink.gam.mussel<- fam.gam.mussel$linkinv
ilink.gam.mussel

mod.mussel<-gam.nb.mussel
ndata.mussel <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.mussel <- add_column(ndata.mussel, fit = predict(mod.mussel, newdata = ndata.mussel, type = 'response'))

predict(mod.mussel, newdata = ndata.mussel, type = 'response')
ndata.mussel <- bind_cols(ndata.mussel, setNames(as_tibble(predict(mod.mussel, ndata.mussel, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.mussel <- mutate(ndata.mussel,
                          fit_resp  = ilink.gam.mussel(fit_link),
                          right_upr = ilink.gam.mussel(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.mussel(fit_link - (2 * se_link)))

ndata.mussel$min.10.pH.unscaled<-ndata.mussel$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.mussel <- ggplot(ndata.mussel, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = mussel, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Mytilus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.mussel,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.mussel
ggsave("C:Graphs April 2020//mussel_pred.png")


# GAM negbin barnacles / gam.nb.num.barn -----------------------------------------------------------


nbinom12.barn <- fitdistr(invasion.exp.data_zscores$num.barn, "Negative Binomial")
qqp(invasion.exp.data_zscores$num.barn, "nbinom", size = nbinom12.barn.alive$estimate[[1]], mu = nbinom12.barn$estimate[[2]])

#negative binomial 
gam.nb.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.barn$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.barn.1<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.nb.num.barn, gam.nb.num.barn.1, gam.poisson.num.barn)

plot(gam.poisson.num.barn, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.poisson.num.barn)
qq_plot(gam.poisson.num.barn, method = 'simulate')
#looks really good
#k_check(gam.poisson.num.barn)
summary(gam.poisson.num.barn)


#a few outside the area
#appraise a bit funnelly

gam.poisson.num.barn.unordered<- gam(num.barn ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")


fam.gam.num.barn <- family(gam.poisson.num.barn)
fam.gam.num.barn
str(fam.gam.num.barn)
ilink.gam.num.barn<- fam.gam.num.barn$linkinv
ilink.gam.num.barn


mod.num.barn<-gam.poisson.num.barn
ndata.num.barn <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

ndata.num.barn


## add the fitted values by predicting from the model for the new data
ndata.num.barn <- add_column(ndata.num.barn, fit = predict(mod.num.barn, newdata = ndata.num.barn, type = 'response'))

predict(mod.num.barn, newdata = ndata.num.barn, type = 'response')
ndata.num.barn <- bind_cols(ndata.num.barn, setNames(as_tibble(predict(mod.num.barn, ndata.num.barn, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.barn <- mutate(ndata.num.barn,
                                  fit_resp  = ilink.gam.num.barn(fit_link),
                                  right_upr = ilink.gam.num.barn(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.num.barn(fit_link - (2 * se_link)))


ndata.num.barn$min.10.pH.unscaled<-ndata.num.barn$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.num.barn <- ggplot(ndata.num.barn, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.barn, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.barn
ggsave("C:Graphs April 2020//num.barn_pred.png")



# GAM negbin disporella / gam.nb.disporella ----------------------------------------------------------

nbinom12.disporella <- fitdistr(invasion.exp.data_zscores$disporella, "Negative Binomial")
qqp(invasion.exp.data_zscores$disporella, "nbinom", size = nbinom12.disporella$estimate[[1]], mu = nbinom12.disporella$estimate[[2]])
#getting theta

gam.nb.disporella<- gam(disporella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.disporella.1<- gam(disporella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.disporella<- gam(disporella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.nb.disporella.1,gam.nb.disporella,gam.poisson.disporella)
#used estimated theta

plot(gam.nb.disporella, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.disporella)
#not bad!
#qq_plot(gam.nb.disporella, method = 'simulate')
#k_check(gam.nb.disporella)
summary(gam.nb.disporella)

#a few outside the area
#appraise a bit funnelly

gam.nb.disporella.unordered<- gam(disporella ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")


fam.gam.disporella <- family(gam.nb.disporella)
fam.gam.disporella
str(fam.gam.disporella)
ilink.gam.disporella<- fam.gam.disporella$linkinv
ilink.gam.disporella


mod.disporella<-gam.nb.disporella
ndata.disporella <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

ndata.disporella

str(invasion.exp.data_zscores)
str(ndata.disporella)

## add the fitted values by predicting from the model for the new data
ndata.disporella <- add_column(ndata.disporella, fit = predict(mod.disporella, newdata = ndata.disporella, type = 'response'))

predict(mod.disporella, newdata = ndata.disporella, type = 'response')
ndata.disporella <- bind_cols(ndata.disporella, setNames(as_tibble(predict(mod.disporella, ndata.disporella, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.disporella <- mutate(ndata.disporella,
                               fit_resp  = ilink.gam.disporella(fit_link),
                               right_upr = ilink.gam.disporella(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.disporella(fit_link - (2 * se_link)))

ndata.disporella$min.10.pH.unscaled<-ndata.disporella$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.disporella <- ggplot(ndata.disporella, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = disporella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.disporella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.disporella
ggsave("C:Graphs April 2020//disporella_pred.png")


# GAM negbin num.red.bryo / gam.nb.num.red.bryo --------------------------------------------------------------

nbinom12.num.red.bryo <- fitdistr(invasion.exp.data_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data_zscores$num.red.bryo, "nbinom", size = nbinom12.num.red.bryo$estimate[[1]], mu = nbinom12.num.red.bryo$estimate[[2]])
#getting theta

gam.nb.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.red.bryo, gam.nb.num.red.bryo.1, gam.poisson.num.red.bryo)


plot(gam.poisson.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.poisson.num.red.bryo)
#looks pretty good
qq_plot(gam.poisson.num.red.bryo, method = 'simulate')
#k_check(gam.poisson.num.red.bryo)
summary(gam.poisson.num.red.bryo)

gam.poisson.num.red.bryo.unordered<- gam(num.red.bryo ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
summary(gam.poisson.num.red.bryo.unordered)

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

fam.gam.num.red.bryo <- family(gam.poisson.num.red.bryo)
fam.gam.num.red.bryo
str(fam.gam.num.red.bryo)
ilink.gam.num.red.bryo<- fam.gam.num.red.bryo$linkinv
ilink.gam.num.red.bryo

mod.num.red.bryo<-gam.poisson.num.red.bryo
ndata.num.red.bryo <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                   length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

ndata.num.red.bryo

str(invasion.exp.data_zscores)
str(ndata.num.red.bryo)

## add the fitted values by predicting from the model for the new data
ndata.num.red.bryo <- add_column(ndata.num.red.bryo, fit = predict(mod.num.red.bryo, newdata = ndata.num.red.bryo, type = 'response'))

predict(mod.num.red.bryo, newdata = ndata.num.red.bryo, type = 'response')
ndata.num.red.bryo <- bind_cols(ndata.num.red.bryo, setNames(as_tibble(predict(mod.num.red.bryo, ndata.num.red.bryo, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.red.bryo <- mutate(ndata.num.red.bryo,
                           fit_resp  = ilink.gam.num.red.bryo(fit_link),
                           right_upr = ilink.gam.num.red.bryo(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.num.red.bryo(fit_link - (2 * se_link)))


ndata.num.red.bryo$min.10.pH.unscaled<-ndata.num.red.bryo$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.red.bryo <- ggplot(ndata.num.red.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.red.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.red.bryo
ggsave("C:Graphs April 2020//num.red.bryo_pred.png")


# GAM poisson num nudi / gam.poisson.num.nudi  ------------------------------------------------------------
nbinom12.num.nudi <- fitdistr(invasion.exp.data_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data_zscores$num.nudi, "nbinom", size = nbinom12.num.nudi$estimate[[1]], mu = nbinom12.num.nudi$estimate[[2]])
#theta

gam.nb.num.nudi.1<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.nb.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.poisson.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.nudi, gam.nb.num.nudi.1, gam.poisson.num.nudi)
###poisson is the best fit by 2 dAIC


#appraise(gam.poisson.num.nudi)
#qq_plot(gam.poisson.num.nudi, method = 'simulate')
#looks quite good!
plot(gam.poisson.num.nudi, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.poisson.num.nudi)
summary(gam.poisson.num.nudi)
#resids a bit funny but same in neg bin

gam.poisson.num.nudi.unordered<- gam(num.nudi ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")
want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

fam.gam.num.nudi <- family(gam.poisson.num.nudi)
ilink.gam.num.nudi<- fam.gam.num.nudi$linkinv

mod.num.nudi<-gam.poisson.num.nudi
ndata.num.nudi <- with(invasion.exp.data_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.num.nudi <- add_column(ndata.num.nudi, fit = predict(mod.num.nudi, newdata = ndata.num.nudi, type = 'response'))

predict(mod.num.nudi, newdata = ndata.num.nudi, type = 'response')
ndata.num.nudi <- bind_cols(ndata.num.nudi, setNames(as_tibble(predict(mod.num.nudi, ndata.num.nudi, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.nudi <- mutate(ndata.num.nudi,
                       fit_resp  = ilink.gam.num.nudi(fit_link),
                       right_upr = ilink.gam.num.nudi(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.num.nudi(fit_link - (2 * se_link)))


ndata.num.nudi$min.10.pH.unscaled<-ndata.num.nudi$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.nudi <- ggplot(ndata.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.nudi, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.num.nudi
ggsave("C:Graphs April 2020//num.nudi_pred.png")


# GAM nb() serpulids / gam.nb.num.serpulid.1 -----------------------------------------------------------

nbinom12.num.serpulid <- fitdistr(invasion.exp.data_zscores$num.serpulid, "Negative Binomial")
qqp(invasion.exp.data_zscores$num.serpulid, "nbinom", size = nbinom12.num.serpulid$estimate[[1]], mu = nbinom12.num.serpulid$estimate[[2]])
#theta

#negative binomial first
gam.nb.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.serpulid.1<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.serpulid<- gam(num.serpulid ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.serpulid, gam.nb.num.serpulid.1, gam.poisson.num.serpulid)

##gam.nb.num.serpulid.1 is best

appraise(gam.nb.num.serpulid.1)
qq_plot(gam.nb.num.serpulid.1, method = 'simulate')
#looks good
plot(gam.nb.num.serpulid.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k_check(gam.nb.num.serpulid.1)
summary(gam.nb.num.serpulid.1)

#residuals a bit weird .... but they are the same in glm as in gam (doesn't improve)
#I think because of the zeros? 

gam.nb.num.serpulid.1.unordered<- gam(num.serpulid ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")


want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

fam.gam.num.serpulid <- family(gam.nb.num.serpulid.1)
fam.gam.num.serpulid
str(fam.gam.num.serpulid)
ilink.gam.num.serpulid<- fam.gam.num.serpulid$linkinv
ilink.gam.num.serpulid

mod.num.serpulid<-gam.nb.num.serpulid.1
ndata.num.serpulid <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.num.serpulid <- add_column(ndata.num.serpulid, fit = predict(mod.num.serpulid, newdata = ndata.num.serpulid, type = 'response'))

predict(mod.num.serpulid, newdata = ndata.num.serpulid, type = 'response')
ndata.num.serpulid <- bind_cols(ndata.num.serpulid, setNames(as_tibble(predict(mod.num.serpulid, ndata.num.serpulid, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.serpulid <- mutate(ndata.num.serpulid,
                         fit_resp  = ilink.gam.num.serpulid(fit_link),
                         right_upr = ilink.gam.num.serpulid(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.num.serpulid(fit_link - (2 * se_link)))


ndata.num.serpulid$min.10.pH.unscaled<-ndata.num.serpulid$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.serpulid
ggsave("C:Graphs April 2020//num.serpulid_pred.png")

# colorset_invasives = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")
# colorset_none = c("High"="#FFFFFF" ,"Low"="#FFFFFF","None"= "#666666")
# colorset_low = c("High"="#FFFFFF" ,"Low"="#439E5F","None"= "#666666")
# 
# ### Plot for powerpoint:
# plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
#   
#   geom_line(aes(colour=oInvasives)) +
#   geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
#   xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
#   scale_color_manual(values=colorset_none)+
#   scale_fill_manual(values=colorset_none)+
#   scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
#   geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
#   theme(legend.position='none')
# plt.num.serpulid
# ggsave("C:Graphs April 2020//num.serpulid_pred_none.png")
# 
# plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
#   
#   geom_line(aes(colour=oInvasives)) +
#   geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
#   xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
#   scale_color_manual(values=colorset_low)+
#   scale_fill_manual(values=colorset_low)+
#   scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
#   geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
#   theme(legend.position='none')
# plt.num.serpulid
# ggsave("C:Graphs April 2020//num.serpulid_pred_low.png")
# 
# GAM negbin orange sponge / gam.nb.orange_sponge -------------------------------------------------------

nbinom12.orange_sponge <- fitdistr(invasion.exp.data_zscores$orange_sponge, "Negative Binomial")
qqp(invasion.exp.data_zscores$orange_sponge, "nbinom", size = nbinom12.orange_sponge$estimate[[1]], mu = nbinom12.orange_sponge$estimate[[2]])
#theta

gam.nb.orange_sponge<- gam(orange_sponge ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")
gam.nb.orange_sponge.1<- gam(orange_sponge ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.orange_sponge<- gam(orange_sponge ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.orange_sponge, gam.nb.orange_sponge.1,  gam.poisson.orange_sponge)

###

#appraise(gam.nb.orange_sponge)
#a bit funnel-y
#qq_plot(gam.nb.orange_sponge, method = 'simulate')
plot(gam.nb.orange_sponge, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.nb.orange_sponge)
summary(gam.nb.orange_sponge)

gam.nb.orange_sponge.unordered<- gam(orange_sponge ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)
fam.gam.orange_sponge <- family(gam.nb.orange_sponge)
fam.gam.orange_sponge
str(fam.gam.orange_sponge)
ilink.gam.orange_sponge<- fam.gam.orange_sponge$linkinv
ilink.gam.orange_sponge

mod.orange_sponge<-gam.nb.orange_sponge
ndata.orange_sponge <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.orange_sponge <- add_column(ndata.orange_sponge, fit = predict(mod.orange_sponge, newdata = ndata.orange_sponge, type = 'response'))

predict(mod.orange_sponge, newdata = ndata.orange_sponge, type = 'response')
ndata.orange_sponge <- bind_cols(ndata.orange_sponge, setNames(as_tibble(predict(mod.orange_sponge, ndata.orange_sponge, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.orange_sponge <- mutate(ndata.orange_sponge,
                             fit_resp  = ilink.gam.orange_sponge(fit_link),
                             right_upr = ilink.gam.orange_sponge(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.orange_sponge(fit_link - (2 * se_link)))


ndata.orange_sponge$min.10.pH.unscaled<-ndata.orange_sponge$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 
plt.orange_sponge <- ggplot(ndata.orange_sponge, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = orange_sponge, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Sponge abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.orange_sponge,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.orange_sponge
ggsave("C:Graphs April 2020//orange_sponge_pred.png")



# GAM negbin corella / gam.nb.num.corella -------------------------------------------------------------

nbinom12.num.corella <- fitdistr(invasion.exp.data_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data_zscores$num.corella, "nbinom", size = nbinom12.num.corella$estimate[[1]], mu = nbinom12.num.corella$estimate[[2]])
#extracting theta

gam.nb.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.corella.1<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.corella, gam.nb.num.corella.1, gam.poisson.num.corella)
#

#appraise(gam.nb.num.corella)
#looks pretty good - slight pattern
#qq_plot(gam.nb.num.corella, method = 'simulate')
plot(gam.nb.num.corella, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.nb.num.corella)
summary(gam.nb.num.corella)


gam.nb.num.corella.unordered<- gam(num.corella ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)
fam.gam.num.corella <- family(gam.nb.num.corella)
ilink.gam.num.corella<- fam.gam.num.corella$linkinv


mod.num.corella<-gam.nb.num.corella
ndata.num.corella <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.num.corella <- add_column(ndata.num.corella, fit = predict(mod.num.corella, newdata = ndata.num.corella, type = 'response'))

predict(mod.num.corella, newdata = ndata.num.corella, type = 'response')
ndata.num.corella <- bind_cols(ndata.num.corella, setNames(as_tibble(predict(mod.num.corella, ndata.num.corella, se.fit = TRUE)[1:2]),
                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.corella <- mutate(ndata.num.corella,
                              fit_resp  = ilink.gam.num.corella(fit_link),
                              right_upr = ilink.gam.num.corella(fit_link + (2 * se_link)),
                              right_lwr = ilink.gam.num.corella(fit_link - (2 * se_link)))


ndata.num.corella$min.10.pH.unscaled<-ndata.num.corella$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

par(lheight=0.2) 
# plot 
plt.num.corella <- ggplot(ndata.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.corella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.num.corella
ggsave("C:Graphs April 2020//num.corella_pred.png")



# GAM poisson clam / gam.poisson.clam ----------------------------------------------------------------

nbinom12.clam <- fitdistr(invasion.exp.data_zscores$clam, "Negative Binomial")
qqp(invasion.exp.data_zscores$clam, "nbinom", size = nbinom12.clam$estimate[[1]], mu = nbinom12.clam$estimate[[2]])
#theta

#negative binomial first
gam.nb.clam<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = negbin(nbinom12.clam$estimate[[1]]), select=TRUE, method="REML")
gam.nb.clam.1<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.clam<- gam(clam ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.clam, gam.nb.clam.1, gam.poisson.clam)

###gam poisson is best

#appraise(gam.poisson.clam)
#qq_plot(gam.poisson.clam, method = 'simulate')
plot(gam.poisson.clam, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.poisson.clam)
summary(gam.poisson.clam)

#residuals a bit patterny as well 
gam.poisson.clam.unordered<- gam(clam ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")
summary(gam.poisson.clam.unordered)

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

fam.gam.clam <- family(gam.poisson.clam)
fam.gam.clam
str(fam.gam.clam)
ilink.gam.clam<- fam.gam.clam$linkinv
ilink.gam.clam

mod.clam<-gam.poisson.clam
ndata.clam <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                    length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the model for the new data
ndata.clam <- add_column(ndata.clam, fit = predict(mod.clam, newdata = ndata.clam, type = 'response'))

predict(mod.clam, newdata = ndata.clam, type = 'response')
ndata.clam <- bind_cols(ndata.clam, setNames(as_tibble(predict(mod.clam, ndata.clam, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.clam <- mutate(ndata.clam,
                            fit_resp  = ilink.gam.clam(fit_link),
                            right_upr = ilink.gam.clam(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.clam(fit_link - (2 * se_link)))


ndata.clam$min.10.pH.unscaled<-ndata.clam$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 
plt.clam <- ggplot(ndata.clam, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = clam, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.clam
ggsave("C:Graphs April 2020//clam_pred.png")


# Fig 2 plot generation ---------------------------------------------------
library(patchwork)
fig.2<-wrap_plots(plt.gam.hydroid,plt.botryllid,plt.folliculina,plt.caprellid.percent,plt.membranipora,plt.didemnum,
          plt.mussel,plt.num.barn,plt.disporella,plt.num.red.bryo,plt.num.nudi,plt.num.serpulid,
          plt.orange_sponge,plt.num.corella,plt.clam, ncol=5)+
          plot_annotation(tag_levels = 'a')

fig.2

ggplot2::ggsave(plot=fig.2, "C:Data//For submission//For resubmission//RESUB2//First look//Fig2_resized.pdf", width=8.75, height=4.75, units="in")

#ggplot2::ggsave("C:Data//For submission//For resubmission//Fig2.png", width=65, height=35, units="cm")


# Pulling model results to a table ----------------------------------------

hydroid.gam<-summary(gam.beta.hydroid.3)
botryllus.gam<-summary(gam.beta.botryllid.3)
caprellid.gam<-summary(gam.beta.caprellid.percent)
folliculina.gam<-summary(gam.beta.folliculina)
membranipora.gam<-summary(gam.beta.membranipora)
didemnum.gam<-summary(gam.beta.didemnum)
mussel.gam<-summary(gam.nb.mussel)
alive.barn.gam<-summary(gam.nb.num.barn)
disporella.gam<-summary(gam.nb.disporella)
num.red.bryo.gam<-summary(gam.nb.num.red.bryo)
num.nudi.gam<-summary(gam.poisson.num.nudi)
num.serpulid.gam<-summary(gam.nb.num.serpulid.1)
clam.gam<-summary(gam.poisson.clam)
corella.gam<-summary(gam.nb.num.corella)
orange_sponge.gam<-summary(gam.nb.orange_sponge)

hydroid.gam.unordered<-summary(gam.beta.hydroid.3.unordered)
botryllus.gam.unordered<-summary(gam.beta.botryllid.3.unordered)
caprellid.gam.unordered<-summary(gam.beta.caprellid.percent.unordered)
folliculina.gam.unordered<-summary(gam.beta.folliculina.unordered)
membranipora.gam.unordered<-summary(gam.beta.membranipora.unordered)
didemnum.gam.unordered<-summary(gam.beta.didemnum.unordered)
mussel.gam.unordered<-summary(gam.nb.mussel.unordered)
alive.barn.gam.unordered<-summary(gam.nb.num.barn.unordered)
disporella.gam.unordered<-summary(gam.nb.disporella.unordered)
num.red.bryo.gam.unordered<-summary(gam.nb.num.red.bryo.unordered)
num.nudi.gam.unordered<-summary(gam.poisson.num.nudi.unordered)
num.serpulid.gam.unordered<-summary(gam.nb.num.serpulid.1.unordered)
clam.gam.unordered<-summary(gam.poisson.clam.unordered)
corella.gam.unordered<-summary(gam.nb.num.corella.unordered)
orange_sponge.gam.unordered<-summary(gam.nb.orange_sponge.unordered)

hydroid.gam.p.table<-as.data.frame(hydroid.gam.unordered$p.table)
hydroid.gam.s.table<-as.data.frame(hydroid.gam$s.table)

botryllus.gam.p.table<-as.data.frame(botryllus.gam.unordered$p.table)
botryllus.gam.s.table<-as.data.frame(botryllus.gam$s.table)

caprellid.gam.p.table<-as.data.frame(caprellid.gam.unordered$p.table)
caprellid.gam.s.table<-as.data.frame(caprellid.gam$s.table)

folliculina.gam.p.table<-as.data.frame(folliculina.gam.unordered$p.table)
folliculina.gam.s.table<-as.data.frame(folliculina.gam$s.table)

membranipora.gam.p.table<-as.data.frame(membranipora.gam.unordered$p.table)
membranipora.gam.s.table<-as.data.frame(membranipora.gam$s.table)

didemnum.gam.p.table<-as.data.frame(didemnum.gam.unordered$p.table)
didemnum.gam.s.table<-as.data.frame(didemnum.gam$s.table)

mussel.gam.p.table<-as.data.frame(mussel.gam.unordered$p.table)
mussel.gam.s.table<-as.data.frame(mussel.gam$s.table)

alive.barn.gam.p.table<-as.data.frame(alive.barn.gam.unordered$p.table)
alive.barn.gam.s.table<-as.data.frame(alive.barn.gam$s.table)

disporella.gam.p.table<-as.data.frame(disporella.gam.unordered$p.table)
disporella.gam.s.table<-as.data.frame(disporella.gam$s.table)

num.red.bryo.gam.p.table<-as.data.frame(num.red.bryo.gam.unordered$p.table)
num.red.bryo.gam.s.table<-as.data.frame(num.red.bryo.gam$s.table)

num.nudi.gam.p.table<-as.data.frame(num.nudi.gam.unordered$p.table)
num.nudi.gam.s.table<-as.data.frame(num.nudi.gam$s.table)

num.serpulid.gam.p.table<-as.data.frame(num.serpulid.gam.unordered$p.table)
num.serpulid.gam.s.table<-as.data.frame(num.serpulid.gam$s.table)

orange_sponge.gam.p.table<-as.data.frame(orange_sponge.gam.unordered$p.table)
orange_sponge.gam.s.table<-as.data.frame(orange_sponge.gam$s.table)

corella.gam.p.table<-as.data.frame(corella.gam.unordered$p.table)
corella.gam.s.table<-as.data.frame(corella.gam$s.table)

clam.gam.p.table<-as.data.frame(clam.gam.unordered$p.table)
clam.gam.s.table<-as.data.frame(clam.gam$s.table)



head(clam.gam.p.table)
#### Building the stats table
ptable<-rbind(hydroid.gam.p.table, 
               botryllus.gam.p.table, 
               caprellid.gam.p.table,
               folliculina.gam.p.table,
               membranipora.gam.p.table,
               didemnum.gam.p.table,
               mussel.gam.p.table,
               alive.barn.gam.p.table,
               disporella.gam.p.table,
               num.red.bryo.gam.p.table,
               num.nudi.gam.p.table,
               num.serpulid.gam.p.table,
               orange_sponge.gam.p.table,
               corella.gam.p.table,
               clam.gam.p.table)


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
  group_rows("Disporella", 25, 27) %>% 
  group_rows("num.red.bryoporella", 28, 30) %>% 
  group_rows("Hermissenda", 31, 33) %>% 
  group_rows("Serpulid", 34, 36) %>% 
  group_rows("Sponge", 37, 39) %>% 
  group_rows("Corella", 40, 42) %>% 
  group_rows("Clams", 43, 45) %>% 
save_kable(file = "C:Data//For submission//ptable.html", self_contained = T)


### s table
stable<-rbind(hydroid.gam.s.table, 
              botryllus.gam.s.table, 
              caprellid.gam.s.table,
              folliculina.gam.s.table,
              membranipora.gam.s.table,
              didemnum.gam.s.table,
              mussel.gam.s.table,
              alive.barn.gam.s.table,
              disporella.gam.s.table,
              num.red.bryo.gam.s.table,
              num.nudi.gam.s.table,
              num.serpulid.gam.s.table,
              orange_sponge.gam.s.table,
              corella.gam.s.table,
              clam.gam.s.table)


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
  #                   "disporella",
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
  group_rows("Disporella", 25, 27) %>% 
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
  group_rows("Disporella, negative binomial", 25, 27) %>% 
  group_rows("num.red.bryoporella, negative binomial", 28, 30) %>% 
  group_rows("Hermissenda, negative binomial", 31, 33) %>% 
  group_rows("Serpulid, negative binomial", 34, 36) %>% 
  group_rows("Sponge, negative binomial", 37, 39) %>% 
  group_rows("Corella, negative binomial", 40, 42) %>% 
  group_rows("Clams, poisson", 43, 45) %>% 
  save_kable(file = "C:Data//For submission//For resubmission//RESUB2//First look//pstable.html", self_contained = T)


# Richness ----------------------------------------------------------------
gam.nb.richness.1<- gam(richness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.richness<- gam(richness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.richness.1, gam.poisson.richness)

plot(gam.poisson.richness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.poisson.richness)
#okay but qq plot not the best on ends
#qq_plot(gam.poisson.richness, method = 'simulate')
#k_check(gam.poisson.richness)
summary(gam.poisson.richness)
#a few outside the QQ plot on both ends
#low p value for k - but NS and edf is not super close to k-index

gam.poisson.richness.unordered<- gam(richness ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = poisson, select=TRUE, method="REML")
fam.gam.richness <- family(gam.poisson.richness)
ilink.gam.richness<- fam.gam.richness$linkinv
ilink.gam.richness


mod.richness<-gam.poisson.richness
ndata.richness <- with(invasion.exp.data_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.richness <- add_column(ndata.richness, fit = predict(mod.richness, newdata = ndata.richness, type = 'response'))

predict(mod.richness, newdata = ndata.richness, type = 'response')
ndata.richness <- bind_cols(ndata.richness, setNames(as_tibble(predict(mod.richness, ndata.richness, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.richness <- mutate(ndata.richness,
                               fit_resp  = ilink.gam.richness(fit_link),
                               right_upr = ilink.gam.richness(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.richness(fit_link - (2 * se_link)))


ndata.richness$min.10.pH.unscaled<-ndata.richness$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.richness <- ggplot(ndata.richness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = richness, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species richness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.richness,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,20)
plt.richness
ggsave("C:Graphs April 2020//richness_pred.png")







# Evenness ----------------------------------------------------------------

gam.lm.evenness<- gam(evenness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.evenness.1<- gam(evenness ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")

AICtab( gam.lm.evenness, gam.gamma.evenness.1)

plot(gam.lm.evenness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.evenness)
#looks very good
#qq_plot(gam.lm.evenness, method = 'simulate')
#k_check(gam.lm.evenness)
summary(gam.lm.evenness)

gam.lm.evenness.unordered<- gam(evenness ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")

fam.gam.evenness <- family(gam.lm.evenness)
ilink.gam.evenness<- fam.gam.evenness$linkinv
mod.evenness<-gam.lm.evenness
ndata.evenness <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.evenness <- add_column(ndata.evenness, fit = predict(mod.evenness, newdata = ndata.evenness, type = 'response'))

predict(mod.evenness, newdata = ndata.evenness, type = 'response')
ndata.evenness <- bind_cols(ndata.evenness, setNames(as_tibble(predict(mod.evenness, ndata.evenness, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.evenness <- mutate(ndata.evenness,
                         fit_resp  = ilink.gam.evenness(fit_link),
                         right_upr = ilink.gam.evenness(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.evenness(fit_link - (2 * se_link)))


ndata.evenness$min.10.pH.unscaled<-ndata.evenness$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.evenness <- ggplot(ndata.evenness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = evenness, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species evenness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.evenness
ggsave("C:Graphs April 2020//evenness_pred.png")


# Occupied space ----------------------------------------------------------

gam.lm.occupied.space<- gam(occupied.space ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.occupied.space.1<- gam(occupied.space*0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.binomial.occupied.space<- gam(formula = cbind(occupied.space, 100-occupied.space)~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = binomial, select=TRUE, method="REML")

gam.beta.occupied.space<- gam(occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.occupied.space.1<- gam(occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.occupied.space.2<- gam(occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.occupied.space.3<- gam(occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.lm.occupied.space, gam.beta.occupied.space.3, gam.gamma.occupied.space.1, gam.beta.occupied.space, gam.beta.occupied.space.1, gam.beta.occupied.space.2, gam.binomial.occupied.space)
#beta cauchit is best 

plot(gam.beta.occupied.space.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.occupied.space.3)
#a bit blocky
#qq_plot(gam.beta.occupied.space.3, method = 'simulate')
#k_check(gam.beta.occupied.space.3)
#k gettinga bit low but ns
summary(gam.beta.occupied.space.3)
gam.beta.occupied.space.3.unordered<- gam(occupied.space.001~ s(min.10.pH, k=15)+ Invasives + s(min.10.pH, by=oInvasives, k=15), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.occupied.space.3.unordered)

fam.gam.occupied.space <- family(gam.beta.occupied.space.3)
fam.gam.occupied.space
str(fam.gam.occupied.space)
ilink.gam.occupied.space<- fam.gam.occupied.space$linkinv
ilink.gam.occupied.space


mod.occupied.space<-gam.beta.occupied.space.3
ndata.occupied.space <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.occupied.space <- add_column(ndata.occupied.space, fit = predict(mod.occupied.space, newdata = ndata.occupied.space, type = 'response'))

predict(mod.occupied.space, newdata = ndata.occupied.space, type = 'response')
ndata.occupied.space <- bind_cols(ndata.occupied.space, setNames(as_tibble(predict(mod.occupied.space, ndata.occupied.space, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.occupied.space <- mutate(ndata.occupied.space,
                         fit_resp  = ilink.gam.occupied.space(fit_link),
                         right_upr = ilink.gam.occupied.space(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.occupied.space(fit_link - (2 * se_link)))


ndata.occupied.space$min.10.pH.unscaled<-ndata.occupied.space$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.occupied.space <- ggplot(ndata.occupied.space, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = occupied.space.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.occupied.space
ggsave("C:Graphs April 2020//occupied.space_pred.png")



# Everything wet weight ---------------------------------------------------

gam.lm.everything.wet.weight<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.1, gam.lm.log.everything.wet.weight, gam.tweedie.everything.wet.weight.1,  gam.lm.everything.wet.weight, gam.gamma.everything.wet.weight.1)

#gam.lm.log.everything.wet.weight is best 


plot(gam.lm.log.everything.wet.weight, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.everything.wet.weight)
#Looks quite good
#qq_plot(gam.lm.log.everything.wet.weight, method = 'simulate')
#k_check(gam.lm.log.everything.wet.weight)
summary(gam.lm.log.everything.wet.weight)

gam.lm.log.everything.wet.weight.unordered<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")

fam.gam.everything.wet.weight <- family(gam.lm.log.everything.wet.weight)
fam.gam.everything.wet.weight
ilink.gam.everything.wet.weight<- fam.gam.everything.wet.weight$linkinv
ilink.gam.everything.wet.weight


mod.everything.wet.weight<-gam.lm.log.everything.wet.weight
ndata.everything.wet.weight <- with(invasion.exp.data_zscores, 
                                    data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                    length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.everything.wet.weight <- add_column(ndata.everything.wet.weight, fit = predict(mod.everything.wet.weight, newdata = ndata.everything.wet.weight, type = 'response'))

predict(mod.everything.wet.weight, newdata = ndata.everything.wet.weight, type = 'response')
ndata.everything.wet.weight <- bind_cols(ndata.everything.wet.weight, setNames(as_tibble(predict(mod.everything.wet.weight, ndata.everything.wet.weight, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.everything.wet.weight <- mutate(ndata.everything.wet.weight,
                         fit_resp  = ilink.gam.everything.wet.weight(fit_link),
                         right_upr = ilink.gam.everything.wet.weight(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.everything.wet.weight(fit_link - (2 * se_link)))


ndata.everything.wet.weight$min.10.pH.unscaled<-ndata.everything.wet.weight$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight <- ggplot(ndata.everything.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(everything.wet.weight), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total wet biomass per mesocosm\n(g, log scale)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.everything.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight
ggsave("C:Graphs April 2020//everything.wet.weight_pred.png")


# Everything wet weight per 1 ---------------------------------------------

gam.lm.everything.wet.weight.per.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight.per.1<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.per.1.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.per.1.1, gam.lm.log.everything.wet.weight.per.1, gam.tweedie.everything.wet.weight.per.1.1,gam.lm.everything.wet.weight.per.1, gam.gamma.everything.wet.weight.per.1.1)

#gam.lm.log.everything.wet.weight.per.1 is best by far


plot(gam.lm.log.everything.wet.weight.per.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.everything.wet.weight.per.1)
#good but mabe a bit funnelly
#qq_plot(gam.lm.log.everything.wet.weight.per.1, method = 'simulate')
#k_check(gam.lm.log.everything.wet.weight.per.1)
summary(gam.lm.log.everything.wet.weight.per.1)

gam.lm.log.everything.wet.weight.per.1.unordered<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")

fam.gam.everything.wet.weight.per.1 <- family(gam.lm.log.everything.wet.weight.per.1)
fam.gam.everything.wet.weight.per.1
ilink.gam.everything.wet.weight.per.1<- fam.gam.everything.wet.weight.per.1$linkinv
ilink.gam.everything.wet.weight.per.1


mod.everything.wet.weight.per.1<-gam.lm.log.everything.wet.weight.per.1
ndata.everything.wet.weight.per.1 <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.everything.wet.weight.per.1 <- add_column(ndata.everything.wet.weight.per.1, fit = predict(mod.everything.wet.weight.per.1, newdata = ndata.everything.wet.weight.per.1, type = 'response'))

predict(mod.everything.wet.weight.per.1, newdata = ndata.everything.wet.weight.per.1, type = 'response')
ndata.everything.wet.weight.per.1 <- bind_cols(ndata.everything.wet.weight.per.1, setNames(as_tibble(predict(mod.everything.wet.weight.per.1, ndata.everything.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.everything.wet.weight.per.1 <- mutate(ndata.everything.wet.weight.per.1,
                                      fit_resp  = ilink.gam.everything.wet.weight.per.1(fit_link),
                                      right_upr = ilink.gam.everything.wet.weight.per.1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.everything.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.everything.wet.weight.per.1$min.10.pH.unscaled<-ndata.everything.wet.weight.per.1$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight.per.1 <- ggplot(ndata.everything.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(everything.wet.weight.per.1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Log Biomass per 1 % cover \n(wet weight)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.everything.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight.per.1
ggsave("C:Graphs April 2020//everything.wet.weight.per.1_pred.png")




# Dry weight --------------------------------------------------------------
gam.lm.total_dry_biomass<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass<- gam(log(total_dry_biomass) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass.1, gam.lm.log.total_dry_biomass, gam.tweedie.total_dry_biomass.1, gam.lm.total_dry_biomass, gam.gamma.total_dry_biomass.1)

#gam.lm.loglink is best but lm is only 0.1 so going with simplest.... 


plot(gam.lm.total_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.total_dry_biomass )
#Looks v good
#qq_plot(gam.lm.total_dry_biomass , method = 'simulate')
#k_check(gam.lm.total_dry_biomass )
summary(gam.lm.total_dry_biomass )

gam.lm.total_dry_biomass.unordered<- gam(total_dry_biomass ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
summary(gam.lm.total_dry_biomass.unordered)

fam.gam.total_dry_biomass <- family(gam.lm.total_dry_biomass)
fam.gam.total_dry_biomass
str(fam.gam.total_dry_biomass)
ilink.gam.total_dry_biomass<- fam.gam.total_dry_biomass$linkinv
ilink.gam.total_dry_biomass


mod.total_dry_biomass<-gam.lm.total_dry_biomass
ndata.total_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.total_dry_biomass <- add_column(ndata.total_dry_biomass, fit = predict(mod.total_dry_biomass, newdata = ndata.total_dry_biomass, type = 'response'))

predict(mod.total_dry_biomass, newdata = ndata.total_dry_biomass, type = 'response')
ndata.total_dry_biomass <- bind_cols(ndata.total_dry_biomass, setNames(as_tibble(predict(mod.total_dry_biomass, ndata.total_dry_biomass, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.total_dry_biomass <- mutate(ndata.total_dry_biomass,
                                      fit_resp  = ilink.gam.total_dry_biomass(fit_link),
                                      right_upr = ilink.gam.total_dry_biomass(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.total_dry_biomass(fit_link - (2 * se_link)))


ndata.total_dry_biomass$min.10.pH.unscaled<-ndata.total_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass <- ggplot(ndata.total_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (total_dry_biomass), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.total_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//total_dry_biomass_pred.png")


# Dry weight per 1 % cover --------------------------------------------------------------

gam.lm.total_dry_biomass_per1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass_per1.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass_per1<- gam(log(total_dry_biomass_per1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass_per1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass_per1.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass_per1.1, gam.lm.log.total_dry_biomass_per1, gam.tweedie.total_dry_biomass_per1, gam.lm.total_dry_biomass_per1, gam.gamma.total_dry_biomass_per1.1)

#tweedie is best by 1.7

plot(gam.tweedie.total_dry_biomass_per1 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.tweedie.total_dry_biomass_per1 )
#not available for tweedie
#qq_plot(gam.tweedie.total_dry_biomass_per1 , method = 'simulate')
#looks good
#k_check(gam.tweedie.total_dry_biomass_per1 )
summary(gam.tweedie.total_dry_biomass_per1 )

gam.tweedie.total_dry_biomass_per1.unordered<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=Invasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
summary(gam.tweedie.total_dry_biomass_per1.unordered)

fam.gam.total_dry_biomass_per1 <- family(gam.tweedie.total_dry_biomass_per1)
fam.gam.total_dry_biomass_per1
str(fam.gam.total_dry_biomass_per1)
ilink.gam.total_dry_biomass_per1<- fam.gam.total_dry_biomass_per1$linkinv
ilink.gam.total_dry_biomass_per1


mod.total_dry_biomass_per1<-gam.tweedie.total_dry_biomass_per1
ndata.total_dry_biomass_per1 <- with(invasion.exp.data_zscores, 
                                     data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.total_dry_biomass_per1 <- add_column(ndata.total_dry_biomass_per1, fit = predict(mod.total_dry_biomass_per1, newdata = ndata.total_dry_biomass_per1, type = 'response'))

predict(mod.total_dry_biomass_per1, newdata = ndata.total_dry_biomass_per1, type = 'response')
ndata.total_dry_biomass_per1 <- bind_cols(ndata.total_dry_biomass_per1, setNames(as_tibble(predict(mod.total_dry_biomass_per1, ndata.total_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.total_dry_biomass_per1 <- mutate(ndata.total_dry_biomass_per1,
                                  fit_resp  = ilink.gam.total_dry_biomass_per1(fit_link),
                                  right_upr = ilink.gam.total_dry_biomass_per1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.total_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.total_dry_biomass_per1$min.10.pH.unscaled<-ndata.total_dry_biomass_per1$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass_per1 <- ggplot(ndata.total_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (total_dry_biomass_per1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.total_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass_per1
ggplot2::ggsave("C:Graphs April 2020//total_dry_biomass_per1_pred.png")

# hydroid biomass ---------------------------------------------------------

gam.lm.hydroid_dry_biomass<- gam(hydroid_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.hydroid_dry_biomass<- gam(log(hydroid_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.hydroid_dry_biomass<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.loglink.hydroid_dry_biomass, gam.lm.log.hydroid_dry_biomass, gam.tweedie.hydroid_dry_biomass,  gam.lm.hydroid_dry_biomass, gam.gamma.hydroid_dry_biomass)

#gamma is best 
plot(gam.gamma.hydroid_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.hydroid_dry_biomass )
#looks good, maybe slightly funnelly
#qq_plot(gam.gamma.hydroid_dry_biomass , method = 'simulate')
#k_check(gam.gamma.hydroid_dry_biomass )
summary(gam.gamma.hydroid_dry_biomass )

gam.gamma.hydroid_dry_biomass.unordered<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.hydroid_dry_biomass <- family(gam.gamma.hydroid_dry_biomass)
fam.gam.hydroid_dry_biomass
str(fam.gam.hydroid_dry_biomass)
ilink.gam.hydroid_dry_biomass<- fam.gam.hydroid_dry_biomass$linkinv
ilink.gam.hydroid_dry_biomass


mod.hydroid_dry_biomass<-gam.gamma.hydroid_dry_biomass
ndata.hydroid_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydroid_dry_biomass <- add_column(ndata.hydroid_dry_biomass, fit = predict(mod.hydroid_dry_biomass, newdata = ndata.hydroid_dry_biomass, type = 'response'))

predict(mod.hydroid_dry_biomass, newdata = ndata.hydroid_dry_biomass, type = 'response')
ndata.hydroid_dry_biomass <- bind_cols(ndata.hydroid_dry_biomass, setNames(as_tibble(predict(mod.hydroid_dry_biomass, ndata.hydroid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydroid_dry_biomass <- mutate(ndata.hydroid_dry_biomass,
                                  fit_resp  = ilink.gam.hydroid_dry_biomass(fit_link),
                                  right_upr = ilink.gam.hydroid_dry_biomass(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.hydroid_dry_biomass(fit_link - (2 * se_link)))


ndata.hydroid_dry_biomass$min.10.pH.unscaled<-ndata.hydroid_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.hydroid_dry_biomass <- ggplot(ndata.hydroid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (hydroid_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Obelia") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.hydroid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.hydroid_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//hydroid_dry_biomass_pred.png")




# botryllus biomass -------------------------------------------------------
invasion.exp.data_zscores$tunicate_dry_biomass[invasion.exp.data_zscores$tunicate_dry_biomass<0]<-0

gam.lm.tunicate_dry_biomass<- gam(tunicate_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.tunicate_dry_biomass<- gam(log(tunicate_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.tunicate_dry_biomass<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.tunicate_dry_biomass, gam.lm.log.tunicate_dry_biomass, gam.tweedie.tunicate_dry_biomass,gam.lm.tunicate_dry_biomass, gam.gamma.tunicate_dry_biomass)

#gamma is the best

plot(gam.gamma.tunicate_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.tunicate_dry_biomass )
#a bit patterny
#qq_plot(gam.gamma.tunicate_dry_biomass , method = 'simulate')
#k_check(gam.gamma.tunicate_dry_biomass )
summary(gam.gamma.tunicate_dry_biomass )

gam.gamma.tunicate_dry_biomass.unordered<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.tunicate_dry_biomass <- family(gam.gamma.tunicate_dry_biomass)
fam.gam.tunicate_dry_biomass
str(fam.gam.tunicate_dry_biomass)
ilink.gam.tunicate_dry_biomass<- fam.gam.tunicate_dry_biomass$linkinv
ilink.gam.tunicate_dry_biomass


mod.tunicate_dry_biomass<-gam.gamma.tunicate_dry_biomass
ndata.tunicate_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.tunicate_dry_biomass <- add_column(ndata.tunicate_dry_biomass, fit = predict(mod.tunicate_dry_biomass, newdata = ndata.tunicate_dry_biomass, type = 'response'))

predict(mod.tunicate_dry_biomass, newdata = ndata.tunicate_dry_biomass, type = 'response')
ndata.tunicate_dry_biomass <- bind_cols(ndata.tunicate_dry_biomass, setNames(as_tibble(predict(mod.tunicate_dry_biomass, ndata.tunicate_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.tunicate_dry_biomass <- mutate(ndata.tunicate_dry_biomass,
                                    fit_resp  = ilink.gam.tunicate_dry_biomass(fit_link),
                                    right_upr = ilink.gam.tunicate_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.tunicate_dry_biomass(fit_link - (2 * se_link)))


ndata.tunicate_dry_biomass$min.10.pH.unscaled<-ndata.tunicate_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.tunicate_dry_biomass <- ggplot(ndata.tunicate_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (tunicate_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.tunicate_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.tunicate_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//tunicate_dry_biomass_pred.png")



# caprellid biomass -------------------------------------------------------
gam.lm.caprellid_dry_biomass<- gam(caprellid_dry_biomass ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass<- gam(log(caprellid_dry_biomass+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass, gam.lm.log.caprellid_dry_biomass, gam.tweedie.caprellid_dry_biomass, gam.lm.caprellid_dry_biomass, gam.gamma.caprellid_dry_biomass)
#gamma by 1.2

plot(gam.gamma.caprellid_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.caprellid_dry_biomass )
# look good
#qq_plot(gam.gamma.caprellid_dry_biomass , method = 'simulate')
#k_check(gam.gamma.caprellid_dry_biomass )
summary(gam.gamma.caprellid_dry_biomass )
gam.gamma.caprellid_dry_biomass.unordered<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.caprellid_dry_biomass <- family(gam.gamma.caprellid_dry_biomass)
ilink.gam.caprellid_dry_biomass<- fam.gam.caprellid_dry_biomass$linkinv
mod.caprellid_dry_biomass<-gam.gamma.caprellid_dry_biomass
ndata.caprellid_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.caprellid_dry_biomass <- add_column(ndata.caprellid_dry_biomass, fit = predict(mod.caprellid_dry_biomass, newdata = ndata.caprellid_dry_biomass, type = 'response'))

predict(mod.caprellid_dry_biomass, newdata = ndata.caprellid_dry_biomass, type = 'response')
ndata.caprellid_dry_biomass <- bind_cols(ndata.caprellid_dry_biomass, setNames(as_tibble(predict(mod.caprellid_dry_biomass, ndata.caprellid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid_dry_biomass <- mutate(ndata.caprellid_dry_biomass,
                                    fit_resp  = ilink.gam.caprellid_dry_biomass(fit_link),
                                    right_upr = ilink.gam.caprellid_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.caprellid_dry_biomass(fit_link - (2 * se_link)))


ndata.caprellid_dry_biomass$min.10.pH.unscaled<-ndata.caprellid_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass <- ggplot(ndata.caprellid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (caprellid_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//caprellid_dry_biomass_pred.png")



# caprellid biomass per invididual -------------------------------------------------------

invasion.exp.data_zscores$caprellid_dry_biomass_per1<-invasion.exp.data_zscores$caprellid_dry_biomass/(food.caprellid.data_zscores$total.caprellids)

gam.lm.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass_per1<- gam(log(caprellid_dry_biomass_per1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass_per1<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass_per1, gam.lm.log.caprellid_dry_biomass_per1, gam.tweedie.caprellid_dry_biomass_per1,  gam.lm.caprellid_dry_biomass_per1, gam.gamma.caprellid_dry_biomass_per1)

#gamma is the best by 2.2
plot(gam.gamma.caprellid_dry_biomass_per1 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.caprellid_dry_biomass_per1 )
#a bit funnelly and qq right tail
#qq_plot(gam.gamma.caprellid_dry_biomass_per1 , method = 'simulate')
#k_check(gam.gamma.caprellid_dry_biomass_per1 )
summary(gam.gamma.caprellid_dry_biomass_per1 )
gam.gamma.caprellid_dry_biomass_per1.unordered<- gam(caprellid_dry_biomass_per1+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.caprellid_dry_biomass_per1 <- family(gam.gamma.caprellid_dry_biomass_per1)
fam.gam.caprellid_dry_biomass_per1
str(fam.gam.caprellid_dry_biomass_per1)
ilink.gam.caprellid_dry_biomass_per1<- fam.gam.caprellid_dry_biomass_per1$linkinv
ilink.gam.caprellid_dry_biomass_per1


mod.caprellid_dry_biomass_per1<-gam.gamma.caprellid_dry_biomass_per1
ndata.caprellid_dry_biomass_per1 <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.caprellid_dry_biomass_per1 <- add_column(ndata.caprellid_dry_biomass_per1, fit = predict(mod.caprellid_dry_biomass_per1, newdata = ndata.caprellid_dry_biomass_per1, type = 'response'))

predict(mod.caprellid_dry_biomass_per1, newdata = ndata.caprellid_dry_biomass_per1, type = 'response')
ndata.caprellid_dry_biomass_per1 <- bind_cols(ndata.caprellid_dry_biomass_per1, setNames(as_tibble(predict(mod.caprellid_dry_biomass_per1, ndata.caprellid_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid_dry_biomass_per1 <- mutate(ndata.caprellid_dry_biomass_per1,
                                      fit_resp  = ilink.gam.caprellid_dry_biomass_per1(fit_link),
                                      right_upr = ilink.gam.caprellid_dry_biomass_per1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.caprellid_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.caprellid_dry_biomass_per1$min.10.pH.unscaled<-ndata.caprellid_dry_biomass_per1$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass_per1 <- ggplot(ndata.caprellid_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (caprellid_dry_biomass_per1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,0.012)
plt.caprellid_dry_biomass_per1
ggplot2::ggsave("C:Graphs April 2020//caprellid_dry_biomass_per1_pred.png")


# rest biomass ------------------------------------------------------------

#kcheck was significant so increase k from 10 to 11
gam.lm.rest_dry_biomass<- gam(rest_dry_biomass ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.rest_dry_biomass<- gam(log(rest_dry_biomass+0.1) ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.rest_dry_biomass<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oInvasives + s(min.10.pH, by=oInvasives, k=12),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.rest_dry_biomass, gam.lm.log.rest_dry_biomass, gam.tweedie.rest_dry_biomass, gam.lm.rest_dry_biomass, gam.gamma.rest_dry_biomass)

#tweedie the best 


plot(gam.tweedie.rest_dry_biomass , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#qq_plot(gam.tweedie.rest_dry_biomass , method = 'simulate')
#k_check(gam.tweedie.rest_dry_biomass )
summary(gam.tweedie.rest_dry_biomass)

gam.tweedie.rest_dry_biomass.unordered<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
summary(gam.tweedie.rest_dry_biomass.unordered)

fam.gam.rest_dry_biomass <- family(gam.tweedie.rest_dry_biomass)
fam.gam.rest_dry_biomass
str(fam.gam.rest_dry_biomass)
ilink.gam.rest_dry_biomass<- fam.gam.rest_dry_biomass$linkinv
ilink.gam.rest_dry_biomass


mod.rest_dry_biomass<-gam.tweedie.rest_dry_biomass
ndata.rest_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))
## add the fitted values by predicting from the model for the new data
ndata.rest_dry_biomass <- add_column(ndata.rest_dry_biomass, fit = predict(mod.rest_dry_biomass, newdata = ndata.rest_dry_biomass, type = 'response'))

predict(mod.rest_dry_biomass, newdata = ndata.rest_dry_biomass, type = 'response')
ndata.rest_dry_biomass <- bind_cols(ndata.rest_dry_biomass, setNames(as_tibble(predict(mod.rest_dry_biomass, ndata.rest_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.rest_dry_biomass <- mutate(ndata.rest_dry_biomass,
                                    fit_resp  = ilink.gam.rest_dry_biomass(fit_link),
                                    right_upr = ilink.gam.rest_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.rest_dry_biomass(fit_link - (2 * se_link)))


ndata.rest_dry_biomass$min.10.pH.unscaled<-ndata.rest_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.rest_dry_biomass <- ggplot(ndata.rest_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = (rest_dry_biomass+0.01), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression("Remaining dry weight per tile (g)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.rest_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.rest_dry_biomass
ggplot2::ggsave("C:Graphs April 2020//rest_dry_biomass_pred.png")




# Mussel wet weight -------------------------------------------------------

gam.lm.Mussel.wet.weight<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.1<- gam(Mussel.wet.weight+0.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.Mussel.wet.weight.1, gam.lm.log.Mussel.wet.weight, gam.tweedie.Mussel.wet.weight.1,gam.lm.Mussel.wet.weight, gam.gamma.Mussel.wet.weight.1)

#gam.lm.log.Mussel.wet.weight is best 

plot(gam.lm.log.Mussel.wet.weight, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.Mussel.wet.weight)
#look pretty good, esp qq
#qq_plot(gam.lm.log.Mussel.wet.weight, method = 'simulate')
#k_check(gam.lm.log.Mussel.wet.weight)
summary(gam.lm.log.Mussel.wet.weight)

#residuals many in a straight line but otherwise good 

gam.lm.log.Mussel.wet.weight.unordered<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
fam.gam.Mussel.wet.weight <- family(gam.lm.log.Mussel.wet.weight)
fam.gam.Mussel.wet.weight
ilink.gam.Mussel.wet.weight<- fam.gam.Mussel.wet.weight$linkinv
ilink.gam.Mussel.wet.weight


mod.Mussel.wet.weight<-gam.lm.log.Mussel.wet.weight
ndata.Mussel.wet.weight <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.Mussel.wet.weight <- add_column(ndata.Mussel.wet.weight, fit = predict(mod.Mussel.wet.weight, newdata = ndata.Mussel.wet.weight, type = 'response'))

predict(mod.Mussel.wet.weight, newdata = ndata.Mussel.wet.weight, type = 'response')
ndata.Mussel.wet.weight <- bind_cols(ndata.Mussel.wet.weight, setNames(as_tibble(predict(mod.Mussel.wet.weight, ndata.Mussel.wet.weight, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.Mussel.wet.weight <- mutate(ndata.Mussel.wet.weight,
                                      fit_resp  = ilink.gam.Mussel.wet.weight(fit_link),
                                      right_upr = ilink.gam.Mussel.wet.weight(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.Mussel.wet.weight(fit_link - (2 * se_link)))


ndata.Mussel.wet.weight$min.10.pH.unscaled<-ndata.Mussel.wet.weight$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight <- ggplot(ndata.Mussel.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = log(Mussel.wet.weight+0.1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Mytilus") ~"wet weight (g, log scale)"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.Mussel.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight
ggsave("C:Graphs April 2020//Mussel.wet.weight_pred.png")






# Mussel wet weight per individual mussel ---------------------------------
#big outlier .... 5 mussels weigh 105.5 g?? so 17 g per mussel? 

invasion.exp.data_zscores$Mussel.wet.weight[invasion.exp.data_zscores$Mesocosm==8]<-14.0
invasion.exp.data_zscores$Mussel.wet.weight.per.1<-(invasion.exp.data_zscores$Mussel.wet.weight)/(invasion.exp.data_zscores$mussel+1)

gam.lm.Mussel.wet.weight.per.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.per.1<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight.per.1<- gam(log(Mussel.wet.weight.per.1+0.01) ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.per.1.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.per.1.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.loglink.Mussel.wet.weight.per.1.1, gam.lm.log.Mussel.wet.weight.per.1, gam.tweedie.Mussel.wet.weight.per.1.1,  gam.lm.Mussel.wet.weight.per.1, gam.gamma.Mussel.wet.weight.per.1)

#gamma is best 
plot(gam.gamma.Mussel.wet.weight.per.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.Mussel.wet.weight.per.1)
# good but a bit funnelly
#qq_plot(gam.gamma.Mussel.wet.weight.per.1, method = 'simulate')
#k_check(gam.gamma.Mussel.wet.weight.per.1)
summary(gam.gamma.Mussel.wet.weight.per.1)

gam.gamma.Mussel.wet.weight.per.1.unordered<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.Mussel.wet.weight.per.1 <- family(gam.gamma.Mussel.wet.weight.per.1)
fam.gam.Mussel.wet.weight.per.1
str(fam.gam.Mussel.wet.weight.per.1)
ilink.gam.Mussel.wet.weight.per.1<- fam.gam.Mussel.wet.weight.per.1$linkinv
ilink.gam.Mussel.wet.weight.per.1


mod.Mussel.wet.weight.per.1<-gam.gamma.Mussel.wet.weight.per.1
ndata.Mussel.wet.weight.per.1 <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.Mussel.wet.weight.per.1 <- add_column(ndata.Mussel.wet.weight.per.1, fit = predict(mod.Mussel.wet.weight.per.1, newdata = ndata.Mussel.wet.weight.per.1, type = 'response'))

predict(mod.Mussel.wet.weight.per.1, newdata = ndata.Mussel.wet.weight.per.1, type = 'response')
ndata.Mussel.wet.weight.per.1 <- bind_cols(ndata.Mussel.wet.weight.per.1, setNames(as_tibble(predict(mod.Mussel.wet.weight.per.1, ndata.Mussel.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.Mussel.wet.weight.per.1 <- mutate(ndata.Mussel.wet.weight.per.1,
                                  fit_resp  = ilink.gam.Mussel.wet.weight.per.1(fit_link),
                                  right_upr = ilink.gam.Mussel.wet.weight.per.1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.Mussel.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.Mussel.wet.weight.per.1$min.10.pH.unscaled<-ndata.Mussel.wet.weight.per.1$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight.per.1 <- ggplot(ndata.Mussel.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = Mussel.wet.weight.per.1, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Biomass per mussel \n(wet weight)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.Mussel.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight.per.1
ggsave("C:Graphs April 2020//Mussel.wet.weight.per.1_pred.png")




# Hydtobot ----------------------------------------------------------------
invasion.exp.data_zscores$hydtobot[invasion.exp.data_zscores$hydtobot==1]<-0.99
invasion.exp.data_zscores$hydtobot[invasion.exp.data_zscores$hydtobot==0]<-0.01

gam.lm.hydtobot<- gam(log(hydtobot+0.01)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.tweedie.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family =tw(), select=TRUE, method="REML")

gam.beta.hydtobot<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot.1<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot.2<- gam(hydtobot~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")


AICtab( gam.tweedie.hydtobot ,gam.beta.hydtobot, gam.lm.hydtobot, gam.log.lm.hydtobot, gam.beta.hydtobot.1, gam.beta.hydtobot.2)
### beta logit is the best 

plot(gam.beta.hydtobot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.hydtobot)
#qq a bit wiggly
#qq_plot(gam.beta.hydtobot, method = 'simulate')
#k_check(gam.beta.hydtobot)
gam.check(gam.tweedie.hydtobot)
summary(gam.beta.hydtobot)
#not the best qq plot
#resids are in rows

gam.beta.hydtobot.unordered<- gam(hydtobot~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
fam.gam.hydtobot <- family(gam.beta.hydtobot)
fam.gam.hydtobot 

ilink.gam.hydtobot <- fam.gam.hydtobot$linkinv
ilink.gam.hydtobot


invasion.exp.data_zscores$min.10.pH.unscaled <-invasion.exp.data_zscores$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data_zscores)

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

mod.hydtobot<-gam.beta.hydtobot
ndata.hydtobot <- with(invasion.exp.data_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                               length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydtobot <- add_column(ndata.hydtobot, fit = predict(mod.hydtobot, newdata = ndata.hydtobot, type = 'response'))


ndata.hydtobot <- bind_cols(ndata.hydtobot, setNames(as_tibble(predict(mod.hydtobot, ndata.hydtobot, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydtobot <- mutate(ndata.hydtobot,
                        fit_resp  = ilink.gam.hydtobot(fit_link),
                        right_upr = ilink.gam.hydtobot(fit_link + (2 * se_link)),
                        right_lwr = ilink.gam.hydtobot(fit_link - (2 * se_link)))




ndata.hydtobot$min.10.pH.unscaled<-ndata.hydtobot$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.hydtobot <- ggplot(ndata.hydtobot, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydtobot, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "cover ratio"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.hydtobot,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.2, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=6), legend.title = element_text(size=7))+ 
  geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot 


# Hydtobot by weight ------------------------------------------------------
invasion.exp.data_zscores$hydtobot_dry_biomass[invasion.exp.data_zscores$hydtobot_dry_biomass==1]<-0.99
invasion.exp.data_zscores$hydtobot_dry_biomass[invasion.exp.data_zscores$hydtobot_dry_biomass==0]<-0.01

gam.lm.hydtobot_dry_biomass<- gam(log(hydtobot_dry_biomass+0.01)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.tweedie.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family =tw(), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.1<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.2<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")

AICtab(gam.tweedie.hydtobot_dry_biomass ,gam.beta.hydtobot_dry_biomass, gam.lm.hydtobot_dry_biomass, gam.log.lm.hydtobot_dry_biomass, gam.beta.hydtobot_dry_biomass.1, gam.beta.hydtobot_dry_biomass.2, lm.hydtobot_dry_biomass, lm.log.hydtobot_dry_biomass)
### logit is the best 

plot(gam.beta.hydtobot_dry_biomass, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.hydtobot_dry_biomass)
#look ok, qq is bloacky
#qq_plot(gam.beta.hydtobot_dry_biomass, method = 'simulate')
#k_check(gam.beta.hydtobot_dry_biomass)
summary(gam.beta.hydtobot_dry_biomass)
gam.beta.hydtobot_dry_biomass.unordered<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.hydtobot_dry_biomass <- family(gam.beta.hydtobot_dry_biomass)
fam.gam.hydtobot_dry_biomass 
str(fam.gam.hydtobot_dry_biomass )
ilink.gam.hydtobot_dry_biomass <- fam.gam.hydtobot_dry_biomass$linkinv
ilink.gam.hydtobot_dry_biomass

invasion.exp.data_zscores$min.10.pH.unscaled <-invasion.exp.data_zscores$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')
head(invasion.exp.data_zscores)

want <- seq(1, nrow(invasion.exp.data_zscores), length.out = 100)

mod.hydtobot_dry_biomass<-gam.beta.hydtobot_dry_biomass
ndata.hydtobot_dry_biomass <- with(invasion.exp.data_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydtobot_dry_biomass <- add_column(ndata.hydtobot_dry_biomass, fit = predict(mod.hydtobot_dry_biomass, newdata = ndata.hydtobot_dry_biomass, type = 'response'))


ndata.hydtobot_dry_biomass <- bind_cols(ndata.hydtobot_dry_biomass, setNames(as_tibble(predict(mod.hydtobot_dry_biomass, ndata.hydtobot_dry_biomass, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydtobot_dry_biomass <- mutate(ndata.hydtobot_dry_biomass,
                         fit_resp  = ilink.gam.hydtobot_dry_biomass(fit_link),
                         right_upr = ilink.gam.hydtobot_dry_biomass(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.hydtobot_dry_biomass(fit_link - (2 * se_link)))




ndata.hydtobot_dry_biomass$min.10.pH.unscaled<-ndata.hydtobot_dry_biomass$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.hydtobot_dry_biomass <- ggplot(ndata.hydtobot_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = hydtobot_dry_biomass, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "biomass ratio"))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.hydtobot_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)#+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot_dry_biomass 


#mesocosm #8 didn't have either tunicates or hydroids weight - 0% hydroids, 2% tunicates


# CAP1 --------------------------------------------------------------------
#negative values so can't do gamma
gam.lm.CAP1<- gam(CAP1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.loglink.CAP1.1<- gam(CAP1 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( gam.loglink.CAP1.1, gam.lm.CAP1)
#gam.lm.CAP1

plot(gam.lm.CAP1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.CAP1)
#look good
#qq_plot(gam.lm.CAP1, method = 'simulate')
#k_check(gam.lm.CAP1)
summary(gam.lm.CAP1)
gam.lm.CAP1.unordered<- gam(CAP1 ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data_zscores, select=TRUE, method="REML")
summary(gam.lm.CAP1.unordered)

fam.gam.CAP1 <- family(gam.lm.CAP1)
fam.gam.CAP1
str(fam.gam.CAP1)
ilink.gam.CAP1<- fam.gam.CAP1$linkinv
ilink.gam.CAP1


mod.CAP1<-gam.lm.CAP1
ndata.CAP1 <- with(invasion.exp.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.CAP1 <- add_column(ndata.CAP1, fit = predict(mod.CAP1, newdata = ndata.CAP1, type = 'response'))

predict(mod.CAP1, newdata = ndata.CAP1, type = 'response')
ndata.CAP1 <- bind_cols(ndata.CAP1, setNames(as_tibble(predict(mod.CAP1, ndata.CAP1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.CAP1 <- mutate(ndata.CAP1,
                                      fit_resp  = ilink.gam.CAP1(fit_link),
                                      right_upr = ilink.gam.CAP1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.CAP1(fit_link - (2 * se_link)))


ndata.CAP1$min.10.pH.unscaled<-ndata.CAP1$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.CAP1 <- ggplot(ndata.CAP1, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(CAP1), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(36% of constrained variation)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.1, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.CAP1
ggsave("C:Graphs April 2020//CAP1_pred.png")




# Distances ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.lm.distances<- gam(distcentroid ~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data_zscores, select=TRUE, method="REML")
gam.loglink.distances.1<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.gamma.distances<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.loglink.distances.1,  gam.lm.distances, gam.gamma.distances)
#gam.lm.distances although both are equal


plot(gam.lm.distances, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.distances)
#looks good
#qq_plot(gam.lm.distances, method = 'simulate')
#k_check(gam.lm.distances)
#good now
summary(gam.lm.distances)

gam.lm.distances.unordered<- gam(distcentroid~ s(min.10.pH, k=11)+ Invasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data_zscores, select=TRUE, method="REML")
summary(gam.lm.distances.unordered)

fam.gam.distances <- family(gam.lm.distances)
fam.gam.distances
str(fam.gam.distances)
ilink.gam.distances<- fam.gam.distances$linkinv
ilink.gam.distances


mod.distances<-gam.lm.distances
ndata.distances <- with(invasion.exp.data_zscores, tibble(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                             length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.distances <- add_column(ndata.distances, fit = predict(mod.distances, newdata = ndata.distances, type = 'response'))

predict(mod.distances, newdata = ndata.distances, type = 'response')
ndata.distances <- bind_cols(ndata.distances, setNames(as_tibble(predict(mod.distances, ndata.distances, se.fit = TRUE)[1:2]),
                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.distances <- mutate(ndata.distances,
                     fit_resp  = ilink.gam.distances(fit_link),
                     right_upr = ilink.gam.distances(fit_link + (2 * se_link)),
                     right_lwr = ilink.gam.distances(fit_link - (2 * se_link)))


ndata.distances$min.10.pH.unscaled<-ndata.distances$min.10.pH * attr(invasion.exp.data_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data_zscores$min.10.pH, 'scaled:center')


# plot 
plt.distances <- ggplot(ndata.distances, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(distcentroid), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
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
fig.s3.hydtobot<-wrap_plots(plt.gam.hydtobot, plt.gam.hydtobot_dry_biomass, ncol=2) + plot_annotation(tag_levels = 'a')

fig.s3.hydtobot

ggplot2::ggsave("C:Data//For submission//For resubmission//RESUB2//First look//Fig.S3.hydotobot.png", width=6, height=3, units="in", dpi=600)



# Community level tables --------------------------------------------------

richness.gam<- summary(gam.poisson.richness)
richness.gam.unordered<- summary(gam.poisson.richness.unordered)

evenness.gam<-summary(gam.lm.evenness)
evenness.gam.unordered<-summary(gam.lm.evenness.unordered)


occupied.space.gam<- summary(gam.beta.occupied.space.3)
occupied.space.gam.unordered<- summary(gam.beta.occupied.space.3.unordered)

distances.gam<- summary(gam.lm.distances)
distances.gam.unordered<- summary(gam.lm.distances.unordered)


CAP1.gam <- summary(gam.lm.CAP1)
CAP1.gam.unordered <- summary(gam.lm.CAP1.unordered)

#dry biomass
total_dry_biomass.gam <- summary(gam.lm.total_dry_biomass)
total_dry_biomass.gam.unordered <- summary(gam.lm.total_dry_biomass.unordered)

hydroid_dry_biomass.gam.unordered <- summary(gam.gamma.hydroid_dry_biomass.unordered) 
hydroid_dry_biomass.gam<- summary(gam.gamma.hydroid_dry_biomass) 

caprellid_dry_biomass.gam.unordered <- summary(gam.gamma.caprellid_dry_biomass.unordered)
caprellid_dry_biomass.gam <- summary(gam.gamma.caprellid_dry_biomass)

tunicate_dry_biomass.gam <- summary(gam.gamma.tunicate_dry_biomass)
tunicate_dry_biomass.gam.unordered <- summary(gam.gamma.tunicate_dry_biomass.unordered)

rest_dry_biomass.gam.unordered <- summary(gam.tweedie.rest_dry_biomass.unordered)
rest_dry_biomass.gam <- summary(gam.tweedie.rest_dry_biomass)

hydtobot_dry_biomass.gam.unordered<-summary(gam.beta.hydtobot_dry_biomass.unordered)
hydtobot_dry_biomass.gam<-summary(gam.beta.hydtobot_dry_biomass)

#wet biomass
everything.wet.weight.gam <-summary(gam.lm.log.everything.wet.weight)
everything.wet.weight.gam.unordered <- summary(gam.lm.log.everything.wet.weight.unordered)

everything.wet.weight.per.1.gam <-summary(gam.lm.log.everything.wet.weight.per.1)
everything.wet.weight.per.1.gam.unordered <- summary(gam.lm.log.everything.wet.weight.per.1.unordered)


Mussel.wet.weight.gam <- summary(gam.lm.log.Mussel.wet.weight)
Mussel.wet.weight.gam.unordered <-summary(gam.lm.log.Mussel.wet.weight.unordered)

Mussel.wet.weight.per.1.gam <- summary(gam.gamma.Mussel.wet.weight.per.1)
Mussel.wet.weight.per.1.gam.unordered <-summary(gam.gamma.Mussel.wet.weight.per.1.unordered)


#competition metric 
hydtobot.gam <- summary(gam.beta.hydtobot)
hydtobot.gam.unordered <- summary(gam.beta.hydtobot.unordered)

#ptable building
richness.gam.p.table<-as.data.frame(richness.gam.unordered$p.table)
richness.gam.s.table<-as.data.frame(richness.gam$s.table)

evenness.gam.p.table<-as.data.frame(evenness.gam.unordered$p.table)
evenness.gam.s.table<-as.data.frame(evenness.gam$s.table)

occupied.space.gam.p.table<-as.data.frame(occupied.space.gam.unordered$p.table)
occupied.space.gam.s.table<-as.data.frame(occupied.space.gam$s.table)

distances.gam.p.table<-as.data.frame(distances.gam.unordered$p.table)
distances.gam.s.table<-as.data.frame(distances.gam$s.table)

CAP1.gam.p.table<-as.data.frame(CAP1.gam.unordered$p.table)
CAP1.gam.s.table<-as.data.frame(CAP1.gam$s.table)

total_dry_biomass.gam.p.table<-as.data.frame(total_dry_biomass.gam.unordered$p.table)
total_dry_biomass.gam.s.table<-as.data.frame(total_dry_biomass.gam$s.table)

hydroid_dry_biomass.gam.p.table<-as.data.frame(hydroid_dry_biomass.gam.unordered$p.table)
hydroid_dry_biomass.gam.s.table<-as.data.frame(hydroid_dry_biomass.gam$s.table)

caprellid_dry_biomass.gam.p.table<-as.data.frame(caprellid_dry_biomass.gam.unordered$p.table)
caprellid_dry_biomass.gam.s.table<-as.data.frame(caprellid_dry_biomass.gam$s.table)

tunicate_dry_biomass.gam.p.table<-as.data.frame(tunicate_dry_biomass.gam.unordered$p.table)
tunicate_dry_biomass.gam.s.table<-as.data.frame(tunicate_dry_biomass.gam$s.table)

rest_dry_biomass.gam.p.table<-as.data.frame(rest_dry_biomass.gam.unordered$p.table)
rest_dry_biomass.gam.s.table<-as.data.frame(rest_dry_biomass.gam$s.table)

everything.wet.weight.gam.p.table<-as.data.frame(everything.wet.weight.gam.unordered$p.table)
everything.wet.weight.gam.s.table<-as.data.frame(everything.wet.weight.gam$s.table)

everything.wet.weight.per.1.gam.p.table<-as.data.frame(everything.wet.weight.per.1.gam.unordered$p.table)
everything.wet.weight.per.1.gam.s.table<-as.data.frame(everything.wet.weight.per.1.gam$s.table)

Mussel.wet.weight.gam.p.table<-as.data.frame(Mussel.wet.weight.gam.unordered$p.table)
Mussel.wet.weight.gam.s.table<-as.data.frame(Mussel.wet.weight.gam$s.table)

Mussel.wet.weight.per.1.gam.p.table<-as.data.frame(Mussel.wet.weight.per.1.gam.unordered$p.table)
Mussel.wet.weight.per.1.gam.s.table<-as.data.frame(Mussel.wet.weight.per.1.gam$s.table)

hydtobot.gam.p.table<-as.data.frame(hydtobot.gam.unordered$p.table)
hydtobot.gam.s.table<-as.data.frame(hydtobot.gam$s.table)

hydtobot_dry_biomass.gam.p.table<-as.data.frame(hydtobot_dry_biomass.gam.unordered$p.table)
hydtobot_dry_biomass.gam.s.table<-as.data.frame(hydtobot_dry_biomass.gam$s.table)

hydtobot_dry_biomass.gam.p.table
hydtobot_dry_biomass.gam.s.table

#richness.gam.p.table and  hydtobot.gam.p.table, is with z value 
colnames(richness.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot_dry_biomass.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(occupied.space.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.community.t<-rbind(richness.gam.p.table,
              evenness.gam.p.table,
              occupied.space.gam.p.table,
              total_dry_biomass.gam.p.table,
              CAP1.gam.p.table,
              distances.gam.p.table,
              hydtobot.gam.p.table, 
              hydtobot_dry_biomass.gam.p.table,
              everything.wet.weight.gam.p.table
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
#richness.gam.p.table and  hydtobot.gam.p.table, is with Chisq
colnames(richness.gam.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(hydtobot.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")
colnames(hydtobot_dry_biomass.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")

colnames(occupied.space.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.community.f<-rbind(richness.gam.s.table,
                          evenness.gam.s.table,
                          occupied.space.gam.s.table,
                          total_dry_biomass.gam.s.table,
                          CAP1.gam.s.table,
                          distances.gam.s.table,
                          hydtobot.gam.s.table, 
                          hydtobot_dry_biomass.gam.s.table,
                          everything.wet.weight.gam.s.table
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
