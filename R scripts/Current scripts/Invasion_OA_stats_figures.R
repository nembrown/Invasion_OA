
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
library(patchwork)
library(codyn)


# dataframe management ----------------------------------------------------

invasion.exp.data<-read.csv("C:Biological data/allweeks_cover_counts_with_pres.csv",stringsAsFactors = FALSE, na.strings = c("NA","") )
head(invasion.exp.data)
#ordered and unordered factors
invasion.exp.data$oTreatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=TRUE)
invasion.exp.data$Treatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=FALSE)

invasion.exp.data$oInvasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=TRUE)
invasion.exp.data$Invasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=FALSE)

invasion.exp.data$oCO2.Treatment<-factor(invasion.exp.data$CO2.Treatment, levels=c("AIR", "CO2"), ordered=TRUE)
invasion.exp.data$CO2.Treatment<-factor(invasion.exp.data$CO2.Treatment, levels=c("AIR", "CO2"), ordered=FALSE)

#New variables
invasion.exp.data$num.nudi<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.eggs+invasion.exp.data$nudi.hatched
invasion.exp.data$hydroid.001<-(0.01*(invasion.exp.data$hydroid))+0.01
invasion.exp.data$botryllid.001<-(0.01*(invasion.exp.data$botryllid))+0.01
invasion.exp.data$bot.eaten.001<-(0.01*(invasion.exp.data$bot.eaten))+0.01
invasion.exp.data$mem.eaten.001<-(0.01*(invasion.exp.data$mem.eaten))+0.01

invasion.exp.data$membranipora.001<-(0.01*(invasion.exp.data$membranipora))+0.01
invasion.exp.data$mussel.001<-(0.01*(invasion.exp.data$mussel))+0.01
invasion.exp.data$didemnum<-invasion.exp.data$white.bryo
invasion.exp.data$folliculina<-invasion.exp.data$protozoa
invasion.exp.data$folliculina.001<-(0.01*(invasion.exp.data$folliculina))+0.01
invasion.exp.data$didemnum.001<-(0.01*(invasion.exp.data$didemnum))+0.01
invasion.exp.data$occupied.space<-(100 - invasion.exp.data$bare)
invasion.exp.data$occupied.space.001<-(0.01*(invasion.exp.data$occupied.space))+0.01
invasion.exp.data$native.occupied.space<-(100 - invasion.exp.data$botryllid  -invasion.exp.data$bot.eaten -  invasion.exp.data$bare)
invasion.exp.data$native.occupied.space.001<-(0.01*(invasion.exp.data$native.occupied.space))+0.01


head(invasion.exp.data)
#splitting into mid and end
invasion.exp.data.16<-invasion.exp.data %>% filter(Week==16)
invasion.exp.data.8<-invasion.exp.data %>% filter(Week==8)

  
# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
invasion.exp.data.16.community<-read.csv("C:Biological data/invasion.exp.data.16.community.csv")
invasion.exp.data.16<-merge(invasion.exp.data.16, invasion.exp.data.16.community, by="Tile.ID")
invasion.exp.data.16_zscores<-invasion.exp.data.16
#invasion.exp.data.16_zscores$hydrogen.concentration<-scale(invasion.exp.data.16$hydrogen.concentration, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$av.pH<-scale(invasion.exp.data.16$av.pH, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$min.10.pH<-scale(invasion.exp.data.16$min.10.pH, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$Mesocosm <- as.factor(invasion.exp.data.16$Mesocosm)
invasion.exp.data.16_zscores$av.pH.unscaled <-invasion.exp.data.16_zscores$av.pH * attr(invasion.exp.data.16_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$av.pH, 'scaled:center')
invasion.exp.data.16_zscores$min.10.pH.unscaled <-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# for visualizing histograms of pH to use as continuous vs. discrete variable
invasion.exp.data.16_pres<-invasion.exp.data.16 %>% filter(Invasives=="Present")
invasion.exp.data.16_abs<-invasion.exp.data.16 %>% filter(Invasives=="Absent")


invasion.exp.data.8.community<-read.csv("C:Biological data/invasion.exp.data.8.community.csv")
invasion.exp.data.8<-merge(invasion.exp.data.8, invasion.exp.data.8.community, by="Tile.ID")
# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
invasion.exp.data.8_zscores<-invasion.exp.data.8
#invasion.exp.data.8_zscores$hydrogen.concentration<-scale(invasion.exp.data.8$hydrogen.concentration, center=TRUE, scale=TRUE)
invasion.exp.data.8_zscores$av.pH<-scale(invasion.exp.data.8$av.pH, center=TRUE, scale=TRUE)
invasion.exp.data.8_zscores$min.10.pH<-scale(invasion.exp.data.8$min.10.pH, center=TRUE, scale=TRUE)
invasion.exp.data.8_zscores$Mesocosm <- as.factor(invasion.exp.data.8$Mesocosm)
invasion.exp.data.8_zscores$av.pH.unscaled <-invasion.exp.data.8_zscores$av.pH * attr(invasion.exp.data.8_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$av.pH, 'scaled:center')
invasion.exp.data.8_zscores$min.10.pH.unscaled <-invasion.exp.data.8_zscores$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# Notes on contrasts ------------------------------------------------------

#Notes on ordered factors: the factor Invasives is ordered
#options(contrasts = c("contr.sum", "contr.poly"))
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.


# GAM information and resources---------------------------------------------------------------------
#Gavin Simpson's gratia - used to visualize gams

#We have a varying coefficient model aka ANCOVA, so use "by"
#We use gam(...,select=TRUE), automatic model selection via null space penalization

#Invasives is an ordered factor - meaning "none" is the base level and the other two relate to that....
# see https://www.fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/



# Plotting settings -------------------------------------------------------

colorset_invasives = c("Present"="#0276FD" ,"Absent"="#818392")
theme_set(theme_classic(base_size = 6))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

 
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
##appraise(gam.16.beta.botryllid)
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

plt.inv.botryllid.16 <- ggplot(ndata.16.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.botryllid.16
ggsave("C:Graphs August 2020//botryllid_pred.16.png")


# GAM beta botryllus / gam.8.beta.botryllid.2 ----------------------------------------------------

#binomial first
gam.8.binomial.botryllid<- gam(formula = cbind(botryllid, 100-botryllid)~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.8.beta.botryllid<- gam(botryllid.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.1<- gam(botryllid.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.2<- gam(botryllid.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.8.beta.botryllid.3<- gam(botryllid.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives,k=4), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.8.beta.botryllid, gam.8.beta.botryllid.1, gam.8.beta.botryllid.2, gam.8.beta.botryllid.3, gam.8.binomial.botryllid)
#.3 is best

plot(gam.8.beta.botryllid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.beta.botryllid)
qq_plot(gam.8.beta.botryllid, method = 'simulate')
k.check(gam.8.beta.botryllid)
summary(gam.8.beta.botryllid)

gam.8.beta.botryllid.unordered<- gam(botryllid.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=Invasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)


fam.gam.8.botryllid <- family(gam.8.beta.botryllid)
fam.gam.8.botryllid
ilink.gam.8.botryllid<- fam.gam.8.botryllid$linkinv
ilink.gam.8.botryllid
want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)


mod.botryllid<-gam.8.beta.botryllid
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

plt.inv.botryllid.8 <- ggplot(ndata.8.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.botryllid.8
ggsave("C:Graphs August 2020//botryllid_pred.8.png")



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
##appraise(gam.16.beta.botryllid.3)
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

plt.inv.botryllid <- ggplot(ndata.16.botryllid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = botryllid.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.botryllid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.botryllid
ggsave("C:Graphs August 2020//botryllid_pred.png")



# GAM beta botryllus / gam.16.beta.bot.eaten.2 ----------------------------------------------------

#binomial first
gam.16.binomial.bot.eaten<- gam(formula = cbind(bot.eaten, 100-bot.eaten)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.16.beta.bot.eaten<- gam(bot.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.16.beta.bot.eaten.1<- gam(bot.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.16.beta.bot.eaten.2<- gam(bot.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.16.beta.bot.eaten.3<- gam(bot.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.16.beta.bot.eaten, gam.16.beta.bot.eaten.1, gam.16.beta.bot.eaten.2, gam.16.beta.bot.eaten.3, gam.16.binomial.bot.eaten)


plot(gam.16.beta.bot.eaten, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.bot.eaten)
qq_plot(gam.16.beta.bot.eaten, method = 'simulate')
k.check(gam.16.beta.bot.eaten)
summary(gam.16.beta.bot.eaten)

gam.16.beta.bot.eaten.unordered<- gam(bot.eaten.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)


fam.gam.16.bot.eaten <- family(gam.16.beta.bot.eaten)
fam.gam.16.bot.eaten
ilink.gam.16.bot.eaten<- fam.gam.16.bot.eaten$linkinv
ilink.gam.16.bot.eaten
want <- seq(1, nrow(invasion.exp.data.16_zscores), length.out = 100)


mod.bot.eaten<-gam.16.beta.bot.eaten
ndata.16.bot.eaten <- with(invasion.exp.data.16_zscores, 
                           data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                      length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.bot.eatenel for the new data
ndata.16.bot.eaten <- add_column(ndata.16.bot.eaten, fit = predict(mod.bot.eaten, newdata = ndata.16.bot.eaten, type = 'response'))


ndata.16.bot.eaten <- bind_cols(ndata.16.bot.eaten, setNames(as_tibble(predict(mod.bot.eaten, ndata.16.bot.eaten, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.bot.eaten <- mutate(ndata.16.bot.eaten,
                             fit_resp  = ilink.gam.16.bot.eaten(fit_link),
                             right_upr = ilink.gam.16.bot.eaten(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.16.bot.eaten(fit_link - (2 * se_link)))


invasion.exp.data.16_zscores$min.10.pH.unscaled<-invasion.exp.data.16_zscores$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

ndata.16.bot.eaten$min.10.pH.unscaled<-ndata.16.bot.eaten$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.bot.eaten.16 <- ggplot(ndata.16.bot.eaten, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = bot.eaten.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "proportion eaten")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.bot.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.bot.eaten.16
ggsave("C:Graphs August 2020//bot.eaten_pred.16.png")


# GAM beta botryllus / gam.8.beta.bot.eaten.2 ----------------------------------------------------

#binomial first
gam.8.binomial.bot.eaten<- gam(formula = cbind(bot.eaten, 100-bot.eaten)~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.8.beta.bot.eaten<- gam(bot.eaten.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.8.beta.bot.eaten.1<- gam(bot.eaten.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.8.beta.bot.eaten.2<- gam(bot.eaten.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.8.beta.bot.eaten.3<- gam(bot.eaten.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.8.beta.bot.eaten, gam.8.beta.bot.eaten.1, gam.8.beta.bot.eaten.2, gam.8.beta.bot.eaten.3, gam.8.binomial.bot.eaten)
#.3 is best

plot(gam.8.beta.bot.eaten.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.beta.bot.eaten.3)
qq_plot(gam.8.beta.bot.eaten.3, method = 'simulate')
k.check(gam.8.beta.bot.eaten.3)
summary(gam.8.beta.bot.eaten.3)

gam.8.beta.bot.eaten.3.unordered<- gam(bot.eaten.001~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=Invasives, k=4), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


fam.gam.8.bot.eaten <- family(gam.8.beta.bot.eaten.3)
fam.gam.8.bot.eaten
ilink.gam.8.bot.eaten<- fam.gam.8.bot.eaten$linkinv
ilink.gam.8.bot.eaten
want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)


mod.bot.eaten<-gam.8.beta.bot.eaten.3
ndata.8.bot.eaten <- with(invasion.exp.data.8_zscores, 
                          data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))

## add the fitted values by predicting from the mod.bot.eatenel for the new data
ndata.8.bot.eaten <- add_column(ndata.8.bot.eaten, fit = predict(mod.bot.eaten, newdata = ndata.8.bot.eaten, type = 'response'))


ndata.8.bot.eaten <- bind_cols(ndata.8.bot.eaten, setNames(as_tibble(predict(mod.bot.eaten, ndata.8.bot.eaten, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.bot.eaten <- mutate(ndata.8.bot.eaten,
                            fit_resp  = ilink.gam.8.bot.eaten(fit_link),
                            right_upr = ilink.gam.8.bot.eaten(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.8.bot.eaten(fit_link - (2 * se_link)))


invasion.exp.data.8_zscores$min.10.pH.unscaled<-invasion.exp.data.8_zscores$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

ndata.8.bot.eaten$min.10.pH.unscaled<-ndata.8.bot.eaten$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.bot.eaten.8 <- ggplot(ndata.8.bot.eaten, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = bot.eaten.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "proportion eaten")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.bot.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.bot.eaten.8
ggsave("C:Graphs August 2020//bot.eaten_pred.8.png")


# GAM beta folliculina / gam.16.beta.folliculina -----------------------------------------------------------

gam.16.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")
gam.16.beta.folliculina<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.folliculina.1<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.folliculina.2<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.folliculina.3<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.16.beta.folliculina, gam.16.beta.folliculina.1, gam.16.beta.folliculina.2,gam.16.binomial.folliculina, gam.16.beta.folliculina.3)
#simplest logit


plot(gam.16.beta.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.folliculina)
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

plt.inv.folliculina.16 <- ggplot(ndata.16.folliculina, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = folliculina.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.folliculina.16
ggsave("C:Graphs av.pH//folliculina_pred.16.png")

# GAM beta folliculina / gam.8.beta.folliculina -----------------------------------------------------------

gam.8.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.folliculina<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.folliculina.1<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.folliculina.2<- gam(folliculina.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.folliculina.3<- gam(folliculina.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.8.beta.folliculina, gam.8.beta.folliculina.1, gam.8.beta.folliculina.2,gam.8.binomial.folliculina, gam.8.beta.folliculina.3)
#simplest logit


plot(gam.8.beta.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
##appraise(gam.8.beta.folliculina)
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

plt.inv.folliculina.8 <- ggplot(ndata.8.folliculina, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = folliculina.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.folliculina.8
ggsave("C:Graphs August 2020//folliculina_pred.8.png")


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
#appraise(gam.16.beta.membranipora.3)
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

plt.inv.membranipora.16 <- ggplot(ndata.16.membranipora, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = membranipora.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.membranipora.16
ggsave("C:Graphs August 2020//membranipora_pred.16.png")

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
##appraise(gam.8.beta.membranipora.3)
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

plt.inv.membranipora.8 <- ggplot(ndata.8.membranipora, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = membranipora.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.membranipora.8
ggsave("C:Graphs August 2020//membranipora_pred.8.png")


# GAM beta mem.eaten / gam.16.beta.mem.eaten --------------------------------------------------------

#binomial first
gam.16.binomial.mem.eaten<- gam(formula = cbind(mem.eaten, 100-mem.eaten)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.16.beta.mem.eaten<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.mem.eaten.1<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.mem.eaten.2<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.mem.eaten.3<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.16.beta.mem.eaten, gam.16.beta.mem.eaten.1, gam.16.beta.mem.eaten.2, gam.16.binomial.mem.eaten, gam.16.beta.mem.eaten.3)
#cauchit is best


plot(gam.16.beta.mem.eaten.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.mem.eaten.3)
#not great
qq_plot(gam.16.beta.mem.eaten.3, method = 'simulate')
k.check(gam.16.beta.mem.eaten.3)
summary(gam.16.beta.mem.eaten.3)
vis.gam(gam.16.beta.mem.eaten.3)

gam.16.beta.mem.eaten.3.unordered<- gam(mem.eaten.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


fam.gam.16.mem.eaten <- family(gam.16.beta.mem.eaten.3)
fam.gam.16.mem.eaten
ilink.gam.16.mem.eaten<- fam.gam.16.mem.eaten$linkinv
ilink.gam.16.mem.eaten


mod.mem.eaten<-gam.16.beta.mem.eaten.3
ndata.16.mem.eaten <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.mem.eatenel for the new data
ndata.16.mem.eaten <- add_column(ndata.16.mem.eaten, fit = predict(mod.mem.eaten, newdata = ndata.16.mem.eaten, type = 'response'))


ndata.16.mem.eaten <- bind_cols(ndata.16.mem.eaten, setNames(as_tibble(predict(mod.mem.eaten, ndata.16.mem.eaten, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.mem.eaten <- mutate(ndata.16.mem.eaten,
                                fit_resp  = ilink.gam.16.mem.eaten(fit_link),
                                right_upr = ilink.gam.16.mem.eaten(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.16.mem.eaten(fit_link - (2 * se_link)))

ndata.16.mem.eaten$min.10.pH.unscaled<-ndata.16.mem.eaten$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.mem.eaten.16 <- ggplot(ndata.16.mem.eaten, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = mem.eaten.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("mem.eaten")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.mem.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.mem.eaten.16
ggsave("C:Graphs August 2020//mem.eaten_pred.16.png")

# GAM beta mem.eaten / gam.8.beta.mem.eaten --------------------------------------------------------

#binomial first
gam.8.binomial.mem.eaten<- gam(formula = cbind(mem.eaten, 100-mem.eaten)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.8.beta.mem.eaten<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.mem.eaten.1<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.mem.eaten.2<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.mem.eaten.3<- gam(mem.eaten.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.8.beta.mem.eaten, gam.8.beta.mem.eaten.1, gam.8.beta.mem.eaten.2, gam.8.binomial.mem.eaten, gam.8.beta.mem.eaten.3)
#cauchit is best


plot(gam.8.beta.mem.eaten.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.beta.mem.eaten.3)
qq_plot(gam.8.beta.mem.eaten.3, method = 'simulate')
k.check(gam.8.beta.mem.eaten.3)
summary(gam.8.beta.mem.eaten.3)
vis.gam(gam.8.beta.mem.eaten.3)

gam.8.beta.mem.eaten.3.unordered<- gam(mem.eaten.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


fam.gam.8.mem.eaten <- family(gam.8.beta.mem.eaten.3)
fam.gam.8.mem.eaten
ilink.gam.8.mem.eaten<- fam.gam.8.mem.eaten$linkinv
ilink.gam.8.mem.eaten


mod.mem.eaten<-gam.8.beta.mem.eaten.3
ndata.8.mem.eaten <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                     length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.mem.eatenel for the new data
ndata.8.mem.eaten <- add_column(ndata.8.mem.eaten, fit = predict(mod.mem.eaten, newdata = ndata.8.mem.eaten, type = 'response'))


ndata.8.mem.eaten <- bind_cols(ndata.8.mem.eaten, setNames(as_tibble(predict(mod.mem.eaten, ndata.8.mem.eaten, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.mem.eaten <- mutate(ndata.8.mem.eaten,
                               fit_resp  = ilink.gam.8.mem.eaten(fit_link),
                               right_upr = ilink.gam.8.mem.eaten(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.mem.eaten(fit_link - (2 * se_link)))

ndata.8.mem.eaten$min.10.pH.unscaled<-ndata.8.mem.eaten$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.mem.eaten.8 <- ggplot(ndata.8.mem.eaten, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = mem.eaten.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("mem.eaten")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.mem.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.mem.eaten.8
ggsave("C:Graphs August 2020//mem.eaten_pred.8.png")

 
# GAM beta mussel / gam.16.beta.mussel --------------------------------------------------------

#binomial first
gam.16.binomial.mussel<- gam(formula = cbind(mussel, 100-mussel)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.16.beta.mussel<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.mussel.1<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.mussel.2<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.mussel.3<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.16.beta.mussel, gam.16.beta.mussel.1, gam.16.beta.mussel.2, gam.16.binomial.mussel, gam.16.beta.mussel.3)



plot(gam.16.beta.mussel, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
##appraise(gam.16.beta.mussel)
qq_plot(gam.16.beta.mussel, method = 'simulate')
k.check(gam.16.beta.mussel)
summary(gam.16.beta.mussel)
vis.gam(gam.16.beta.mussel)

gam.16.beta.mussel.unordered<- gam(mussel.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logitt"), select=TRUE, method="REML")


fam.gam.16.mussel <- family(gam.16.beta.mussel.3)
fam.gam.16.mussel
ilink.gam.16.mussel<- fam.gam.16.mussel$linkinv
ilink.gam.16.mussel


mod.mussel<-gam.16.beta.mussel.3
ndata.16.mussel <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.musselel for the new data
ndata.16.mussel <- add_column(ndata.16.mussel, fit = predict(mod.mussel, newdata = ndata.16.mussel, type = 'response'))


ndata.16.mussel <- bind_cols(ndata.16.mussel, setNames(as_tibble(predict(mod.mussel, ndata.16.mussel, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.mussel <- mutate(ndata.16.mussel,
                                fit_resp  = ilink.gam.16.mussel(fit_link),
                                right_upr = ilink.gam.16.mussel(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.16.mussel(fit_link - (2 * se_link)))

ndata.16.mussel$min.10.pH.unscaled<-ndata.16.mussel$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.mussel.16 <- ggplot(ndata.16.mussel, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = mussel.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("mussel")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.mussel,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.mussel.16
ggsave("C:Graphs August 2020//mussel_pred.16.png")

# GAM beta mussel / gam.8.beta.mussel --------------------------------------------------------

#binomial first
gam.8.binomial.mussel<- gam(formula = cbind(mussel, 100-mussel)~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.8.beta.mussel<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.mussel.1<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.mussel.2<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.mussel.3<- gam(mussel.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.8.beta.mussel, gam.8.beta.mussel.1, gam.8.beta.mussel.2, gam.8.binomial.mussel, gam.8.beta.mussel.3)



plot(gam.8.beta.mussel, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.beta.mussel)
qq_plot(gam.8.beta.mussel, method = 'simulate')
k.check(gam.8.beta.mussel)
summary(gam.8.beta.mussel)
vis.gam(gam.8.beta.mussel)

gam.8.beta.mussel.unordered<- gam(mussel.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.8.mussel <- family(gam.8.beta.mussel)
fam.gam.8.mussel
ilink.gam.8.mussel<- fam.gam.8.mussel$linkinv
ilink.gam.8.mussel


mod.mussel<-gam.8.beta.mussel
ndata.8.mussel <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the mod.musselel for the new data
ndata.8.mussel <- add_column(ndata.8.mussel, fit = predict(mod.mussel, newdata = ndata.8.mussel, type = 'response'))


ndata.8.mussel <- bind_cols(ndata.8.mussel, setNames(as_tibble(predict(mod.mussel, ndata.8.mussel, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.mussel <- mutate(ndata.8.mussel,
                          fit_resp  = ilink.gam.8.mussel(fit_link),
                          right_upr = ilink.gam.8.mussel(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.8.mussel(fit_link - (2 * se_link)))

ndata.8.mussel$min.10.pH.unscaled<-ndata.8.mussel$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')

# plot 

plt.inv.mussel.8 <- ggplot(ndata.8.mussel, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = mussel.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("mussel")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.mussel,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.mussel.8
ggsave("C:Graphs August 2020//mussel_pred.8.png")


# GAM negbin barnacles / gam.16.nb.num.barn -----------------------------------------------------------
nbinom.16.barn <- fitdistr(invasion.exp.data.16_zscores$num.barn, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.barn, "nbinom", size = nbinom.16.barn$estimate[[1]], mu = nbinom.16.barn$estimate[[2]])

#negative binomial 
gam.16.nb.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.barn$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.barn.1<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.barn<- gam(num.barn ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.16.nb.num.barn, gam.16.nb.num.barn.1, gam.16.poisson.num.barn)

plot(gam.16.poisson.num.barn, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
##appraise(gam.16.poisson.num.barn)
qq_plot(gam.16.poisson.num.barn, method = 'simulate')
#looks really good
k.check(gam.16.poisson.num.barn)
summary(gam.16.poisson.num.barn)


#a few outside the area
###appraise a bit funnelly

gam.16.poisson.num.barn.unordered<- gam(num.barn ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson(), select=TRUE, method="REML")
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
plt.inv.num.barn.16 <- ggplot(ndata.16.num.barn, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.barn, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.barn.16
ggsave("C:Graphs August 2020//num.barn_pred.16.png")

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
##appraise(gam.8.poisson.num.barn)
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
plt.inv.num.barn.8 <- ggplot(ndata.8.num.barn, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.barn, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.barn.8
ggsave("C:Graphs August 2020//num.barn_pred.8.png")


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
##appraise(gam.16.nb.num.white.bryo)
qq_plot(gam.16.nb.num.white.bryo, method = 'simulate')
k.check(gam.16.nb.num.white.bryo)
summary(gam.16.nb.num.white.bryo)

#a few outside the area
###appraise a bit funnelly

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
plt.inv.num.white.bryo.16 <- ggplot(ndata.16.num.white.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.white.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.white.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.white.bryo.16
ggsave("C:Graphs August 2020//num.white.bryo_pred.16.png")


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
##appraise(gam.8.poisson.num.white.bryo)
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
plt.inv.num.white.bryo.8 <- ggplot(ndata.8.num.white.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.white.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.white.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.white.bryo.8
ggsave("C:Graphs August 2020//num.white.bryo_pred.8.png")


# GAM negbin num.red.bryo / gam.16.nb.num.red.bryo --------------------------------------------------------------
nbinom.16.num.red.bryo <- fitdistr(invasion.exp.data.16_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.red.bryo, "nbinom", size = nbinom.16.num.red.bryo$estimate[[1]], mu = nbinom.16.num.red.bryo$estimate[[2]])
#getting theta

gam.16.nb.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.red.bryo, gam.16.nb.num.red.bryo.1, gam.16.poisson.num.red.bryo)


plot(gam.16.poisson.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.poisson.num.red.bryo)
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
plt.inv.num.red.bryo.16 <- ggplot(ndata.16.num.red.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.red.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.red.bryo.16
ggsave("C:Graphs August 2020//num.red.bryo_pred.16.png")


# GAM negbin num.red.bryo / gam.8.nb.num.red.bryo --------------------------------------------------------------
nbinom.8.num.red.bryo <- fitdistr(invasion.exp.data.8_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.red.bryo, "nbinom", size = nbinom.8.num.red.bryo$estimate[[1]], mu = nbinom.8.num.red.bryo$estimate[[2]])
#getting theta
invasion.exp.data.8_zscores$num.red.bryo


gam.8.nb.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.red.bryo<- gam(num.red.bryo ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.red.bryo, gam.8.nb.num.red.bryo.1, gam.8.poisson.num.red.bryo)


plot(gam.8.poisson.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
##appraise(gam.8.poisson.num.red.bryo)
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
plt.inv.num.red.bryo.8 <- ggplot(ndata.8.num.red.bryo, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.red.bryo, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,2.5)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.red.bryo.8
ggsave("C:Graphs August 2020//num.red.bryo_pred.8.png")


# GAM poisson num nudi / gam.16.poisson.num.nudi  ------------------------------------------------------------
nbinom.16.num.nudi <- fitdistr(invasion.exp.data.16_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.nudi, "nbinom", size = nbinom.16.num.nudi$estimate[[1]], mu = nbinom.16.num.nudi$estimate[[2]])
#theta

gam.16.nb.num.nudi.1<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.nb.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.16.poisson.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
AICtab(gam.16.nb.num.nudi, gam.16.nb.num.nudi.1, gam.16.poisson.num.nudi)

#appraise(gam.16.poisson.num.nudi)
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
plt.inv.num.nudi.16 <- ggplot(ndata.16.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.nudi, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.inv.num.nudi.16
ggsave("C:Graphs August 2020//num.nudi_pred.16.png")


# GAM poisson num nudi / gam.8.poisson.num.nudi  ------------------------------------------------------------
nbinom.8.num.nudi <- fitdistr(invasion.exp.data.8_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.nudi, "nbinom", size = nbinom.8.num.nudi$estimate[[1]], mu = nbinom.8.num.nudi$estimate[[2]])
#theta

gam.8.nb.num.nudi.1<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.nb.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.16.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.8.poisson.num.nudi<- gam(num.nudi ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")
AICtab(gam.8.nb.num.nudi, gam.8.nb.num.nudi.1)


#appraise(gam.8.nb.num.nudi.1)
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
plt.inv.num.nudi.8 <- ggplot(ndata.8.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.nudi, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,4)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.inv.num.nudi.8
ggsave("C:Graphs August 2020//num.nudi_pred.8.png")

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

##appraise(gam.16.nb.num.serpulid.1)
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
plt.inv.num.serpulid.16 <- ggplot(ndata.16.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.serpulid.16
ggsave("C:Graphs August 2020//num.serpulid_pred.16.png")

# GAM nb() serpulids / gam.8.nb.num.serpulid.1 -----------------------------------------------------------
nbinom.8.num.serpulid <- fitdistr(invasion.exp.data.8_zscores$num.serpulid, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.serpulid, "nbinom", size = nbinom.8.num.serpulid$estimate[[1]], mu = nbinom.8.num.serpulid$estimate[[2]])
#theta

#negative binomial first
gam.8.nb.num.serpulid<- gam(num.serpulid ~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.serpulid.1<- gam(num.serpulid ~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.serpulid<- gam(num.serpulid ~ s(min.10.pH, k=4)+ oInvasives + s(min.10.pH, by=oInvasives, k=4),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.serpulid, gam.8.nb.num.serpulid.1, gam.8.poisson.num.serpulid)
##gam.8.poisson is best

##appraise(gam.8.nb.num.serpulid.1)
qq_plot(gam.8.nb.num.serpulid.1, method = 'simulate')
plot(gam.8.nb.num.serpulid.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.8.nb.num.serpulid.1)
summary(gam.8.nb.num.serpulid.1)

gam.8.nb.num.serpulid.1.unordered<- gam(num.serpulid ~ s(min.10.pH, k=4)+ Invasives + s(min.10.pH, by=oInvasives, k=4),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")


want <- seq(1, nrow(invasion.exp.data.8_zscores), length.out = 100)
fam.gam.8.num.serpulid <- family(gam.8.nb.num.serpulid.1)
fam.gam.8.num.serpulid
ilink.gam.8.num.serpulid<- fam.gam.8.num.serpulid$linkinv
ilink.gam.8.num.serpulid
mod.num.serpulid<-gam.8.nb.num.serpulid.1
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
plt.inv.num.serpulid.8 <- ggplot(ndata.8.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.serpulid, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+ylim(0,6)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.serpulid.8
ggsave("C:Graphs August 2020//num.serpulid_pred.8.png")



# GAM negbin corella / gam.16.nb.num.corella -------------------------------------------------------------

nbinom.16.num.corella <- fitdistr(invasion.exp.data.16_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data.16_zscores$num.corella, "nbinom", size = nbinom.16.num.corella$estimate[[1]], mu = nbinom.16.num.corella$estimate[[2]])
#extracting theta

gam.16.nb.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = negbin(nbinom.16.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.16.nb.num.corella.1<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.corella, gam.16.nb.num.corella.1, gam.16.poisson.num.corella)


#appraise(gam.16.nb.num.corella)
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
plt.inv.num.corella.16 <- ggplot(ndata.16.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.corella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.corella.16
ggsave("C:Graphs August 2020//num.corella_pred.16.png")

# GAM negbin corella / gam.8.nb.num.corella -------------------------------------------------------------

nbinom.8.num.corella <- fitdistr(invasion.exp.data.8_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.corella, "nbinom", size = nbinom.8.num.corella$estimate[[1]], mu = nbinom.8.num.corella$estimate[[2]])
#extracting theta

gam.8.nb.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.corella.1<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.corella<- gam(num.corella ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.corella, gam.8.nb.num.corella.1, gam.8.poisson.num.corella)
#poisson

##appraise(gam.8.poisson.num.corella)
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
plt.inv.num.corella.8 <- ggplot(ndata.8.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.corella, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.num.corella.8
ggsave("C:Graphs August 2020//num.corella_pred.8.png")



# Fig 2 plot generation ---------------------------------------------------
fig.week.16<-wrap_plots(plt.inv.botryllid.16,plt.inv.folliculina.16,plt.inv.membranipora.16,
          plt.inv.mussel.16,plt.inv.num.barn.16,plt.inv.num.white.bryo.16,plt.inv.num.red.bryo.16,plt.inv.num.nudi.16,plt.inv.num.serpulid.16,
          plt.inv.num.corella.16, ncol=5)+
          plot_annotation(tag_levels = 'A')

theme_set(theme_classic(base_size = 5))
#theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

fig.week.16
ggplot2::ggsave(plot=fig.week.16, "C:Graphs August 2020//Fig_Invasion_pres_wk_16.pdf", width=18, height=8, units="cm")



fig.week.8<-wrap_plots(plt.inv.botryllid.8,plt.inv.folliculina.8,plt.inv.membranipora.8,
                        plt.inv.mussel.8,plt.inv.num.barn.8,plt.inv.num.white.bryo.8,plt.inv.num.red.bryo.8,plt.inv.num.nudi.8,plt.inv.num.serpulid.8,
                        plt.inv.num.corella.8, ncol=5)+
                        plot_annotation(tag_levels = 'A')

fig.week.8

ggplot2::ggsave(plot=fig.week.8, "C:Graphs August 2020//Fig_Invasion_pres_wk_8.pdf", width=18, height=8, units="cm")

head(invasion.exp.data.8_zscores)

# Pulling model 16 results to a table ----------------------------------------

botryllus.gam.16<-summary(gam.16.beta.botryllid)
botryllus.eaten.gam.16<-summary(gam.16.beta.bot.eaten)
folliculina.gam.16<-summary(gam.16.beta.folliculina)
membranipora.gam.16<-summary(gam.16.beta.membranipora.3)
membranipora.eaten.gam.16<-summary(gam.16.beta.mem.eaten.3)
mussel.gam.16<-summary(gam.16.beta.mussel.3)
alive.barn.gam.16<-summary(gam.16.poisson.num.barn)
num.white.bryo.gam.16<-summary(gam.16.nb.num.white.bryo)
num.red.bryo.gam.16<-summary(gam.16.poisson.num.red.bryo)
num.nudi.gam.16<-summary(gam.16.poisson.num.nudi)
num.serpulid.gam.16<-summary(gam.16.nb.num.serpulid.1)
corella.gam.16<-summary(gam.16.nb.num.corella)


botryllus.gam.16.unordered<-summary(gam.16.beta.botryllid.unordered)
botryllus.eaten.gam.16.unordered<-summary(gam.16.beta.bot.eaten.unordered)
folliculina.gam.16.unordered<-summary(gam.16.beta.folliculina.unordered)
membranipora.gam.16.unordered<-summary(gam.16.beta.membranipora.3.unordered)
membranipora.eaten.gam.16.unordered<-summary(gam.16.beta.mem.eaten.3.unordered)
mussel.gam.16.unordered<-summary(gam.16.beta.mussel.3.unordered)
alive.barn.gam.16.unordered<-summary(gam.16.poisson.num.barn.unordered)
num.white.bryo.gam.16.unordered<-summary(gam.16.nb.num.white.bryo.unordered)
num.red.bryo.gam.16.unordered<-summary(gam.16.poisson.num.red.bryo.unordered)
num.nudi.gam.16.unordered<-summary(gam.16.poisson.num.nudi.unordered)
num.serpulid.gam.16.unordered<-summary(gam.16.nb.num.serpulid.1.unordered)
corella.gam.16.unordered<-summary(gam.16.nb.num.corella.unordered)


botryllus.gam.16.p.table<-as.data.frame(botryllus.gam.16.unordered$p.table)
botryllus.gam.16.s.table<-as.data.frame(botryllus.gam.16$s.table)

botryllus.eaten.gam.16.p.table<-as.data.frame(botryllus.eaten.gam.16.unordered$p.table)
botryllus.eaten.gam.16.s.table<-as.data.frame(botryllus.eaten.gam.16$s.table)

folliculina.gam.16.p.table<-as.data.frame(folliculina.gam.16.unordered$p.table)
folliculina.gam.16.s.table<-as.data.frame(folliculina.gam.16$s.table)

membranipora.gam.16.p.table<-as.data.frame(membranipora.gam.16.unordered$p.table)
membranipora.gam.16.s.table<-as.data.frame(membranipora.gam.16$s.table)

membranipora.eaten.gam.16.p.table<-as.data.frame(membranipora.eaten.gam.16.unordered$p.table)
membranipora.eaten.gam.16.s.table<-as.data.frame(membranipora.eaten.gam.16$s.table)

mussel.gam.16.p.table<-as.data.frame(mussel.gam.16.unordered$p.table)
mussel.gam.16.s.table<-as.data.frame(mussel.gam.16$s.table)

alive.barn.gam.16.p.table<-as.data.frame(alive.barn.gam.16.unordered$p.table)
alive.barn.gam.16.s.table<-as.data.frame(alive.barn.gam.16$s.table)

num.white.bryo.gam.16.p.table<-as.data.frame(num.white.bryo.gam.16.unordered$p.table)
num.white.bryo.gam.16.s.table<-as.data.frame(num.white.bryo.gam.16$s.table)

num.red.bryo.gam.16.p.table<-as.data.frame(num.red.bryo.gam.16.unordered$p.table)
num.red.bryo.gam.16.s.table<-as.data.frame(num.red.bryo.gam.16$s.table)

num.nudi.gam.16.p.table<-as.data.frame(num.nudi.gam.16.unordered$p.table)
num.nudi.gam.16.s.table<-as.data.frame(num.nudi.gam.16$s.table)

num.serpulid.gam.16.p.table<-as.data.frame(num.serpulid.gam.16.unordered$p.table)
num.serpulid.gam.16.s.table<-as.data.frame(num.serpulid.gam.16$s.table)

corella.gam.16.p.table<-as.data.frame(corella.gam.16.unordered$p.table)
corella.gam.16.s.table<-as.data.frame(corella.gam.16$s.table)

#### Building the stats table
ptable.16<-rbind( botryllus.gam.16.p.table, 
               botryllus.eaten.gam.16.p.table,
               folliculina.gam.16.p.table,
               membranipora.gam.16.p.table, 
               membranipora.eaten.gam.16.p.table,
               mussel.gam.16.p.table,
               alive.barn.gam.16.p.table,
               num.white.bryo.gam.16.p.table,
               num.red.bryo.gam.16.p.table,
               num.nudi.gam.16.p.table,
               num.serpulid.gam.16.p.table,
               corella.gam.16.p.table)


colnames(ptable.16) <- c("Estimate", "SE", "z", "p")
ptable.16$Factor<-rep(c("Intercept", "Invasives Present"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.16 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 18) %>% 
  group_rows("Hermissenda", 19, 20) %>% 
  group_rows("Serpulid", 21, 22) %>% 
  group_rows("Corella", 23, 24) %>% 
save_kable(file = "C:Biological data//ptable.16.html", self_contained = T)


### s table
stable.16<-rbind(botryllus.gam.16.s.table, 
              botryllus.eaten.gam.16.s.table,
              folliculina.gam.16.s.table,
              membranipora.gam.16.s.table,
              membranipora.eaten.gam.16.s.table,
              mussel.gam.16.s.table,
              alive.barn.gam.16.s.table,
              num.white.bryo.gam.16.s.table,
              num.red.bryo.gam.16.s.table,
              num.nudi.gam.16.s.table,
              num.serpulid.gam.16.s.table,
              corella.gam.16.s.table)


colnames(stable.16) <- c("Estimated_df", "Reference_df", "Chi_squared", "p_smooth")
stable.16$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Invasives present"))

stable.16 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 18) %>% 
  group_rows("Hermissenda", 19, 20) %>% 
  group_rows("Serpulid", 21, 22) %>% 
  group_rows("Corella", 23, 24) %>% 
  save_kable(file = "C:Biological data//stable.16.html", self_contained = T)

  
pstable.16<-cbind(ptable.16, stable.16)

pstable.16 %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth, Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 18) %>% 
  group_rows("Hermissenda", 19, 20) %>% 
  group_rows("Serpulid", 21, 22) %>% 
  group_rows("Corella", 23, 24) %>% 
  save_kable(file = "C:Biological data//pstable.16.html", self_contained = T)



# Pulling model 8  results to a table ----------------------------------------

botryllus.gam.8<-summary(gam.8.beta.botryllid.3)
botryllus.eaten.gam.8<-summary(gam.8.beta.bot.eaten.3)
folliculina.gam.8<-summary(gam.8.beta.folliculina)
membranipora.gam.8<-summary(gam.8.beta.membranipora.3)
membranipora.eaten.gam.8<-summary(gam.8.beta.mem.eaten.3)
mussel.gam.8<-summary(gam.8.beta.mussel.3)
alive.barn.gam.8<-summary(gam.8.poisson.num.barn)
num.white.bryo.gam.8<-summary(gam.8.poisson.num.white.bryo)
num.nudi.gam.8<-summary(gam.8.nb.num.nudi.1)
num.serpulid.gam.8<-summary(gam.8.poisson.num.serpulid)
corella.gam.8<-summary(gam.8.poisson.num.corella)

botryllus.gam.8.unordered<-summary(gam.8.beta.botryllid.3.unordered)
botryllus.eaten.gam.8.unordered<-summary(gam.8.beta.bot.eaten.3.unordered)
folliculina.gam.8.unordered<-summary(gam.8.beta.folliculina.unordered)
membranipora.gam.8.unordered<-summary(gam.8.beta.membranipora.3.unordered)
membranipora.eaten.gam.8.unordered<-summary(gam.8.beta.mem.eaten.3.unordered)
mussel.gam.8.unordered<-summary(gam.8.beta.mussel.3.unordered)
alive.barn.gam.8.unordered<-summary(gam.8.poisson.num.barn.unordered)
num.white.bryo.gam.8.unordered<-summary(gam.8.poisson.num.white.bryo.unordered)
num.nudi.gam.8.unordered<-summary(gam.8.nb.num.nudi.1.unordered)
num.serpulid.gam.8.unordered<-summary(gam.8.poisson.num.serpulid.unordered)
corella.gam.8.unordered<-summary(gam.8.poisson.num.corella.unordered)

botryllus.gam.8.p.table<-as.data.frame(botryllus.gam.8.unordered$p.table)
botryllus.gam.8.s.table<-as.data.frame(botryllus.gam.8$s.table)

botryllus.eaten.gam.8.p.table<-as.data.frame(botryllus.eaten.gam.8.unordered$p.table)
botryllus.eaten.gam.8.s.table<-as.data.frame(botryllus.eaten.gam.8$s.table)

folliculina.gam.8.p.table<-as.data.frame(folliculina.gam.8.unordered$p.table)
folliculina.gam.8.s.table<-as.data.frame(folliculina.gam.8$s.table)

membranipora.gam.8.p.table<-as.data.frame(membranipora.gam.8.unordered$p.table)
membranipora.gam.8.s.table<-as.data.frame(membranipora.gam.8$s.table)

membranipora.eaten.gam.8.p.table<-as.data.frame(membranipora.eaten.gam.8.unordered$p.table)
membranipora.eaten.gam.8.s.table<-as.data.frame(membranipora.eaten.gam.8$s.table)

mussel.gam.8.p.table<-as.data.frame(mussel.gam.8.unordered$p.table)
mussel.gam.8.s.table<-as.data.frame(mussel.gam.8$s.table)

alive.barn.gam.8.p.table<-as.data.frame(alive.barn.gam.8.unordered$p.table)
alive.barn.gam.8.s.table<-as.data.frame(alive.barn.gam.8$s.table)

num.white.bryo.gam.8.p.table<-as.data.frame(num.white.bryo.gam.8.unordered$p.table)
num.white.bryo.gam.8.s.table<-as.data.frame(num.white.bryo.gam.8$s.table)

num.nudi.gam.8.p.table<-as.data.frame(num.nudi.gam.8.unordered$p.table)
num.nudi.gam.8.s.table<-as.data.frame(num.nudi.gam.8$s.table)

num.serpulid.gam.8.p.table<-as.data.frame(num.serpulid.gam.8.unordered$p.table)
num.serpulid.gam.8.s.table<-as.data.frame(num.serpulid.gam.8$s.table)

corella.gam.8.p.table<-as.data.frame(corella.gam.8.unordered$p.table)
corella.gam.8.s.table<-as.data.frame(corella.gam.8$s.table)


#### Building the stats table
ptable.8<-rbind( botryllus.gam.8.p.table, 
                  botryllus.eaten.gam.8.p.table,
                  folliculina.gam.8.p.table,
                  membranipora.gam.8.p.table, 
                  membranipora.eaten.gam.8.p.table,
                  mussel.gam.8.p.table,
                  alive.barn.gam.8.p.table,
                  num.white.bryo.gam.8.p.table,
                  num.nudi.gam.8.p.table,
                  num.serpulid.gam.8.p.table,
                  corella.gam.8.p.table)


colnames(ptable.8) <- c("Estimate", "SE", "z", "p")
ptable.8$Factor<-rep(c("Intercept", "Invasives Present"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.8 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 18) %>% 
  group_rows("Serpulid", 19, 20) %>% 
  group_rows("Corella", 21, 22) %>% 
  save_kable(file = "C:Biological data//ptable.8.html", self_contained = T)


### s table
stable.8<-rbind(botryllus.gam.8.s.table, 
                 botryllus.eaten.gam.8.s.table,
                 folliculina.gam.8.s.table,
                 membranipora.gam.8.s.table,
                 membranipora.eaten.gam.8.s.table,
                 mussel.gam.8.s.table,
                 alive.barn.gam.8.s.table,
                 num.white.bryo.gam.8.s.table,
                 num.nudi.gam.8.s.table,
                 num.serpulid.gam.8.s.table,
                 corella.gam.8.s.table)


colnames(stable.8) <- c("Estimated_df", "Reference_df", "Chi_squared", "p_smooth")
stable.8$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Invasives present"))

stable.8 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 18) %>% 
  group_rows("Serpulid", 19, 20) %>% 
  group_rows("Corella", 21, 22) %>% 
  save_kable(file = "C:Biological data//stable.8.html", self_contained = T)


pstable.8<-cbind(ptable.8, stable.8)

pstable.8 %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth, Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 6) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 18) %>% 
  group_rows("Serpulid", 19, 20) %>% 
  group_rows("Corella", 21, 22) %>% 
  save_kable(file = "C:Biological data//pstable.8.html", self_contained = T)



# num.species.no.bot 16 ----------------------------------------------------------------
gam.16.nb.num.species.no.bot.1<- gam(num.species.no.bot ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.species.no.bot<- gam(num.species.no.bot ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.species.no.bot.1, gam.16.poisson.num.species.no.bot)

plot(gam.16.poisson.num.species.no.bot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.poisson.num.species.no.bot)
#okay but qq plot not the best on ends
qq_plot(gam.16.poisson.num.species.no.bot, method = 'simulate')
k.check(gam.16.poisson.num.species.no.bot)
summary(gam.16.poisson.num.species.no.bot)
#a few outside the QQ plot on both ends
#low p value for k - but NS and edf is not super close to k-index

gam.16.poisson.num.species.no.bot.unordered<- gam(num.species.no.bot ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")
fam.gam.16.num.species.no.bot <- family(gam.16.poisson.num.species.no.bot)
ilink.gam.16.num.species.no.bot<- fam.gam.16.num.species.no.bot$linkinv
ilink.gam.16.num.species.no.bot


mod.num.species.no.bot<-gam.16.poisson.num.species.no.bot
ndata.16.num.species.no.bot <- with(invasion.exp.data.16_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.num.species.no.bot <- add_column(ndata.16.num.species.no.bot, fit = predict(mod.num.species.no.bot, newdata = ndata.16.num.species.no.bot, type = 'response'))

predict(mod.num.species.no.bot, newdata = ndata.16.num.species.no.bot, type = 'response')
ndata.16.num.species.no.bot <- bind_cols(ndata.16.num.species.no.bot, setNames(as_tibble(predict(mod.num.species.no.bot, ndata.16.num.species.no.bot, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.num.species.no.bot <- mutate(ndata.16.num.species.no.bot,
                               fit_resp  = ilink.gam.16.num.species.no.bot(fit_link),
                               right_upr = ilink.gam.16.num.species.no.bot(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.16.num.species.no.bot(fit_link - (2 * se_link)))


ndata.16.num.species.no.bot$min.10.pH.unscaled<-ndata.16.num.species.no.bot$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.num.species.no.bot.16 <- ggplot(ndata.16.num.species.no.bot, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.species.no.bot, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Native species richness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.species.no.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,20)
plt.inv.num.species.no.bot.16
ggsave("C:Graphs August 2020//native_richness.16.png")


# num.species.no.bot 8 ----------------------------------------------------------------
gam.8.nb.num.species.no.bot.1<- gam(num.species.no.bot ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.species.no.bot<- gam(num.species.no.bot ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.species.no.bot.1, gam.8.poisson.num.species.no.bot)

plot(gam.8.poisson.num.species.no.bot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.poisson.num.species.no.bot)
#okay but qq plot not the best on ends
qq_plot(gam.8.poisson.num.species.no.bot, method = 'simulate')
k.check(gam.8.poisson.num.species.no.bot)
summary(gam.8.poisson.num.species.no.bot)
#a few outside the QQ plot on both ends
#low p value for k - but NS and edf is not super close to k-index

gam.8.poisson.num.species.no.bot.unordered<- gam(num.species.no.bot ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")
fam.gam.8.num.species.no.bot <- family(gam.8.poisson.num.species.no.bot)
ilink.gam.8.num.species.no.bot<- fam.gam.8.num.species.no.bot$linkinv
ilink.gam.8.num.species.no.bot


mod.num.species.no.bot<-gam.8.poisson.num.species.no.bot
ndata.8.num.species.no.bot <- with(invasion.exp.data.8_zscores, 
                                    data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                               length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.8.num.species.no.bot <- add_column(ndata.8.num.species.no.bot, fit = predict(mod.num.species.no.bot, newdata = ndata.8.num.species.no.bot, type = 'response'))

predict(mod.num.species.no.bot, newdata = ndata.8.num.species.no.bot, type = 'response')
ndata.8.num.species.no.bot <- bind_cols(ndata.8.num.species.no.bot, setNames(as_tibble(predict(mod.num.species.no.bot, ndata.8.num.species.no.bot, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.num.species.no.bot <- mutate(ndata.8.num.species.no.bot,
                                      fit_resp  = ilink.gam.8.num.species.no.bot(fit_link),
                                      right_upr = ilink.gam.8.num.species.no.bot(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.8.num.species.no.bot(fit_link - (2 * se_link)))


ndata.8.num.species.no.bot$min.10.pH.unscaled<-ndata.8.num.species.no.bot$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.num.species.no.bot.8 <- ggplot(ndata.8.num.species.no.bot, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = num.species.no.bot, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Native species richness")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.species.no.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,20)
plt.inv.num.species.no.bot.8
ggsave("C:Graphs August 2020//native_richness.8.png")


# Occupied space 16 ----------------------------------------------------------

gam.16.lm.native.occupied.space<- gam(native.occupied.space ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.gamma.native.occupied.space.1<- gam(native.occupied.space*0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")
gam.16.binomial.native.occupied.space<- gam(formula = cbind(native.occupied.space, 100-native.occupied.space)~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

gam.16.beta.native.occupied.space<- gam(native.occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.16.beta.native.occupied.space.1<- gam(native.occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.16.beta.native.occupied.space.2<- gam(native.occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.16.beta.native.occupied.space.3<- gam(native.occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.16_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.16.lm.native.occupied.space, gam.16.beta.native.occupied.space.3, gam.16.gamma.native.occupied.space.1, gam.16.beta.native.occupied.space, gam.16.beta.native.occupied.space.1, gam.16.beta.native.occupied.space.2, gam.16.binomial.native.occupied.space)


plot(gam.16.beta.native.occupied.space.3 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.beta.native.occupied.space.3 )
#looks good
qq_plot(gam.16.beta.native.occupied.space.3 , method = 'simulate')
k.check(gam.16.beta.native.occupied.space.3 )
summary(gam.16.beta.native.occupied.space.3 )
gam.16.beta.native.occupied.space.3.unordered<- gam(native.occupied.space.001~ s(min.10.pH, k=15)+ Invasives + s(min.10.pH, by=oInvasives, k=15), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.16.beta.native.occupied.space.3.unordered)

fam.gam.16.native.occupied.space <- family(gam.16.beta.native.occupied.space.3 )
fam.gam.16.native.occupied.space
str(fam.gam.16.native.occupied.space)
ilink.gam.16.native.occupied.space<- fam.gam.16.native.occupied.space$linkinv
ilink.gam.16.native.occupied.space


mod.native.occupied.space<-gam.16.beta.native.occupied.space.3 
ndata.16.native.occupied.space <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.native.occupied.space <- add_column(ndata.16.native.occupied.space, fit = predict(mod.native.occupied.space, newdata = ndata.16.native.occupied.space, type = 'response'))

predict(mod.native.occupied.space, newdata = ndata.16.native.occupied.space, type = 'response')
ndata.16.native.occupied.space <- bind_cols(ndata.16.native.occupied.space, setNames(as_tibble(predict(mod.native.occupied.space, ndata.16.native.occupied.space, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.native.occupied.space <- mutate(ndata.16.native.occupied.space,
                         fit_resp  = ilink.gam.16.native.occupied.space(fit_link),
                         right_upr = ilink.gam.16.native.occupied.space(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.16.native.occupied.space(fit_link - (2 * se_link)))


ndata.16.native.occupied.space$min.10.pH.unscaled<-ndata.16.native.occupied.space$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.native.occupied.space.16 <- ggplot(ndata.16.native.occupied.space, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = native.occupied.space.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.native.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.native.occupied.space.16
ggsave("C:Graphs August 2020//native.occupied.space_pred.16.png")




# Occupied space 8 ----------------------------------------------------------

gam.8.lm.native.occupied.space<- gam(native.occupied.space ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
gam.8.gamma.native.occupied.space.1<- gam(native.occupied.space*0.01 ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, family = Gamma, select=TRUE, method="REML")
gam.8.binomial.native.occupied.space<- gam(formula = cbind(native.occupied.space, 100-native.occupied.space)~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

gam.8.beta.native.occupied.space<- gam(native.occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.8.beta.native.occupied.space.1<- gam(native.occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.8.beta.native.occupied.space.2<- gam(native.occupied.space.001~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.8.beta.native.occupied.space.3<- gam(native.occupied.space.001~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives), data = invasion.exp.data.8_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.8.lm.native.occupied.space, gam.8.beta.native.occupied.space.3, gam.8.gamma.native.occupied.space.1, gam.8.beta.native.occupied.space, gam.8.beta.native.occupied.space.1, gam.8.beta.native.occupied.space.2, gam.8.binomial.native.occupied.space)
#beta is best 

plot(gam.8.beta.native.occupied.space, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.beta.native.occupied.space)
#looks good
qq_plot(gam.8.beta.native.occupied.space, method = 'simulate')
k.check(gam.8.beta.native.occupied.space)
summary(gam.8.beta.native.occupied.space)
gam.8.beta.native.occupied.space.unordered<- gam(native.occupied.space.001~ s(min.10.pH, k=15)+ Invasives + s(min.10.pH, by=oInvasives, k=15), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.8.beta.native.occupied.space.unordered)

fam.gam.8.native.occupied.space <- family(gam.8.beta.native.occupied.space)
fam.gam.8.native.occupied.space
str(fam.gam.8.native.occupied.space)
ilink.gam.8.native.occupied.space<- fam.gam.8.native.occupied.space$linkinv
ilink.gam.8.native.occupied.space


mod.native.occupied.space<-gam.8.beta.native.occupied.space
ndata.8.native.occupied.space <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.8.native.occupied.space <- add_column(ndata.8.native.occupied.space, fit = predict(mod.native.occupied.space, newdata = ndata.8.native.occupied.space, type = 'response'))

predict(mod.native.occupied.space, newdata = ndata.8.native.occupied.space, type = 'response')
ndata.8.native.occupied.space <- bind_cols(ndata.8.native.occupied.space, setNames(as_tibble(predict(mod.native.occupied.space, ndata.8.native.occupied.space, se.fit = TRUE)[1:2]),
                                                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.native.occupied.space <- mutate(ndata.8.native.occupied.space,
                                         fit_resp  = ilink.gam.8.native.occupied.space(fit_link),
                                         right_upr = ilink.gam.8.native.occupied.space(fit_link + (2 * se_link)),
                                         right_lwr = ilink.gam.8.native.occupied.space(fit_link - (2 * se_link)))


ndata.8.native.occupied.space$min.10.pH.unscaled<-ndata.8.native.occupied.space$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.native.occupied.space.8 <- ggplot(ndata.8.native.occupied.space, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y = native.occupied.space.001, shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.native.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.native.occupied.space.8
ggsave("C:Graphs August 2020//native.occupied.space_pred.8.png")


# CAP1 - 16 --------------------------------------------------------------------
#negative values so can't do gamma
gam.16.lm.CAP1.inv<- gam(CAP1.inv ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.CAP1.inv.1<- gam(CAP1.inv ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( gam.16.loglink.CAP1.inv.1, gam.16.lm.CAP1.inv)
#gam.16.lm.CAP1.inv

plot(gam.16.lm.CAP1.inv, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.CAP1.inv)
#look good
qq_plot(gam.16.lm.CAP1.inv, method = 'simulate')
k.check(gam.16.lm.CAP1.inv)
summary(gam.16.lm.CAP1.inv)
gam.16.lm.CAP1.inv.unordered<- gam(CAP1.inv ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
summary(gam.16.lm.CAP1.inv.unordered)

fam.gam.16.CAP1.inv <- family(gam.16.lm.CAP1.inv)
fam.gam.16.CAP1.inv
str(fam.gam.16.CAP1.inv)
ilink.gam.16.CAP1.inv<- fam.gam.16.CAP1.inv$linkinv
ilink.gam.16.CAP1.inv


mod.CAP1.inv<-gam.16.lm.CAP1.inv
ndata.16.CAP1.inv <- with(invasion.exp.data.16_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.16.CAP1.inv <- add_column(ndata.16.CAP1.inv, fit = predict(mod.CAP1.inv, newdata = ndata.16.CAP1.inv, type = 'response'))

predict(mod.CAP1.inv, newdata = ndata.16.CAP1.inv, type = 'response')
ndata.16.CAP1.inv <- bind_cols(ndata.16.CAP1.inv, setNames(as_tibble(predict(mod.CAP1.inv, ndata.16.CAP1.inv, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.16.CAP1.inv <- mutate(ndata.16.CAP1.inv,
                                      fit_resp  = ilink.gam.16.CAP1.inv(fit_link),
                                      right_upr = ilink.gam.16.CAP1.inv(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.16.CAP1.inv(fit_link - (2 * se_link)))


ndata.16.CAP1.inv$min.10.pH.unscaled<-ndata.16.CAP1.inv$min.10.pH * attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.CAP1.inv.16 <- ggplot(ndata.16.CAP1.inv, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(CAP1.inv), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(58% of constrained variation)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.CAP1.inv,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.1, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.inv.CAP1.inv.16
ggsave("C:Graphs August 2020//CAP1.inv_pred.16.png")

# CAP1.inv - 8 --------------------------------------------------------------------
#negative values so can't do gamma
gam.8.lm.CAP1.inv<- gam(CAP1.inv ~ s(min.10.pH)+ oInvasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")

AICtab( gam.8.loglink.CAP1.inv.1, gam.8.lm.CAP1.inv)
#gam.8.lm.CAP1.inv

plot(gam.8.lm.CAP1.inv, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.lm.CAP1.inv)
#look good
qq_plot(gam.8.lm.CAP1.inv, method = 'simulate')
k.check(gam.8.lm.CAP1.inv)
summary(gam.8.lm.CAP1.inv)
gam.8.lm.CAP1.inv.unordered<- gam(CAP1.inv ~ s(min.10.pH)+ Invasives + s(min.10.pH, by=oInvasives),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
summary(gam.8.lm.CAP1.inv.unordered)

fam.gam.8.CAP1.inv <- family(gam.8.lm.CAP1.inv)
fam.gam.8.CAP1.inv
str(fam.gam.8.CAP1.inv)
ilink.gam.8.CAP1.inv<- fam.gam.8.CAP1.inv$linkinv
ilink.gam.8.CAP1.inv


mod.CAP1.inv<-gam.8.lm.CAP1.inv
ndata.8.CAP1.inv <- with(invasion.exp.data.8_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                               length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.8.CAP1.inv <- add_column(ndata.8.CAP1.inv, fit = predict(mod.CAP1.inv, newdata = ndata.8.CAP1.inv, type = 'response'))

predict(mod.CAP1.inv, newdata = ndata.8.CAP1.inv, type = 'response')
ndata.8.CAP1.inv <- bind_cols(ndata.8.CAP1.inv, setNames(as_tibble(predict(mod.CAP1.inv, ndata.8.CAP1.inv, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.CAP1.inv <- mutate(ndata.8.CAP1.inv,
                        fit_resp  = ilink.gam.8.CAP1.inv(fit_link),
                        right_upr = ilink.gam.8.CAP1.inv(fit_link + (2 * se_link)),
                        right_lwr = ilink.gam.8.CAP1.inv(fit_link - (2 * se_link)))


ndata.8.CAP1.inv$min.10.pH.unscaled<-ndata.8.CAP1.inv$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.CAP1.inv.8 <- ggplot(ndata.8.CAP1.inv, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(CAP1.inv), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(61% of constrained variation)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.CAP1.inv,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.1, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.inv.CAP1.inv.8
ggsave("C:Graphs August 2020//CAP1.inv_pred.8.png")



# Distances 16 ---------------------------------------------------------------

gam.16.lm.distances<- gam(distcentroid ~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.distances.1<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.16.gamma.distances<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.16.loglink.distances.1,  gam.16.lm.distances, gam.16.gamma.distances)
#gam.16.lm.distances although both are equal


plot(gam.16.lm.distances, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.16.lm.distances)
#looks good
qq_plot(gam.16.lm.distances, method = 'simulate')
k.check(gam.16.lm.distances)
#good
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
plt.inv.distances.16 <- ggplot(ndata.16.distances, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(distcentroid), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.16_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.16.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.distances.16
ggsave("C:Graphs August 2020//distances_pred.16.png")

# Distances 8 ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.8.lm.distances<- gam(distcentroid ~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
gam.8.loglink.distances.1<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.8_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.8.gamma.distances<- gam(distcentroid~ s(min.10.pH, k=11)+ oInvasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.8_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.8.loglink.distances.1,  gam.8.lm.distances, gam.8.gamma.distances)
#gam.8.lm.distances although both are equal


plot(gam.8.lm.distances, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.8.lm.distances)
#looks good
qq_plot(gam.8.lm.distances, method = 'simulate')
k.check(gam.8.lm.distances)
#good
summary(gam.8.lm.distances)

gam.8.lm.distances.unordered<- gam(distcentroid~ s(min.10.pH, k=11)+ Invasives + s(min.10.pH, by=oInvasives, k=11),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
summary(gam.8.lm.distances.unordered)

fam.gam.8.distances <- family(gam.8.lm.distances)
fam.gam.8.distances
str(fam.gam.8.distances)
ilink.gam.8.distances<- fam.gam.8.distances$linkinv
ilink.gam.8.distances


mod.distances<-gam.8.lm.distances
ndata.8.distances <- with(invasion.exp.data.8_zscores, tibble(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oInvasives = oInvasives[want],  CO2.Treatment= CO2.Treatment[want]))


## add the fitted values by predicting from the model for the new data
ndata.8.distances <- add_column(ndata.8.distances, fit = predict(mod.distances, newdata = ndata.8.distances, type = 'response'))

predict(mod.distances, newdata = ndata.8.distances, type = 'response')
ndata.8.distances <- bind_cols(ndata.8.distances, setNames(as_tibble(predict(mod.distances, ndata.8.distances, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.8.distances <- mutate(ndata.8.distances,
                             fit_resp  = ilink.gam.8.distances(fit_link),
                             right_upr = ilink.gam.8.distances(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.8.distances(fit_link - (2 * se_link)))


ndata.8.distances$min.10.pH.unscaled<-ndata.8.distances$min.10.pH * attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$min.10.pH, 'scaled:center')


# plot 
plt.inv.distances.8 <- ggplot(ndata.8.distances, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oInvasives)) +
  geom_point(aes(y =(distcentroid), shape=CO2.Treatment, colour=oInvasives), data = invasion.exp.data.8_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_invasives, guide = guide_legend(title="Invasives", title.position = "top"))+
  scale_fill_manual(values=colorset_invasives, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH Invasives", title.position = "top"))+
  geom_ribbon(data = ndata.8.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oInvasives), alpha = 0.10)+
  theme(legend.position='none')
plt.inv.distances.8
ggsave("C:Graphs August 2020//distances_pred.8.png")



# Community plotting ------------------------------------------------------

#### community fig week 16
fig.community.week.16<-wrap_plots(plt.inv.native.occupied.space.16,
                            plt.inv.num.species.no.bot.16,
                            plt.inv.CAP1.inv.16, plt.inv.distances.16, ncol=2)+
                            plot_annotation(tag_levels = 'a')

fig.community.week.16

ggplot2::ggsave(plot=fig.community.week.16, "C:Graphs August 2020//Fig.inv.community.week16.pdf", width=4, height=3, units="in")



#### community fig week 8
fig.community.week.8<-wrap_plots(plt.inv.native.occupied.space.8,
                                  plt.inv.num.species.no.bot.8,
                                  plt.inv.CAP1.inv.8, plt.inv.distances.8, ncol=2)+
  plot_annotation(tag_levels = 'a')

fig.community.week.8

ggplot2::ggsave(plot=fig.community.week.8, "C:Graphs August 2020//Fig.inv.community.week8.pdf", width=4, height=3, units="in")



# Community level 16 tables --------------------------------------------------

num.species.no.bot.gam.16<- summary(gam.16.poisson.num.species.no.bot)
num.species.no.bot.gam.16.unordered<- summary(gam.16.poisson.num.species.no.bot.unordered)

occupied.space.gam.16<- summary(gam.16.beta.native.occupied.space.3)
occupied.space.gam.16.unordered<- summary(gam.16.beta.native.occupied.space.3.unordered)

distances.gam.16<- summary(gam.16.lm.distances)
distances.gam.16.unordered<- summary(gam.16.lm.distances.unordered)

CAP1.gam.16 <- summary(gam.16.lm.CAP1)
CAP1.gam.16.unordered <- summary(gam.16.lm.CAP1.unordered)

#ptable building
num.species.no.bot.gam.16.p.table<-as.data.frame(num.species.no.bot.gam.16.unordered$p.table)
num.species.no.bot.gam.16.s.table<-as.data.frame(num.species.no.bot.gam.16$s.table)

occupied.space.gam.16.p.table<-as.data.frame(occupied.space.gam.16.unordered$p.table)
occupied.space.gam.16.s.table<-as.data.frame(occupied.space.gam.16$s.table)

distances.gam.16.p.table<-as.data.frame(distances.gam.16.unordered$p.table)
distances.gam.16.s.table<-as.data.frame(distances.gam.16$s.table)

CAP1.gam.16.p.table<-as.data.frame(CAP1.gam.16.unordered$p.table)
CAP1.gam.16.s.table<-as.data.frame(CAP1.gam.16$s.table)

#num.species.no.bot.gam.16.p.table and  hydtobot.gam.16.p.table, is with z value 
colnames(num.species.no.bot.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(occupied.space.gam.16.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.community.t.16<-rbind(num.species.no.bot.gam.16.p.table,
              occupied.space.gam.16.p.table,
              CAP1.gam.16.p.table,
              distances.gam.16.p.table
              )


colnames(ptable.community.t.16) <- c("Estimate", "SE", "t", "p")
ptable.community.t.16$Factor<-rep(c("Intercept", "Invasives Present"))

ptable.community.t.16 %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (z)", 1,2) %>% 
  group_rows("Occupied space, beta (z)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//ptable.community.t.16.html", self_contained = T)


#num.species.no.bot.gam.16.p.table and  hydtobot.gam.16.p.table, is with Chisq
colnames(num.species.no.bot.gam.16.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(occupied.space.gam.16.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.community.f.16<-rbind(num.species.no.bot.gam.16.s.table,
                          occupied.space.gam.16.s.table,
                          CAP1.gam.16.s.table,
                          distances.gam.16.s.table)


colnames(stable.community.f.16) <- c("Estimated_df", "Reference_df", "F", "p_smooth")
stable.community.f.16$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Invasives Present"))

stable.community.f.16 %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//stable.community.f.16.html", self_contained = T)

pstable.community.16<-cbind(ptable.community.t.16, stable.community.f.16)


pstable.community.16 %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.051, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth, Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//pstable.community.16.html ", self_contained = T)

# Community level 8 tables --------------------------------------------------

num.species.no.bot.gam.8<- summary(gam.8.poisson.num.species.no.bot)
num.species.no.bot.gam.8.unordered<- summary(gam.8.poisson.num.species.no.bot.unordered)

occupied.space.gam.8<- summary(gam.8.beta.native.occupied.space)
occupied.space.gam.8.unordered<- summary(gam.8.beta.native.occupied.space.unordered)

distances.gam.8<- summary(gam.8.lm.distances)
distances.gam.8.unordered<- summary(gam.8.lm.distances.unordered)

CAP1.gam.8 <- summary(gam.8.lm.CAP1)
CAP1.gam.8.unordered <- summary(gam.8.lm.CAP1.unordered)

#ptable building
num.species.no.bot.gam.8.p.table<-as.data.frame(num.species.no.bot.gam.8.unordered$p.table)
num.species.no.bot.gam.8.s.table<-as.data.frame(num.species.no.bot.gam.8$s.table)

occupied.space.gam.8.p.table<-as.data.frame(occupied.space.gam.8.unordered$p.table)
occupied.space.gam.8.s.table<-as.data.frame(occupied.space.gam.8$s.table)

distances.gam.8.p.table<-as.data.frame(distances.gam.8.unordered$p.table)
distances.gam.8.s.table<-as.data.frame(distances.gam.8$s.table)

CAP1.gam.8.p.table<-as.data.frame(CAP1.gam.8.unordered$p.table)
CAP1.gam.8.s.table<-as.data.frame(CAP1.gam.8$s.table)

#num.species.no.bot.gam.8.p.table and  hydtobot.gam.8.p.table, is with z value 
colnames(num.species.no.bot.gam.8.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(occupied.space.gam.8.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.community.t.8<-rbind(num.species.no.bot.gam.8.p.table,
                             occupied.space.gam.8.p.table,
                             CAP1.gam.8.p.table,
                             distances.gam.8.p.table)


colnames(ptable.community.t.8) <- c("Estimate", "SE", "t", "p")
ptable.community.t.8$Factor<-rep(c("Intercept", "Invasives Present"))

ptable.community.t.8 %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (z)", 1,2) %>% 
  group_rows("Occupied space, beta (z)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//ptable.community.t.8.html", self_contained = T)


#num.species.no.bot.gam.8.p.table and  hydtobot.gam.8.p.table, is with Chisq
colnames(num.species.no.bot.gam.8.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(occupied.space.gam.8.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.community.f.8<-rbind(num.species.no.bot.gam.8.s.table,
                             occupied.space.gam.8.s.table,
                             CAP1.gam.8.s.table,
                             distances.gam.8.s.table)


colnames(stable.community.f.8) <- c("Estimated_df", "Reference_df", "F", "p_smooth")
stable.community.f.8$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Invasives Present"))

stable.community.f.8 %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//stable.community.f.8.html", self_contained = T)

pstable.community.8<-cbind(ptable.community.t.8, stable.community.f.8)


pstable.community.8 %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.051, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth, Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,6) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//pstable.community.8.html ", self_contained = T)
