
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
# invasion.exp.data$num.nudi<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.eggs+invasion.exp.data$nudi.hatched

invasion.exp.data$bot.total<-invasion.exp.data$botryllid + invasion.exp.data$bot.eaten
invasion.exp.data$mem.total<-invasion.exp.data$membranipora + invasion.exp.data$mem.eaten + invasion.exp.data$mem.dead
invasion.exp.data$corella.total<-invasion.exp.data$corella + invasion.exp.data$dead.corella

invasion.exp.data$prop.mem.dead<-invasion.exp.data$mem.dead/invasion.exp.data$membranipora
invasion.exp.data$prop.mem.eaten<-invasion.exp.data$mem.eaten/invasion.exp.data$membranipora


invasion.exp.data$hydroid.001<-(0.01*(invasion.exp.data$hydroid))+0.01
invasion.exp.data$botryllid.001<-(0.01*(invasion.exp.data$botryllid))+0.01
invasion.exp.data$bot.eaten.001<-(0.01*(invasion.exp.data$bot.eaten))+0.01
invasion.exp.data$mem.eaten.001<-(0.01*(invasion.exp.data$mem.eaten))+0.01
invasion.exp.data$mem.total.001<-(0.01*(invasion.exp.data$mem.total))+0.01

invasion.exp.data$membranipora.001<-(0.01*(invasion.exp.data$membranipora))+0.01
invasion.exp.data$mussel.001<-(0.01*(invasion.exp.data$mussel))+0.01
invasion.exp.data$didemnum<-invasion.exp.data$white.bryo+invasion.exp.data$fan.bryo
invasion.exp.data$folliculina<-invasion.exp.data$protozoa
invasion.exp.data$folliculina.001<-(0.01*(invasion.exp.data$folliculina))+0.01
invasion.exp.data$didemnum.001<-(0.01*(invasion.exp.data$didemnum))+0.01
invasion.exp.data$occupied.space<-(100 - invasion.exp.data$bare)
invasion.exp.data$occupied.space.001<-(0.01*(invasion.exp.data$occupied.space))+0.01
invasion.exp.data$native.occupied.space<-(100 - invasion.exp.data$botryllid  -invasion.exp.data$bot.eaten -  invasion.exp.data$bare)
invasion.exp.data$native.occupied.space.001<-(0.01*(invasion.exp.data$native.occupied.space))+0.01

invasion.exp.data$prop.mem.dead.001<-(0.01*invasion.exp.data$prop.mem.dead)+0.01
invasion.exp.data$prop.mem.eaten.001<-(0.01*invasion.exp.data$prop.mem.eaten)+0.01


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
#invasion.exp.data.16_zscores$botryllid.001<-scale(invasion.exp.data.16$botryllid.001, center=TRUE, scale=TRUE)
invasion.exp.data.16_zscores$Mesocosm <- as.factor(invasion.exp.data.16$Mesocosm)
invasion.exp.data.16_zscores$av.pH.unscaled <-invasion.exp.data.16_zscores$av.pH * attr(invasion.exp.data.16_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.16_zscores$av.pH, 'scaled:center')
#invasion.exp.data.16_zscores$botryllid.001.unscaled <-invasion.exp.data.16_zscores$botryllid.001 * attr(invasion.exp.data.16_zscores$botryllid.001, 'scaled:scale') + attr(invasion.exp.data.16_zscores$botryllid.001, 'scaled:center')


# for visualizing histograms of pH to use as continuous vs. discrete variable
invasion.exp.data.16_pres<-invasion.exp.data.16 %>% filter(Invasives=="Present")
invasion.exp.data.16_abs<-invasion.exp.data.16 %>% filter(Invasives=="Absent")


invasion.exp.data.8.community<-read.csv("C:Biological data/invasion.exp.data.8.community.csv")
invasion.exp.data.8<-merge(invasion.exp.data.8, invasion.exp.data.8.community, by="Tile.ID")
# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
invasion.exp.data.8_zscores<-invasion.exp.data.8
#invasion.exp.data.8_zscores$hydrogen.concentration<-scale(invasion.exp.data.8$hydrogen.concentration, center=TRUE, scale=TRUE)
invasion.exp.data.8_zscores$av.pH<-scale(invasion.exp.data.8$av.pH, center=TRUE, scale=TRUE)
#invasion.exp.data.8_zscores$botryllid.001<-scale(invasion.exp.data.8$botryllid.001, center=TRUE, scale=TRUE)
invasion.exp.data.8_zscores$Mesocosm <- as.factor(invasion.exp.data.8$Mesocosm)
invasion.exp.data.8_zscores$av.pH.unscaled <-invasion.exp.data.8_zscores$av.pH * attr(invasion.exp.data.8_zscores$av.pH, 'scaled:scale') + attr(invasion.exp.data.8_zscores$av.pH, 'scaled:center')
#invasion.exp.data.8_zscores$botryllid.001.unscaled <-invasion.exp.data.8_zscores$botryllid.001 * attr(invasion.exp.data.8_zscores$botryllid.001, 'scaled:scale') + attr(invasion.exp.data.8_zscores$botryllid.001, 'scaled:center')


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

#Used code from 

# Plotting settings -------------------------------------------------------

colorset_invasives = c("Present"="#A20228" ,"Absent"="#818392")
colorset_CO2.Treatment = c("CO2"="#A20228" ,"AIR"="#818392")

theme_set(theme_classic(base_size = 12))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

 

# GAM beta folliculina -----------------------------------------------------------

gam.8.binomial.folliculina<- gam(formula = cbind(folliculina, 100-folliculina)~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.folliculina<- gam(folliculina.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#beta explains more deviance

mod.folliculina<-gam.8.beta.folliculina
plot(mod.folliculina, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.folliculina)
qq_plot(mod.folliculina, method = 'simulate')
k.check(mod.folliculina)
summary(mod.folliculina)

mod.folliculina.unordered<- gam(folliculina.001~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")

x1<-list(seq(min(invasion.exp.data.8_zscores$bot.total),max(invasion.exp.data.8_zscores$bot.total),length=100))
x2<-list(unique(invasion.exp.data.8_zscores$CO2.Treatment))
all_var<-x1
all_var<-c(x1,x2)
#expand.grid on it
all_var<-expand.grid(all_var)
names(all_var)[1]<-"bot.total"
names(all_var)[2]<-"CO2.Treatment"
all_var$oCO2.Treatment<-all_var$CO2.Treatment
ndata.8.folliculina<-as_tibble(all_var)

fam.gam.8.folliculina <- family(mod.folliculina)
ilink.gam.8.folliculina<- fam.gam.8.folliculina$linkinv

## add the fitted values by predicting from the mod.folliculinael for the new data
ndata.8.folliculina <- add_column(ndata.8.folliculina, fit = predict(mod.folliculina, newdata = ndata.8.folliculina, type = 'response'))
ndata.8.folliculina <- bind_cols(ndata.8.folliculina, setNames(as_tibble(predict(mod.folliculina, ndata.8.folliculina, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))
ndata.8.folliculina <- mutate(ndata.8.folliculina,
                               fit_resp  = ilink.gam.8.folliculina(fit_link),
                               right_upr = ilink.gam.8.folliculina(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.folliculina(fit_link - (2 * se_link)))

plt.folliculina.8 <- ggplot(ndata.8.folliculina, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = folliculina.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.folliculina,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.folliculina.8
ggsave("C:Graphs August 2020//folliculina_pred.8.png")


# GAM binomial membranipora --------------------------------------------------------
gam.8.binomial.membranipora<- gam(formula = cbind(membranipora, 100-membranipora)~ s(bot.total, k=6)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment, k=6), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.membranipora<- gam(membranipora.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial explains more deviance 40 vs. 12%

summary(gam.8.beta.membranipora)

mod.membranipora<-gam.8.binomial.membranipora
plot(mod.membranipora, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.membranipora)
qq_plot(mod.membranipora, method = 'simulate')
k.check(mod.membranipora)
summary(mod.membranipora)

mod.membranipora.unordered<-gam(formula = cbind(membranipora, 100-membranipora)~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

ndata.8.membranipora<-as_tibble(all_var)
fam.gam.8.membranipora <- family(mod.membranipora)
ilink.gam.8.membranipora<- fam.gam.8.membranipora$linkinv

## add the fitted values by predicting from the mod.membraniporael for the new data
ndata.8.membranipora <- add_column(ndata.8.membranipora, fit = predict(mod.membranipora, newdata = ndata.8.membranipora, type = 'response'))
ndata.8.membranipora <- bind_cols(ndata.8.membranipora, setNames(as_tibble(predict(mod.membranipora, ndata.8.membranipora, se.fit = TRUE)[1:2]),
                                                                   c('fit_link','se_link')))

ndata.8.membranipora <- mutate(ndata.8.membranipora,
                                fit_resp  = ilink.gam.8.membranipora(fit_link),
                                right_upr = ilink.gam.8.membranipora(fit_link + (2 * se_link)),
                                right_lwr = ilink.gam.8.membranipora(fit_link - (2 * se_link)))

plt.membranipora.8 <- ggplot(ndata.8.membranipora, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = membranipora.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.membranipora,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.membranipora.8
ggsave("C:Graphs August 2020//membranipora_pred.8.png")


# GAM binomial mem.total --------------------------------------------------------
gam.8.binomial.mem.total<- gam(formula = cbind(mem.total, 100-mem.total)~ s(bot.total, k=6)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment, k=6), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.mem.total<- gam(mem.total.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial explains more deviance 40 vs. 12%

summary(gam.8.beta.mem.total)

mod.mem.total<-gam.8.binomial.mem.total
plot(mod.mem.total, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.mem.total)
qq_plot(mod.mem.total, method = 'simulate')
k.check(mod.mem.total)
summary(mod.mem.total)

mod.mem.total.unordered<-gam(formula = cbind(mem.total, 100-mem.total)~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

ndata.8.mem.total<-as_tibble(all_var)
fam.gam.8.mem.total <- family(mod.mem.total)
ilink.gam.8.mem.total<- fam.gam.8.mem.total$linkinv

## add the fitted values by predicting from the mod.mem.totalel for the new data
ndata.8.mem.total <- add_column(ndata.8.mem.total, fit = predict(mod.mem.total, newdata = ndata.8.mem.total, type = 'response'))
ndata.8.mem.total <- bind_cols(ndata.8.mem.total, setNames(as_tibble(predict(mod.mem.total, ndata.8.mem.total, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.mem.total <- mutate(ndata.8.mem.total,
                               fit_resp  = ilink.gam.8.mem.total(fit_link),
                               right_upr = ilink.gam.8.mem.total(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.mem.total(fit_link - (2 * se_link)))

plt.mem.total.8 <- ggplot(ndata.8.mem.total, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = mem.total.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("mem.total")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.mem.total,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.mem.total.8
ggsave("C:Graphs August 2020//mem.total_pred.8.png")



# Prop mem eaten ----------------------------------------------------------
gam.8.binomial.prop.mem.eaten<- gam(formula = cbind(prop.mem.eaten, 100-prop.mem.eaten)~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.prop.mem.eaten<- gam(prop.mem.eaten.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial explains more deviance 40 vs. 12%

mod.prop.mem.eaten<-gam.8.binomial.prop.mem.eaten
plot(mod.prop.mem.eaten, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.prop.mem.eaten)
qq_plot(mod.prop.mem.eaten, method = 'simulate')
k.check(mod.prop.mem.eaten)
summary(mod.prop.mem.eaten)

mod.prop.mem.eaten.unordered<-gam(formula = cbind(prop.mem.eaten, 100-prop.mem.eaten)~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

ndata.8.prop.mem.eaten<-as_tibble(all_var)
fam.gam.8.prop.mem.eaten <- family(mod.prop.mem.eaten)
ilink.gam.8.prop.mem.eaten<- fam.gam.8.prop.mem.eaten$linkinv

## add the fitted values by predicting from the mod.prop.mem.eatenel for the new data
ndata.8.prop.mem.eaten <- add_column(ndata.8.prop.mem.eaten, fit = predict(mod.prop.mem.eaten, newdata = ndata.8.prop.mem.eaten, type = 'response'))
ndata.8.prop.mem.eaten <- bind_cols(ndata.8.prop.mem.eaten, setNames(as_tibble(predict(mod.prop.mem.eaten, ndata.8.prop.mem.eaten, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

ndata.8.prop.mem.eaten <- mutate(ndata.8.prop.mem.eaten,
                            fit_resp  = ilink.gam.8.prop.mem.eaten(fit_link),
                            right_upr = ilink.gam.8.prop.mem.eaten(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.8.prop.mem.eaten(fit_link - (2 * se_link)))

plt.prop.mem.eaten.8 <- ggplot(ndata.8.prop.mem.eaten, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = prop.mem.eaten.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("prop.mem.eaten")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.prop.mem.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.prop.mem.eaten.8
ggsave("C:Graphs August 2020//prop.mem.eaten_pred.8.png")


# Week 16 prop mem eaten --------------------------------------------------

gam.16.binomial.prop.mem.eaten<- gam(formula = cbind(prop.mem.eaten, 100-prop.mem.eaten)~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")
gam.16.beta.prop.mem.eaten<- gam(prop.mem.eaten.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")

mod.prop.mem.eaten<-gam.16.beta.prop.mem.eaten
plot(mod.prop.mem.eaten, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.prop.mem.eaten)
qq_plot(mod.prop.mem.eaten, method = 'simulate')
k.check(mod.prop.mem.eaten)
summary(mod.prop.mem.eaten)

mod.prop.mem.eaten.unordered<-gam(formula = cbind(prop.mem.eaten, 100-prop.mem.eaten)~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")

ndata.16.prop.mem.eaten<-as_tibble(all_var)
fam.gam.16.prop.mem.eaten <- family(mod.prop.mem.eaten)
ilink.gam.16.prop.mem.eaten<- fam.gam.16.prop.mem.eaten$linkinv

## add the fitted values by predicting from the mod.prop.mem.eatenel for the new data
ndata.16.prop.mem.eaten <- add_column(ndata.16.prop.mem.eaten, fit = predict(mod.prop.mem.eaten, newdata = ndata.16.prop.mem.eaten, type = 'response'))
ndata.16.prop.mem.eaten <- bind_cols(ndata.16.prop.mem.eaten, setNames(as_tibble(predict(mod.prop.mem.eaten, ndata.16.prop.mem.eaten, se.fit = TRUE)[1:2]),
                                                                     c('fit_link','se_link')))

ndata.16.prop.mem.eaten <- mutate(ndata.16.prop.mem.eaten,
                                 fit_resp  = ilink.gam.16.prop.mem.eaten(fit_link),
                                 right_upr = ilink.gam.16.prop.mem.eaten(fit_link + (2 * se_link)),
                                 right_lwr = ilink.gam.16.prop.mem.eaten(fit_link - (2 * se_link)))

plt.prop.mem.eaten.16 <- ggplot(ndata.16.prop.mem.eaten, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = prop.mem.eaten.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("prop.mem.eaten")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.prop.mem.eaten,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.prop.mem.eaten.16
ggsave("C:Graphs August 2020//prop.mem.eaten_pred.16.png")




# GAM binomial mussel --------------------------------------------------------
gam.8.binomial.mussel<- gam(formula = cbind(mussel, 100-mussel)~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.mussel<- gam(mussel.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial explains more deviance 4 vs. 1%

mod.mussel<-gam.8.binomial.mussel
plot(mod.mussel, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.mussel)
qq_plot(mod.mussel, method = 'simulate')
k.check(mod.mussel)
summary(mod.mussel)

mod.mussel.unordered<-gam(formula = cbind(mussel, 100-mussel)~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")

ndata.8.mussel<-as_tibble(all_var)
fam.gam.8.mussel <- family(mod.mussel)
ilink.gam.8.mussel<- fam.gam.8.mussel$linkinv

## add the fitted values by predicting from the mod.musselel for the new data
ndata.8.mussel <- add_column(ndata.8.mussel, fit = predict(mod.mussel, newdata = ndata.8.mussel, type = 'response'))
ndata.8.mussel <- bind_cols(ndata.8.mussel, setNames(as_tibble(predict(mod.mussel, ndata.8.mussel, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.mussel <- mutate(ndata.8.mussel,
                               fit_resp  = ilink.gam.8.mussel(fit_link),
                               right_upr = ilink.gam.8.mussel(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.mussel(fit_link - (2 * se_link)))

plt.mussel.8 <- ggplot(ndata.8.mussel, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = mussel.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("mussel")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.mussel,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.mussel.8
ggsave("C:Graphs August 2020//mussel_pred.8.png")


# GAM poisson barnacles -----------------------------------------------------------
nbinom.8.num.barn <- fitdistr(invasion.exp.data.8_zscores$num.barn, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.barn, "nbinom", size = nbinom.8.num.barn$estimate[[1]], mu = nbinom.8.num.barn$estimate[[2]])

#negative binomial 
gam.8.nb.num.barn<- gam(num.barn ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.barn$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.barn.1<- gam(num.barn ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.barn<- gam(num.barn ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.barn, gam.8.nb.num.barn.1, gam.8.poisson.num.barn)
#poisson is better

mod.num.barn<-gam.8.poisson.num.barn
plot(mod.num.barn, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.barn)
qq_plot(mod.num.barn, method = 'simulate')
k.check(mod.num.barn)
summary(mod.num.barn)

gam.8.poisson.num.barn.unordered<- gam(num.barn ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.barn<-as_tibble(all_var)
fam.gam.8.num.barn <- family(mod.num.barn)
ilink.gam.8.num.barn<- fam.gam.8.num.barn$linkinv

## add the fitted values by predicting from the mod.num.barnel for the new data
ndata.8.num.barn <- add_column(ndata.8.num.barn, fit = predict(mod.num.barn, newdata = ndata.8.num.barn, type = 'response'))
ndata.8.num.barn <- bind_cols(ndata.8.num.barn, setNames(as_tibble(predict(mod.num.barn, ndata.8.num.barn, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.barn <- mutate(ndata.8.num.barn,
                               fit_resp  = ilink.gam.8.num.barn(fit_link),
                               right_upr = ilink.gam.8.num.barn(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.barn(fit_link - (2 * se_link)))

plt.num.barn.8 <- ggplot(ndata.8.num.barn, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.barn, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Barnacle abundance"), textstyle("(# individuals)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.barn,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.barn.8
ggsave("C:Graphs August 2020//num.barn_pred.8.png")

# GAM negbin num.white.bryo / gam.8.nb.num.white.bryo ----------------------------------------------------------
nbinom.8.num.white.bryo <- fitdistr(invasion.exp.data.8_zscores$num.white.bryo, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.white.bryo, "nbinom", size = nbinom.8.num.white.bryo$estimate[[1]], mu = nbinom.8.num.white.bryo$estimate[[2]])

#negative binomial 
gam.8.nb.num.white.bryo<- gam(num.white.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.white.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.white.bryo.1<- gam(num.white.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.white.bryo<- gam(num.white.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.white.bryo, gam.8.nb.num.white.bryo.1, gam.8.poisson.num.white.bryo)


mod.num.white.bryo<-gam.8.poisson.num.white.bryo
plot(mod.num.white.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.white.bryo)
qq_plot(mod.num.white.bryo, method = 'simulate')
k.check(mod.num.white.bryo)
summary(mod.num.white.bryo)

gam.8.poisson.num.white.bryo.unordered<- gam(num.white.bryo ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.white.bryo<-as_tibble(all_var)
fam.gam.8.num.white.bryo <- family(mod.num.white.bryo)
ilink.gam.8.num.white.bryo<- fam.gam.8.num.white.bryo$linkinv

## add the fitted values by predicting from the mod.num.white.bryoel for the new data
ndata.8.num.white.bryo <- add_column(ndata.8.num.white.bryo, fit = predict(mod.num.white.bryo, newdata = ndata.8.num.white.bryo, type = 'response'))
ndata.8.num.white.bryo <- bind_cols(ndata.8.num.white.bryo, setNames(as_tibble(predict(mod.num.white.bryo, ndata.8.num.white.bryo, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))

ndata.8.num.white.bryo <- mutate(ndata.8.num.white.bryo,
                           fit_resp  = ilink.gam.8.num.white.bryo(fit_link),
                           right_upr = ilink.gam.8.num.white.bryo(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.8.num.white.bryo(fit_link - (2 * se_link)))

plt.num.white.bryo.8 <- ggplot(ndata.8.num.white.bryo, aes(x = bot.total, y = fit)) + 
    geom_line(aes(colour=oCO2.Treatment)) +
    geom_jitter(aes(y = num.white.bryo, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
    xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Disporella abundance"), textstyle("(# colonies)")))))+  
    scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
    scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ylim(0,3)+
    scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
    geom_ribbon(data = ndata.8.num.white.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
    theme(legend.position='none')
plt.num.white.bryo.8
ggsave("C:Graphs August 2020//num.white.bryo_pred.8.png")

# GAM negbin num.red.bryo / gam.8.nb.num.red.bryo --------------------------------------------------------------
nbinom.8.num.red.bryo <- fitdistr(invasion.exp.data.8_zscores$num.red.bryo, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.red.bryo, "nbinom", size = nbinom.8.num.red.bryo$estimate[[1]], mu = nbinom.8.num.red.bryo$estimate[[2]])

#negative binomial 
gam.8.nb.num.red.bryo<- gam(num.red.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.red.bryo$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.red.bryo.1<- gam(num.red.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.red.bryo<- gam(num.red.bryo ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.red.bryo, gam.8.nb.num.red.bryo.1, gam.8.poisson.num.red.bryo)

mod.num.red.bryo<-gam.8.poisson.num.red.bryo
plot(mod.num.red.bryo, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.red.bryo)
qq_plot(mod.num.red.bryo, method = 'simulate')
k.check(mod.num.red.bryo)
summary(mod.num.red.bryo)

gam.8.poisson.num.red.bryo.unordered<- gam(num.red.bryo ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.red.bryo<-as_tibble(all_var)
fam.gam.8.num.red.bryo <- family(mod.num.red.bryo)
ilink.gam.8.num.red.bryo<- fam.gam.8.num.red.bryo$linkinv

## add the fitted values by predicting from the mod.num.red.bryoel for the new data
ndata.8.num.red.bryo <- add_column(ndata.8.num.red.bryo, fit = predict(mod.num.red.bryo, newdata = ndata.8.num.red.bryo, type = 'response'))
ndata.8.num.red.bryo <- bind_cols(ndata.8.num.red.bryo, setNames(as_tibble(predict(mod.num.red.bryo, ndata.8.num.red.bryo, se.fit = TRUE)[1:2]),
                                                                     c('fit_link','se_link')))

ndata.8.num.red.bryo <- mutate(ndata.8.num.red.bryo,
                                 fit_resp  = ilink.gam.8.num.red.bryo(fit_link),
                                 right_upr = ilink.gam.8.num.red.bryo(fit_link + (2 * se_link)),
                                 right_lwr = ilink.gam.8.num.red.bryo(fit_link - (2 * se_link)))

plt.num.red.bryo.8 <- ggplot(ndata.8.num.red.bryo, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.red.bryo, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Schizoporella abundance"), textstyle("(# colonies)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.red.bryo,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.red.bryo.8
ggsave("C:Graphs August 2020//num.red.bryo_pred.8.png")


# GAM poisson num nudi / gam.8.poisson.num.nudi  ------------------------------------------------------------
nbinom.8.num.nudi <- fitdistr(invasion.exp.data.8_zscores$num.nudi, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.nudi, "nbinom", size = nbinom.8.num.nudi$estimate[[1]], mu = nbinom.8.num.nudi$estimate[[2]])

#negative binomial 
gam.8.nb.num.nudi<- gam(num.nudi ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.nudi.1<- gam(num.nudi ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.nudi<- gam(num.nudi ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.nudi, gam.8.nb.num.nudi.1, gam.8.poisson.num.nudi)
#nb is same as poisson

mod.num.nudi<-gam.8.poisson.num.nudi
plot(mod.num.nudi, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.nudi)
qq_plot(mod.num.nudi, method = 'simulate')
k.check(mod.num.nudi)
summary(mod.num.nudi)

gam.8.poisson.num.nudi.unordered<- gam(num.nudi ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")
ndata.8.num.nudi<-as_tibble(all_var)
fam.gam.8.num.nudi <- family(mod.num.nudi)
ilink.gam.8.num.nudi<- fam.gam.8.num.nudi$linkinv

## add the fitted values by predicting from the mod.num.nudiel for the new data
ndata.8.num.nudi <- add_column(ndata.8.num.nudi, fit = predict(mod.num.nudi, newdata = ndata.8.num.nudi, type = 'response'))
ndata.8.num.nudi <- bind_cols(ndata.8.num.nudi, setNames(as_tibble(predict(mod.num.nudi, ndata.8.num.nudi, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.nudi <- mutate(ndata.8.num.nudi,
                               fit_resp  = ilink.gam.8.num.nudi(fit_link),
                               right_upr = ilink.gam.8.num.nudi(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.nudi(fit_link - (2 * se_link)))

plt.num.nudi.8 <- ggplot(ndata.8.num.nudi, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.nudi, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Nudibranch abundance"), textstyle("(# nudibranchs)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.nudi.8
ggsave("C:Graphs August 2020//num.nudi_pred.8.png")

# GAM poisson num nudi / gam.8.poisson.num.nudi.eggs  ------------------------------------------------------------
nbinom.8.num.nudi.egg <- fitdistr(invasion.exp.data.8_zscores$num.nudi.egg, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.nudi.egg, "nbinom", size = nbinom.8.num.nudi.egg$estimate[[1]], mu = nbinom.8.num.nudi.egg$estimate[[2]])

#negative binomial 
gam.8.nb.num.nudi.egg<- gam(num.nudi.egg ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.nudi.egg$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.nudi.egg.1<- gam(num.nudi.egg ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.nudi.egg<- gam(num.nudi.egg ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.nudi.egg, gam.8.nb.num.nudi.egg.1, gam.8.poisson.num.nudi.egg)
#nb is same as poisson

mod.num.nudi.egg<-gam.8.poisson.num.nudi.egg
plot(mod.num.nudi.egg, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.nudi.egg)
qq_plot(mod.num.nudi.egg, method = 'simulate')
k.check(mod.num.nudi.egg)
summary(mod.num.nudi.egg)

gam.8.poisson.num.nudi.egg.unordered<- gam(num.nudi.egg ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.nudi.egg<-as_tibble(all_var)
fam.gam.8.num.nudi.egg <- family(mod.num.nudi.egg)
ilink.gam.8.num.nudi.egg<- fam.gam.8.num.nudi.egg$linkinv

## add the fitted values by predicting from the mod.num.nudi.eggel for the new data
ndata.8.num.nudi.egg <- add_column(ndata.8.num.nudi.egg, fit = predict(mod.num.nudi.egg, newdata = ndata.8.num.nudi.egg, type = 'response'))
ndata.8.num.nudi.egg <- bind_cols(ndata.8.num.nudi.egg, setNames(as_tibble(predict(mod.num.nudi.egg, ndata.8.num.nudi.egg, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))

ndata.8.num.nudi.egg <- mutate(ndata.8.num.nudi.egg,
                           fit_resp  = ilink.gam.8.num.nudi.egg(fit_link),
                           right_upr = ilink.gam.8.num.nudi.egg(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.8.num.nudi.egg(fit_link - (2 * se_link)))

plt.num.nudi.egg.8 <- ggplot(ndata.8.num.nudi.egg, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_point(aes(y = num.nudi.egg, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Nudibranch abundance"), textstyle("(# eggs)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.nudi.egg,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.nudi.egg.8
ggsave("C:Graphs August 2020//num.nudi.egg_pred.8.png")

# GAM poisson num nudi / gam.8.poisson.num.nudi.alls  ------------------------------------------------------------
invasion.exp.data.8_zscores$num.nudi.all<- invasion.exp.data.8_zscores$num.nudi+invasion.exp.data.8_zscores$num.nudi.egg

nbinom.8.num.nudi.all <- fitdistr(invasion.exp.data.8_zscores$num.nudi.all, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.nudi.all, "nbinom", size = nbinom.8.num.nudi.all$estimate[[1]], mu = nbinom.8.num.nudi.all$estimate[[2]])

#negative binomial 
gam.8.nb.num.nudi.all<- gam(num.nudi.all ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.nudi.all$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.nudi.all.1<- gam(num.nudi.all ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.nudi.all<- gam(num.nudi.all ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.nudi.all, gam.8.nb.num.nudi.all.1, gam.8.poisson.num.nudi.all)

mod.num.nudi.all<-gam.8.poisson.num.nudi.all
plot(mod.num.nudi.all, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.nudi.all)
qq_plot(mod.num.nudi.all, method = 'simulate')
k.check(mod.num.nudi.all)
summary(mod.num.nudi.all)

gam.8.poisson.num.nudi.all.unordered<- gam(num.nudi.all ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")
ndata.8.num.nudi.all<-as_tibble(all_var)
fam.gam.8.num.nudi.all <- family(mod.num.nudi.all)
ilink.gam.8.num.nudi.all<- fam.gam.8.num.nudi.all$linkinv

## add the fitted values by predicting from the mod.num.nudi.allel for the new data
ndata.8.num.nudi.all <- add_column(ndata.8.num.nudi.all, fit = predict(mod.num.nudi.all, newdata = ndata.8.num.nudi.all, type = 'response'))
ndata.8.num.nudi.all <- bind_cols(ndata.8.num.nudi.all, setNames(as_tibble(predict(mod.num.nudi.all, ndata.8.num.nudi.all, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.nudi.all <- mutate(ndata.8.num.nudi.all,
                               fit_resp  = ilink.gam.8.num.nudi.all(fit_link),
                               right_upr = ilink.gam.8.num.nudi.all(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.nudi.all(fit_link - (2 * se_link)))

plt.num.nudi.all.8 <- ggplot(ndata.8.num.nudi.all, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.nudi.all, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Nudibranch abundance"), textstyle("(# nudibranchs)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.nudi.all,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.nudi.all.8
ggsave("C:Graphs August 2020//num.nudi.all_pred.8.png")



# GAM nb() serpulids / gam.8.nb.num.serpulid.1 -----------------------------------------------------------
nbinom.8.num.serpulid <- fitdistr(invasion.exp.data.8_zscores$num.serpulid, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.serpulid, "nbinom", size = nbinom.8.num.serpulid$estimate[[1]], mu = nbinom.8.num.serpulid$estimate[[2]])

#negative binomial 
gam.8.nb.num.serpulid<- gam(num.serpulid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.serpulid.1<- gam(num.serpulid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.serpulid<- gam(num.serpulid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")


AICtab(gam.8.nb.num.serpulid, gam.8.nb.num.serpulid.1, gam.8.poisson.num.serpulid)
#nb is same as poisson


gam.8.nb.num.serpulid.te<- gam(num.serpulid ~ s(bot.total)+ s(min.10.pH) + te(bot.total, min.10.pH),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.serpulid.1<- gam(num.serpulid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.serpulid<- gam(num.serpulid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")


mod.num.serpulid<-gam.8.nb.num.serpulid
plot(mod.num.serpulid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.serpulid)
qq_plot(mod.num.serpulid, method = 'simulate')
k.check(mod.num.serpulid)
summary(mod.num.serpulid)

gam.8.nb.num.serpulid.unordered<- gam(num.serpulid ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.serpulid$estimate[[1]]), select=TRUE, method="REML")

ndata.8.num.serpulid<-as_tibble(all_var)
fam.gam.8.num.serpulid <- family(mod.num.serpulid)
ilink.gam.8.num.serpulid<- fam.gam.8.num.serpulid$linkinv

## add the fitted values by predicting from the mod.num.serpulidel for the new data
ndata.8.num.serpulid <- add_column(ndata.8.num.serpulid, fit = predict(mod.num.serpulid, newdata = ndata.8.num.serpulid, type = 'response'))
ndata.8.num.serpulid <- bind_cols(ndata.8.num.serpulid, setNames(as_tibble(predict(mod.num.serpulid, ndata.8.num.serpulid, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.serpulid <- mutate(ndata.8.num.serpulid,
                               fit_resp  = ilink.gam.8.num.serpulid(fit_link),
                               right_upr = ilink.gam.8.num.serpulid(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.serpulid(fit_link - (2 * se_link)))

plt.num.serpulid.8 <- ggplot(ndata.8.num.serpulid, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.serpulid, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# individuals)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.serpulid.8
ggsave("C:Graphs August 2020//num.serpulid_pred.8.png")


# GAM negbin corella / gam.16.nb.num.corella -------------------------------------------------------------
nbinom.8.num.corella <- fitdistr(invasion.exp.data.8_zscores$num.corella, "Negative Binomial")
qqp(invasion.exp.data.8_zscores$num.corella, "nbinom", size = nbinom.8.num.corella$estimate[[1]], mu = nbinom.8.num.corella$estimate[[2]])

#negative binomial 
gam.8.nb.num.corella<- gam(num.corella ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = negbin(nbinom.8.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.8.nb.num.corella.1<- gam(num.corella ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.corella<- gam(num.corella ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.8.nb.num.corella, gam.8.nb.num.corella.1, gam.8.poisson.num.corella)
#nb is same as poisson

mod.num.corella<-gam.8.poisson.num.corella
plot(mod.num.corella, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.corella)
qq_plot(mod.num.corella, method = 'simulate')
k.check(mod.num.corella)
summary(mod.num.corella)

gam.8.poisson.num.corella.unordered<- gam(num.corella ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.corella<-as_tibble(all_var)
fam.gam.8.num.corella <- family(mod.num.corella)
ilink.gam.8.num.corella<- fam.gam.8.num.corella$linkinv

## add the fitted values by predicting from the mod.num.corellael for the new data
ndata.8.num.corella <- add_column(ndata.8.num.corella, fit = predict(mod.num.corella, newdata = ndata.8.num.corella, type = 'response'))
ndata.8.num.corella <- bind_cols(ndata.8.num.corella, setNames(as_tibble(predict(mod.num.corella, ndata.8.num.corella, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.corella <- mutate(ndata.8.num.corella,
                               fit_resp  = ilink.gam.8.num.corella(fit_link),
                               right_upr = ilink.gam.8.num.corella(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.corella(fit_link - (2 * se_link)))

plt.num.corella.8 <- ggplot(ndata.8.num.corella, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.corella, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Corella abundance"), textstyle("(# individuals)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.corella.8
ggsave("C:Graphs August 2020//num.corella_pred.8.png")

# Fig 2 plot generation ---------------------------------------------------
fig.week.16<-wrap_plots(plt.folliculina.16,plt.membranipora.16,
          plt.mussel.16,plt.num.barn.16,plt.num.white.bryo.16,plt.num.red.bryo.16,plt.num.nudi.16,plt.num.serpulid.16,
          plt.num.corella.16, ncol=5)+
          plot_annotation(tag_levels = 'A')

theme_set(theme_classic(base_size = 5))
#theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

fig.week.16
ggplot2::ggsave(plot=fig.week.16, "C:Graphs August 2020//Fig_wk_16.pdf", width=16, height=6, units="cm")



fig.week.8<-wrap_plots(plt.folliculina.8,plt.membranipora.8,
                        plt.mussel.8,plt.num.barn.8,plt.num.white.bryo.8,plt.num.red.bryo.8,plt.num.nudi.8,plt.num.serpulid.8,
                        plt.num.corella.8, ncol=5)+
                        plot_annotation(tag_levels = 'A')

fig.week.8

ggplot2::ggsave(plot=fig.week.8, "C:Graphs August 2020//Fig_wk_8.pdf", width=16, height=6, units="cm")

head(invasion.exp.data.8_zscores)

# Pulling model 16 results to a table ----------------------------------------

botryllus.gam.16<-summary(gam.16.beta.bot.total.001)
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


botryllus.gam.16.unordered<-summary(gam.16.beta.bot.total.001.unordered)
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
ptable.16$Factor<-rep(c("Intercept", "CO2.Treatment Present"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.16 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 16) %>% 
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
stable.16$Smooth_terms<-rep(c("smooth bot.total.001", "smooth bot.total.001 * CO2.Treatment present"))

stable.16 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 16) %>% 
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
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Schizoporella", 17, 16) %>% 
  group_rows("Hermissenda", 19, 20) %>% 
  group_rows("Serpulid", 21, 22) %>% 
  group_rows("Corella", 23, 24) %>% 
  save_kable(file = "C:Biological data//pstable.16.html", self_contained = T)



# Pulling model 8  results to a table ----------------------------------------

botryllus.gam.8<-summary(gam.8.beta.bot.total.001.3)
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

botryllus.gam.8.unordered<-summary(gam.8.beta.bot.total.001.3.unordered)
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
ptable.8$Factor<-rep(c("Intercept", "CO2.Treatment Present"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.8 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 16) %>% 
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
stable.8$Smooth_terms<-rep(c("smooth bot.total.001", "smooth bot.total.001 * CO2.Treatment present"))

stable.8 %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Botryllus", 1, 2) %>%
  group_rows("Botryllus eaten", 3, 4) %>% 
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 16) %>% 
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
  group_rows("Folliculina",5, 8) %>% 
  group_rows("Membranipora", 7, 8) %>% 
  group_rows("Membranipora eaten", 9,10) %>% 
  group_rows("Mussels", 11, 12) %>% 
  group_rows("Barnacles", 13, 14) %>% 
  group_rows("Disporella", 15, 16) %>% 
  group_rows("Hermissenda", 17, 16) %>% 
  group_rows("Serpulid", 19, 20) %>% 
  group_rows("Corella", 21, 22) %>% 
  save_kable(file = "C:Biological data//pstable.8.html", self_contained = T)



# num.species.no.bot 8 ----------------------------------------------------------------
gam.8.nb.num.species.no.bot.1<- gam(num.species.no.bot ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = nb(), select=TRUE, method="REML")
gam.8.poisson.num.species.no.bot<- gam(num.species.no.bot ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.8.nb.num.species.no.bot.1, gam.8.poisson.num.species.no.bot)

mod.num.species.no.bot<-gam.8.poisson.num.species.no.bot
plot(mod.num.species.no.bot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.species.no.bot)
qq_plot(mod.num.species.no.bot, method = 'simulate')
k.check(mod.num.species.no.bot)
summary(mod.num.species.no.bot)

gam.8.poisson.num.species.no.bot.unordered<- gam(num.species.no.bot ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = poisson(), select=TRUE, method="REML")

ndata.8.num.species.no.bot<-as_tibble(all_var)
fam.gam.8.num.species.no.bot <- family(mod.num.species.no.bot)
ilink.gam.8.num.species.no.bot<- fam.gam.8.num.species.no.bot$linkinv

## add the fitted values by predicting from the mod.num.species.no.botel for the new data
ndata.8.num.species.no.bot <- add_column(ndata.8.num.species.no.bot, fit = predict(mod.num.species.no.bot, newdata = ndata.8.num.species.no.bot, type = 'response'))
ndata.8.num.species.no.bot <- bind_cols(ndata.8.num.species.no.bot, setNames(as_tibble(predict(mod.num.species.no.bot, ndata.8.num.species.no.bot, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

ndata.8.num.species.no.bot <- mutate(ndata.8.num.species.no.bot,
                               fit_resp  = ilink.gam.8.num.species.no.bot(fit_link),
                               right_upr = ilink.gam.8.num.species.no.bot(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.num.species.no.bot(fit_link - (2 * se_link)))

plt.num.species.no.bot.8 <- ggplot(ndata.8.num.species.no.bot, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.species.no.bot, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Native species richness")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.num.species.no.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.species.no.bot.8
ggsave("C:Graphs August 2020//num.species.no.bot_pred.8.png")

# Occupied space 8 ----------------------------------------------------------

gam.8.binomial.native.occupied.space<- gam(formula = cbind(native.occupied.space, 100-native.occupied.space)~ s(bot.total, k=6)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment, k=6), data = invasion.exp.data.8_zscores, family = binomial, select=TRUE, method="REML")
gam.8.beta.native.occupied.space<- gam(native.occupied.space.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial better

mod.native.occupied.space<-gam.8.binomial.native.occupied.space
plot(mod.native.occupied.space, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.native.occupied.space)
qq_plot(mod.native.occupied.space, method = 'simulate')
k.check(mod.native.occupied.space)
summary(mod.native.occupied.space)

summary(gam.8.binomial.native.occupied.space)

mod.native.occupied.space.unordered<- gam(native.occupied.space.001~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
ndata.8.native.occupied.space<-as_tibble(all_var)
fam.gam.8.native.occupied.space <- family(mod.native.occupied.space)
ilink.gam.8.native.occupied.space<- fam.gam.8.native.occupied.space$linkinv

## add the fitted values by predicting from the mod.native.occupied.spaceel for the new data
ndata.8.native.occupied.space <- add_column(ndata.8.native.occupied.space, fit = predict(mod.native.occupied.space, newdata = ndata.8.native.occupied.space, type = 'response'))
ndata.8.native.occupied.space <- bind_cols(ndata.8.native.occupied.space, setNames(as_tibble(predict(mod.native.occupied.space, ndata.8.native.occupied.space, se.fit = TRUE)[1:2]),
                                                               c('fit_link','se_link')))
ndata.8.native.occupied.space <- mutate(ndata.8.native.occupied.space,
                              fit_resp  = ilink.gam.8.native.occupied.space(fit_link),
                              right_upr = ilink.gam.8.native.occupied.space(fit_link + (2 * se_link)),
                              right_lwr = ilink.gam.8.native.occupied.space(fit_link - (2 * se_link)))

plt.native.occupied.space.8 <- ggplot(ndata.8.native.occupied.space, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = native.occupied.space.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("native.occupied.space")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.native.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.native.occupied.space.8
ggsave("C:Graphs August 2020//native.occupied.space_pred.8.png")




# CAP1 - 8 --------------------------------------------------------------------
#negative values so can't do gamma
gam.8.lm.CAP1<- gam(CAP1 ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")

mod.CAP1<-gam.8.lm.CAP1
plot(mod.CAP1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.CAP1)
qq_plot(mod.CAP1, method = 'simulate')
k.check(mod.CAP1)
summary(mod.CAP1)

gam.8.lm.CAP1.unordered<- gam(CAP1 ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")


ndata.8.CAP1<-as_tibble(all_var)
fam.gam.8.CAP1 <- family(mod.CAP1)
ilink.gam.8.CAP1<- fam.gam.8.CAP1$linkinv

## add the fitted values by predicting from the mod.CAP1el for the new data
ndata.8.CAP1 <- add_column(ndata.8.CAP1, fit = predict(mod.CAP1, newdata = ndata.8.CAP1, type = 'response'))
ndata.8.CAP1 <- bind_cols(ndata.8.CAP1, setNames(as_tibble(predict(mod.CAP1, ndata.8.CAP1, se.fit = TRUE)[1:2]),
                                                                                   c('fit_link','se_link')))
ndata.8.CAP1 <- mutate(ndata.8.CAP1,
                                        fit_resp  = ilink.gam.8.CAP1(fit_link),
                                        right_upr = ilink.gam.8.CAP1(fit_link + (2 * se_link)),
                                        right_lwr = ilink.gam.8.CAP1(fit_link - (2 * se_link)))

plt.CAP1.8 <- ggplot(ndata.8.CAP1, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = CAP1, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Partial-dbRDA axis 1\n(67% of constrained variation)")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.CAP1.8
ggsave("C:Graphs August 2020//CAP1_pred.8.png")




# Distances 8 ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.8.lm.distances<- gam(distcentroid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
gam.8.loglink.distances.1<- gam(distcentroid~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.8.gamma.distances<- gam(distcentroid~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.8.loglink.distances.1,  gam.8.lm.distances, gam.8.gamma.distances)

#gam.8.lm.distances although both are equal
mod.distcentroid<- gam.8.lm.distances
plot(mod.distcentroid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.distcentroid)
qq_plot(mod.distcentroid, method = 'simulate')
k.check(mod.distcentroid)
summary(mod.distcentroid)


summary(gam.8.lm.distances)

gam.8.lm.distances.unordered<- gam(distcentroid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")


ndata.8.distcentroid<-as_tibble(all_var)
fam.gam.8.distcentroid <- family(mod.distcentroid)
ilink.gam.8.distcentroid<- fam.gam.8.distcentroid$linkinv

## add the fitted values by predicting from the mod.distcentroidel for the new data
ndata.8.distcentroid <- add_column(ndata.8.distcentroid, fit = predict(mod.distcentroid, newdata = ndata.8.distcentroid, type = 'response'))
ndata.8.distcentroid <- bind_cols(ndata.8.distcentroid, setNames(as_tibble(predict(mod.distcentroid, ndata.8.distcentroid, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))
ndata.8.distcentroid <- mutate(ndata.8.distcentroid,
                       fit_resp  = ilink.gam.8.distcentroid(fit_link),
                       right_upr = ilink.gam.8.distcentroid(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.8.distcentroid(fit_link - (2 * se_link)))

plt.distcentroid.8 <- ggplot(ndata.8.distcentroid, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = distcentroid, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.distcentroid,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.distcentroid.8
ggsave("C:Graphs August 2020//distcentroid_pred.8.png")


# Evenness ----------------------------------------------------------------


gam.8.lm.evenness<- gam(evenness ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, select=TRUE, method="REML")
gam.8.loglink.evenness.1<- gam(evenness~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.8_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.8.loglink.evenness.1,  gam.8.lm.evenness)

mod.evenness<- gam.8.lm.evenness
plot(mod.evenness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.evenness)
qq_plot(mod.evenness, method = 'simulate')
k.check(mod.evenness)
summary(mod.evenness)

mod.evenness.unordered<- gam(evenness.001~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.8_zscores, family = betar(link="logit"), select=TRUE, method="REML")
ndata.8.evenness<-as_tibble(all_var)
fam.gam.8.evenness <- family(mod.evenness)
ilink.gam.8.evenness<- fam.gam.8.evenness$linkinv

## add the fitted values by predicting from the mod.evennessel for the new data
ndata.8.evenness <- add_column(ndata.8.evenness, fit = predict(mod.evenness, newdata = ndata.8.evenness, type = 'response'))
ndata.8.evenness <- bind_cols(ndata.8.evenness, setNames(as_tibble(predict(mod.evenness, ndata.8.evenness, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))
ndata.8.evenness <- mutate(ndata.8.evenness,
                               fit_resp  = ilink.gam.8.evenness(fit_link),
                               right_upr = ilink.gam.8.evenness(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.8.evenness(fit_link - (2 * se_link)))

plt.evenness.8 <- ggplot(ndata.8.evenness, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = evenness, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.8_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Native species evenness")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.8.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.evenness.8
ggsave("C:Graphs August 2020//evenness_pred.8.png")



###### week 16
# num.species.no.bot 16 ----------------------------------------------------------------
gam.16.nb.num.species.no.bot.1<- gam(num.species.no.bot ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = nb(), select=TRUE, method="REML")
gam.16.poisson.num.species.no.bot<- gam(num.species.no.bot ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.16.nb.num.species.no.bot.1, gam.16.poisson.num.species.no.bot)

mod.num.species.no.bot<-gam.16.poisson.num.species.no.bot
plot(mod.num.species.no.bot, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.num.species.no.bot)
qq_plot(mod.num.species.no.bot, method = 'simulate')
k.check(mod.num.species.no.bot)
summary(mod.num.species.no.bot)

gam.16.poisson.num.species.no.bot.unordered<- gam(num.species.no.bot ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = poisson(), select=TRUE, method="REML")

ndata.16.num.species.no.bot<-as_tibble(all_var)
fam.gam.16.num.species.no.bot <- family(mod.num.species.no.bot)
ilink.gam.16.num.species.no.bot<- fam.gam.16.num.species.no.bot$linkinv

## add the fitted values by predicting from the mod.num.species.no.botel for the new data
ndata.16.num.species.no.bot <- add_column(ndata.16.num.species.no.bot, fit = predict(mod.num.species.no.bot, newdata = ndata.16.num.species.no.bot, type = 'response'))
ndata.16.num.species.no.bot <- bind_cols(ndata.16.num.species.no.bot, setNames(as_tibble(predict(mod.num.species.no.bot, ndata.16.num.species.no.bot, se.fit = TRUE)[1:2]),
                                                                             c('fit_link','se_link')))

ndata.16.num.species.no.bot <- mutate(ndata.16.num.species.no.bot,
                                     fit_resp  = ilink.gam.16.num.species.no.bot(fit_link),
                                     right_upr = ilink.gam.16.num.species.no.bot(fit_link + (2 * se_link)),
                                     right_lwr = ilink.gam.16.num.species.no.bot(fit_link - (2 * se_link)))

plt.num.species.no.bot.16 <- ggplot(ndata.16.num.species.no.bot, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = num.species.no.bot, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle("Native species richness")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.num.species.no.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.num.species.no.bot.16
ggsave("C:Graphs August 2020//num.species.no.bot_pred.16.png")

# Occupied space 16 ----------------------------------------------------------

gam.16.binomial.native.occupied.space<- gam(formula = cbind(native.occupied.space, 100-native.occupied.space)~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = binomial, select=TRUE, method="REML")
gam.16.beta.native.occupied.space<- gam(native.occupied.space.001~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
#binomial better

mod.native.occupied.space<-gam.16.beta.native.occupied.space
plot(mod.native.occupied.space, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.native.occupied.space)
qq_plot(mod.native.occupied.space, method = 'simulate')
k.check(mod.native.occupied.space)
summary(mod.native.occupied.space)

summary(gam.16.binomial.native.occupied.space)

mod.native.occupied.space.unordered<- gam(native.occupied.space.001~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
ndata.16.native.occupied.space<-as_tibble(all_var)
fam.gam.16.native.occupied.space <- family(mod.native.occupied.space)
ilink.gam.16.native.occupied.space<- fam.gam.16.native.occupied.space$linkinv

## add the fitted values by predicting from the mod.native.occupied.spaceel for the new data
ndata.16.native.occupied.space <- add_column(ndata.16.native.occupied.space, fit = predict(mod.native.occupied.space, newdata = ndata.16.native.occupied.space, type = 'response'))
ndata.16.native.occupied.space <- bind_cols(ndata.16.native.occupied.space, setNames(as_tibble(predict(mod.native.occupied.space, ndata.16.native.occupied.space, se.fit = TRUE)[1:2]),
                                                                                   c('fit_link','se_link')))
ndata.16.native.occupied.space <- mutate(ndata.16.native.occupied.space,
                                        fit_resp  = ilink.gam.16.native.occupied.space(fit_link),
                                        right_upr = ilink.gam.16.native.occupied.space(fit_link + (2 * se_link)),
                                        right_lwr = ilink.gam.16.native.occupied.space(fit_link - (2 * se_link)))

plt.native.occupied.space.16 <- ggplot(ndata.16.native.occupied.space, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = native.occupied.space.001, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab(expression(atop(NA,atop(textstyle(italic("native.occupied.space")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.native.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.native.occupied.space.16
ggsave("C:Graphs August 2020//native.occupied.space_pred.16.png")

plt.native.occupied.space.16
plt.native.occupied.space.8



# CAP1 - 16 --------------------------------------------------------------------
#negative values so can't do gamma
gam.16.lm.CAP1<- gam(CAP1 ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")

mod.CAP1<-gam.16.lm.CAP1
plot(mod.CAP1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.CAP1)
qq_plot(mod.CAP1, method = 'simulate')
k.check(mod.CAP1)
summary(mod.CAP1)

gam.16.lm.CAP1.unordered<- gam(CAP1 ~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")


ndata.16.CAP1<-as_tibble(all_var)
fam.gam.16.CAP1 <- family(mod.CAP1)
ilink.gam.16.CAP1<- fam.gam.16.CAP1$linkinv

## add the fitted values by predicting from the mod.CAP1el for the new data
ndata.16.CAP1 <- add_column(ndata.16.CAP1, fit = predict(mod.CAP1, newdata = ndata.16.CAP1, type = 'response'))
ndata.16.CAP1 <- bind_cols(ndata.16.CAP1, setNames(as_tibble(predict(mod.CAP1, ndata.16.CAP1, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))
ndata.16.CAP1 <- mutate(ndata.16.CAP1,
                       fit_resp  = ilink.gam.16.CAP1(fit_link),
                       right_upr = ilink.gam.16.CAP1(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.16.CAP1(fit_link - (2 * se_link)))

plt.CAP1.16 <- ggplot(ndata.16.CAP1, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = CAP1, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Partial-dbRDA axis 1\n(67% of constrained variation)")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.CAP1.16
ggsave("C:Graphs August 2020//CAP1_pred.16.png")




# Distances 16 ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.16.lm.distances<- gam(distcentroid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.distances.1<- gam(distcentroid~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.16.gamma.distances<- gam(distcentroid~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.16.loglink.distances.1,  gam.16.lm.distances, gam.16.gamma.distances)

#gam.16.lm.distances although both are equal
mod.distcentroid<- gam.16.lm.distances
plot(mod.distcentroid, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.distcentroid)
qq_plot(mod.distcentroid, method = 'simulate')
k.check(mod.distcentroid)
summary(mod.distcentroid)


summary(gam.16.lm.distances)

gam.16.lm.distances.unordered<- gam(distcentroid ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")


ndata.16.distcentroid<-as_tibble(all_var)
fam.gam.16.distcentroid <- family(mod.distcentroid)
ilink.gam.16.distcentroid<- fam.gam.16.distcentroid$linkinv

## add the fitted values by predicting from the mod.distcentroidel for the new data
ndata.16.distcentroid <- add_column(ndata.16.distcentroid, fit = predict(mod.distcentroid, newdata = ndata.16.distcentroid, type = 'response'))
ndata.16.distcentroid <- bind_cols(ndata.16.distcentroid, setNames(as_tibble(predict(mod.distcentroid, ndata.16.distcentroid, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))
ndata.16.distcentroid <- mutate(ndata.16.distcentroid,
                               fit_resp  = ilink.gam.16.distcentroid(fit_link),
                               right_upr = ilink.gam.16.distcentroid(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.16.distcentroid(fit_link - (2 * se_link)))

plt.distcentroid.16 <- ggplot(ndata.16.distcentroid, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = distcentroid, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.distcentroid,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.distcentroid.16
ggsave("C:Graphs August 2020//distcentroid_pred.16.png")


# Evenness ----------------------------------------------------------------


gam.16.lm.evenness<- gam(evenness ~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, select=TRUE, method="REML")
gam.16.loglink.evenness.1<- gam(evenness~ s(bot.total)+ oCO2.Treatment + s(bot.total, by=oCO2.Treatment),data = invasion.exp.data.16_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.16.loglink.evenness.1,  gam.16.lm.evenness)

mod.evenness<- gam.16.lm.evenness
plot(mod.evenness, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(mod.evenness)
qq_plot(mod.evenness, method = 'simulate')
k.check(mod.evenness)
summary(mod.evenness)

mod.evenness.unordered<- gam(evenness.001~ s(bot.total)+ CO2.Treatment + s(bot.total, by=oCO2.Treatment), data = invasion.exp.data.16_zscores, family = betar(link="logit"), select=TRUE, method="REML")
ndata.16.evenness<-as_tibble(all_var)
fam.gam.16.evenness <- family(mod.evenness)
ilink.gam.16.evenness<- fam.gam.16.evenness$linkinv

## add the fitted values by predicting from the mod.evennessel for the new data
ndata.16.evenness <- add_column(ndata.16.evenness, fit = predict(mod.evenness, newdata = ndata.16.evenness, type = 'response'))
ndata.16.evenness <- bind_cols(ndata.16.evenness, setNames(as_tibble(predict(mod.evenness, ndata.16.evenness, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))
ndata.16.evenness <- mutate(ndata.16.evenness,
                           fit_resp  = ilink.gam.16.evenness(fit_link),
                           right_upr = ilink.gam.16.evenness(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.16.evenness(fit_link - (2 * se_link)))

plt.evenness.16 <- ggplot(ndata.16.evenness, aes(x = bot.total, y = fit)) + 
  geom_line(aes(colour=oCO2.Treatment)) +
  geom_jitter(aes(y = evenness, shape=Invasives, colour=oCO2.Treatment), data = invasion.exp.data.16_zscores, width=0.05, height=0.01)+
  xlab(expression("Ascidian abundance")) + ylab("Native species evenness")+  
  scale_color_manual(values=colorset_CO2.Treatment, guide = guide_legend(title="CO2.Treatment", title.position = "top"))+
  scale_fill_manual(values=colorset_CO2.Treatment, guide = FALSE)+ 
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH CO2.Treatment", title.position = "top"))+
  geom_ribbon(data = ndata.16.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oCO2.Treatment), alpha = 0.10)+
  theme(legend.position='none')
plt.evenness.16
ggsave("C:Graphs August 2020//evenness_pred.16.png")






# Community plotting ------------------------------------------------------

#### community fig week 16
fig.community.week.16<-wrap_plots(plt.native.occupied.space.16,
                            plt.num.species.no.bot.16,
                            plt.CAP1.16, plt.distances.16, ncol=2)+
                            plot_annotation(tag_levels = 'a')

fig.community.week.16

ggplot2::ggsave(plot=fig.community.week.16, "C:Graphs August 2020//Fig.community.week16.pdf", width=3, height=3, units="in")


colorset_CO2.Treatment = c("CO2"="#A20228" ,"AIR"="#818392")
#### community fig week 8
fig.community.week.8<-wrap_plots(plt.native.occupied.space.8,
                                  plt.num.species.no.bot.8,
                                  plt.CAP1.8, plt.distcentroid.8, ncol=2)+
  plot_annotation(tag_levels = 'a')

fig.community.week.8

ggplot2::ggsave(plot=fig.community.week.8, "C:Graphs August 2020//Fig.community.week8.pdf", width=4, height=3, units="in")



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
ptable.community.t.16$Factor<-rep(c("Intercept", "CO2.Treatment Present"))

ptable.community.t.16 %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (z)", 1,2) %>% 
  group_rows("Occupied space, beta (z)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
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
stable.community.f.16$Smooth_terms<-rep(c("smooth bot.total.001", "smooth bot.total.001 * CO2.Treatment Present"))

stable.community.f.16 %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
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
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
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
ptable.community.t.8$Factor<-rep(c("Intercept", "CO2.Treatment Present"))

ptable.community.t.8 %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (z)", 1,2) %>% 
  group_rows("Occupied space, beta (z)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
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
stable.community.f.8$Smooth_terms<-rep(c("smooth bot.total.001", "smooth bot.total.001 * CO2.Treatment Present"))

stable.community.f.8 %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("num.species.no.bot, poisson (Chi-square)", 1,2) %>% 
  group_rows("Occupied space, beta (Chi-square)", 3,4) %>% 
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
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
  group_rows("Partial dbRDA (1st axis), normal",5,8) %>% 
  group_rows("Heterogeneity of dispersions, normal", 7,8) %>% 
  save_kable(file = "C:Biological Data//pstable.community.8.html ", self_contained = T)




