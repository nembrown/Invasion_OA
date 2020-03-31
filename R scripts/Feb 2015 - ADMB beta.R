repeatedm.stats<-read.csv(file.choose())
#File: All week R input .csv

#If doing continuous non zero data then load the gamma file
#If doing % cover data then load the binary file
repeatedm.stats$mussel<-(repeatedm.stats$mussel+0.0001)
repeatedm.stats$botryllid<-(repeatedm.stats$botryllid+0.0001)
repeatedm.stats$bot.eaten<-(repeatedm.stats$bot.eaten+0.0001)

levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="AIR, Absent"] <- "Control, Invasives Absent"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="CO2, Absent"] <- "Elevated CO2, Invasives Absent"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="AIR, Present"] <- "Control, Invasives Present"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="CO2, Present"] <- "Elevated CO2, Invasives Present"
levels(repeatedm.stats$CO2.Treatment)[levels(repeatedm.stats$CO2.Treatment)=="AIR"] <- "Control"
levels(repeatedm.stats$CO2.Treatment)[levels(repeatedm.stats$CO2.Treatment)=="CO2"] <- "Elevated CO2"


head(repeatedm.stats)
repeatedm.stats<-repeatedm.stats[which(repeatedm.stats$Week>0), ]
repeatedm.stats$Mesocosm<-as.factor(repeatedm.stats$Mesocosm)
repeatedm.stats$Week<-as.factor(repeatedm.stats$Week)
repeatedm.stats$Tile.ID<-as.factor(repeatedm.stats$Tile.ID)
repeatedm.stats$occupied<-1-(repeatedm.stats$bare)



repeatedm.stats16<-repeatedm.stats[repeatedm.stats$Week==16, ]

repeatedm.stats14<-repeatedm.stats[repeatedm.stats$Week==14, ]

repeatedm.stats12<-repeatedm.stats[repeatedm.stats$Week==12, ]

repeatedm.stats10<-repeatedm.stats[repeatedm.stats$Week==10, ]

repeatedm.stats8<-repeatedm.stats[repeatedm.stats$Week==8, ]

repeatedm.stats6<-repeatedm.stats[repeatedm.stats$Week==6, ]

repeatedm.stats4<-repeatedm.stats[repeatedm.stats$Week==4, ]

repeatedm.stats2<-repeatedm.stats[repeatedm.stats$Week==2, ]


library(plyr)

library(ggplot2)

library(doBy)

################################################### Don't run these yet
repeatedm.stats$occupied<-100-(repeatedm.stats$bare)
repeatedm.stats$prop.bot.eaten<-(repeatedm.stats$bot.eaten)/(repeatedm.stats$botryllid + repeatedm.stats$bot.eaten)
repeatedm.stats$bot.eaten[is.na(repeatedm.stats$botryllid)]<-0

repeatedm.stats$botryllid<-(repeatedm.stats$mem.eaten)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$botryllid[is.na(repeatedm.stats$botryllid)]<-0

repeatedm.stats$prop.mem.dead<-(repeatedm.stats$mem.dead)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$prop.mem.dead[is.na(repeatedm.stats$prop.mem.dead)]<-0

repeatedm.stats$total.bot<-(repeatedm.stats$botryllid + repeatedm.stats$bot.eaten)
repeatedm.stats$mussel<-(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$mussel<-(repeatedm.stats$mussel + repeatedm.stats$mem.dead)



#############3Beta distribution:
require(car)
require(MASS)
library(betareg)
library(lmtest)
library(fitdistrplus)
beta<-fitdist((0.01*repeatedm.stats16$botryllid), "beta", start=NULL)
qqp(repeatedm.stats16$botryllid, "beta", shape1 = beta$estimate[[1]], shape2 = beta$estimate[[2]])



### Try glmm ADMB
head(repeatedm.stats16)
repeatedm.stats$Mesocosm<-as.factor(repeatedm.stats$Mesocosm)
repeatedm.stats$Tile<-as.factor(repeatedm.stats$Tile)

library(glmmADMB)
#glmmadmb(formula = NCalls ~ (FoodTreatment + ArrivalTime) * SexParent + offset(logBroodSize) + (1 | Nest), data = Owls, family = "poisson",zeroInflation = TRUE) 
# admb.opts = admbControl(shess = FALSE, noinit = FALSE)
#if you run it with this it should help

?glmmadmb

glmmadmb.botryllid<-glmmadmb(formula = (0.01*botryllid) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats16, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
plot(glmmadmb.botryllid)
summary(glmmadmb.botryllid)
###########
coefplot2(glmmadmb.botryllid)

?glmer
#########
#Model validation
#should it be pearson's residuals? plot(resid(z6, type="pearson")~fitted(z6))
#plot(x, which = c(1:3, 5), 
     caption = list("Residuals vs Fitted", "Normal Q-Q",
                    "Scale-Location", "Cook's distance",
                    "Residuals vs Leverage",
                    expression("Cook's dist vs Leverage  " * h[ii] / (1 - h[ii]))),
     panel = if(add.smooth) panel.smooth else points,
     sub.caption = NULL, main = "",
     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
     ...,
     id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
     qqline = TRUE, cook.levels = c(0.5, 1.0),
     add.smooth = getOption("add.smooth"), label.pos = c(4,2),
     cex.caption = 1)

######################

######## Need to redo

head(repeatedm.stats16)
hist(repeatedm.stats16$membranipora)

#Week 16:
beta.16<-fitdist((0.01*repeatedm.stats16$mussel), "beta", start=NULL)
qqp(repeatedm.stats16$mussel, "beta", shape1 = beta.16$estimate[[1]], shape2 = beta.16$estimate[[2]])
glmmadmb.mussel.16<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats16, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.16)
plot.glmm.admb(glmmadmb.mussel.16)
plot(resid(glmmadmb.mussel.16)~fitted(glmmadmb.mussel.16))


#Week 14:*********FALSE
beta.14<-fitdist((0.01*repeatedm.stats14$mussel), "beta", start=NULL)
qqp(repeatedm.stats14$mussel, "beta", shape1 = beta.14$estimate[[1]], shape2 = beta.14$estimate[[2]])
glmmadmb.mussel.14<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats14, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.14)


#Week 12:
beta.12<-fitdist((0.01*repeatedm.stats16$mussel), "beta", start=NULL)
qqp(repeatedm.stats12$mussel, "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
glmmadmb.mussel.12<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats12, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.12)


#Week 10:
beta.10<-fitdist((0.01*repeatedm.stats10$mussel), "beta", start=NULL)
qqp(repeatedm.stats10$mussel, "beta", shape1 = beta.10$estimate[[1]], shape2 = beta.10$estimate[[2]])
glmmadmb.mussel.10<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats10, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.10)


#Week 8:
beta.8<-fitdist((0.01*repeatedm.stats8$mussel), "beta", start=NULL)
qqp(repeatedm.stats8$mussel, "beta", shape1 = beta.8$estimate[[1]], shape2 = beta.8$estimate[[2]])
glmmadmb.mussel.8<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats8, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.8)


#Week 6: ***FALSE
beta.6<-fitdist((0.01*repeatedm.stats6$mussel), "beta", start=NULL)
qqp(repeatedm.stats6$mussel, "beta", shape1 = beta.6$estimate[[1]], shape2 = beta.6$estimate[[2]])
glmmadmb.mussel.6<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats6, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE, maxfn=10000)) 
summary(glmmadmb.mussel.6)


## try binomial: 
glmmadmb.mussel.6.binom<-glmmadmb(formula = cbind(mussel, (100-mussel)) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats6, family = "binomial") 
summary(glmmadmb.mussel.6.binom)



#Week 4:****** NOT WORKING 
fix(repeatedm.stats4)

beta.4<-fitdist((0.01*repeatedm.stats4$mussel), "beta", start=NULL)
qqp(repeatedm.stats4$mussel, "beta", shape1 = beta.4$estimate[[1]], shape2 = beta.4$estimate[[2]])
glmmadmb.mussel.4<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment * Invasives + (1 | Mesocosm), data = repeatedm.stats4, family = "beta",zeroInflation = TRUE, admb.opts=admbControl(shess = FALSE, noinit = FALSE))
summary(glmmadmb.mussel.4)

#Week 2:
beta.2<-fitdist((0.01*repeatedm.stats2$mussel), "beta", start=NULL)
qqp(repeatedm.stats2$mussel, "beta", shape1 = beta.2$estimate[[1]], shape2 = beta.2$estimate[[2]])
glmmadmb.mussel.2<-glmmadmb(formula = (0.01*mussel) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats2, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.2)


#One sign of overdispersion is that the residual deviance is significantly higher than the residual degrees of freedom

####################
glmmadmb.mussel.6.1<-glmmadmb(formula = cbind((0.01*mussel), (1-(0.01*mussel))) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats6, family = "binomial",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.mussel.6.1)
###############33

repos <- getOption('repos')
repos["CRAN"] <- "http://cran.rstudio.org"
option(repos = repos)
install.packages('UsingR')

install.packages("coda")
install.packages("reshape")
install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")
require("coefplot2")

?formula
#One good way to assess the results of a model fit is to look at a coefficient plot:
#just the results ... can compare to other models!! 
### coefplot2(list(ZIP=fit_zipoiss,ZINB=fit_zinbinom,ZINB1=fit_zinbinom1,ZINB1_brood=fit_zinbinom1_bs),varnames=vn,legend=TRUE)
## see the ben bolker notes and admb files

coefplot2(glmmadmb.botryllid)
dotchart(glmmadmb.botryllid)
?glmmadmb

### other diagnostic plots
qqnorm(residuals(glmmadmb.botryllid));qqline(residuals(glmmadmb.botryllid), col="red")
plot(resid(glmmadmb.botryllid)~fitted(glmmadmb.botryllid))
plot(resid(glmmadmb.botryllid)~repeatedm.stats16$Mesocosm)
plot(resid(glmmadmb.botryllid)~repeatedm.stats16$botryllid)### residuals vs. observed
#residuals vs observed will show correlation. They're more difficult to interpret because of this. 
#Residuals vs fitted shows the best approximation we have to how the errors relate to the population mean, 
#and is somewhat useful for examining the more usual consideration in regression of whether variance is related to mean.

plot(residuals(glmmadmb.mussel.16["within"])~fitted(glmmadmb.mussel.16["within"]))


################
glmmadmb.botryllid2<-glmmadmb(formula = (0.01*botryllid) ~ CO2.Treatment*Invasives, random=~1 | Mesocosm, data = repeatedm.stats16, family = "beta",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.botryllid2)
## exactly the same as above formula

###########################################################################
###Model validation - Graphing distributions
 
#Residuals still have problems, because of week. Is this the right function to plot residuals? 
qqnorm(residuals(glmmadmb.botryllid));qqline(residuals(glmmadmb.botryllid), col="red")
plot(resid(glmmadmb.botryllid)~fitted(glmmadmb.botryllid))
#Makes a pattern

library(boot)
glm.diag.plots(glmmadmb.mussel.16, glm.diag(glmmadmb.mussel.16))
# can't use for glmm, only glm I think









########################## BINOMIAL .... 


################ Glmm admb for binomial model, should use this instead of pql anyways. 
windows()
p4<- ggplot(repeatedm.stats16.Invasives.Only, aes(x=repeatedm.stats16.Invasives.Only$prop.bot.eaten, fill=repeatedm.stats16.Invasives.Only$CO2.Treatment)) + geom_density(alpha=.3) + scale_fill_brewer(palette="Set1")

p4<-p4 + theme_bw() + labs(x="prop.bot.eaten ", y="Kernel Density")+ theme(text = element_text(size=20), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

p4 + theme(legend.text = element_text(colour="black", size = 20))+ theme(legend.title = element_text(colour="white", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

### problem... 

#If you want to fit a binomial model, you need to specify the
#_number_ surviving and the total (denominator), as in

#glmmadmb(cbind(nSurv,nTotal) ~ ....)
#cbind(nSurv,nTotal)
#### don't think it should be total .... shoul dbe just botryllid???
head(repeatedm.stats16)

repeatedm.stats16.Invasives.Only<-repeatedm.stats16[repeatedm.stats16$Invasives=="Present", ]
repeatedm.stats14.Invasives.Only<-repeatedm.stats14[repeatedm.stats14$Invasives=="Present", ]
repeatedm.stats12.Invasives.Only<-repeatedm.stats12[repeatedm.stats12$Invasives=="Present", ]
repeatedm.stats10.Invasives.Only<-repeatedm.stats10[repeatedm.stats10$Invasives=="Present", ]
repeatedm.stats8.Invasives.Only<-repeatedm.stats8[repeatedm.stats8$Invasives=="Present", ]
repeatedm.stats6.Invasives.Only<-repeatedm.stats6[repeatedm.stats6$Invasives=="Present", ]
repeatedm.stats4.Invasives.Only<-repeatedm.stats4[repeatedm.stats4$Invasives=="Present", ]
repeatedm.stats2.Invasives.Only<-repeatedm.stats2[repeatedm.stats2$Invasives=="Present", ]


glmmadmb.prop.bot.eaten.16<-glmmadmb(formula = cbind(bot.eaten, botryllid) ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats16, family = "binomial",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.bot.eaten.16)

#Don't really need glmm if just using invasives present, i.e. one tile from each mesocosm. 

?glm
#do I need zero inflation though? 
###One sign of overdispersion is that the residual deviance is significantly higher than the residual degrees of freedom

hist(repeatedm.stats16.Invasives.Only$prop.bot.eaten)

#Week 16:
glm.prop.bot.eaten.16<-glm(formula = cbind(bot.eaten, botryllid) ~ CO2.Treatment, data = repeatedm.stats16.Invasives.Only, family = "binomial") 
summary(glm.prop.bot.eaten.16)

?glm


1-pchisq(409.67,22)

glm.prop.bot.eaten.16.<-glm(formula = (0.01*prop.bot.eaten) ~ CO2.Treatment, data = repeatedm.stats16.Invasives.Only, family = "binomial") 
summary(glm.prop.bot.eaten.16.)
## it DOES matter if you do the prop or put in the raw data

#Week 16:
glm.prop.bot.eaten.16.binom<-glm(formula = cbind(bot.eaten, botryllid) ~ CO2.Treatment, data = repeatedm.stats16.Invasives.Only, family = "binomial") 
summary(glm.prop.bot.eaten.16.binom)

anova(glm.prop.bot.eaten.16.binom, glm.prop.bot.eaten.16)

glm.prop.bot.eaten.16<-glm(formula = cbind(bot.eaten, botryllid) ~ CO2.Treatment, data = repeatedm.stats16.Invasives.Only, family = "quasibinomial")
plot(glm.prop.bot.eaten.16)
summary(glm.prop.bot.eaten.16)
#overdispersed, wayyy better with quasibinomial, plots 
#AIC 469 vs

#Week 14:
glm.prop.bot.eaten.14<-glm(formula = cbind((bot.eaten), (botryllid)) ~ CO2.Treatment, data = repeatedm.stats14.Invasives.Only, family = "quasibinomial") 
plot(glm.prop.bot.eaten.14)
summary(glm.prop.bot.eaten.14)
1-pchisq(295.547,22)

#Week 12:
glm.prop.bot.eaten.12<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats12.Invasives.Only, family = "quasibinomial") 
plot(glm.prop.bot.eaten.12)
summary(glm.prop.bot.eaten.12)


#Week 10:
glm.prop.bot.eaten.10<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats10.Invasives.Only, family = "quasibinomial") 
summary(glm.prop.bot.eaten.10)
plot(glm.prop.bot.eaten.10)


#Week 8:
glm.prop.bot.eaten.8<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats8.Invasives.Only, family = "quasibinomial") 
summary(glm.prop.bot.eaten.8)
plot(glm.prop.bot.eaten.8)

#Week 6:
glm.prop.bot.eaten.6<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats6.Invasives.Only, family = "quasibinomial") 
summary(glm.prop.bot.eaten.6)
plot(glm.prop.bot.eaten.6)

#Week 4:
glm.prop.bot.eaten.4<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats4.Invasives.Only, family = "binomial") 
summary(glm.prop.bot.eaten.4)
plot(glm.prop.bot.eaten.4)

#Week 2:
glm.prop.bot.eaten.2<-glm(formula = cbind(bot.eaten, (botryllid)) ~ CO2.Treatment  , data = repeatedm.stats2.Invasives.Only, family = "binomial") 
summary(glm.prop.bot.eaten.2)
plot(glm.prop.bot.eaten.2)





############################ Count data poissson - laplace in admb
repeatedm.counts<-read.csv(file.choose())
#File: All week R input .csv
head(repeatedm.counts)
#If doing continuous non zero data then load the gamma file
#If doing % cover data then load the binary file


levels(repeatedm.counts$Treatment)[levels(repeatedm.counts$Treatment)=="AIR, Absent"] <- "Control, Invasives Absent"
levels(repeatedm.counts$Treatment)[levels(repeatedm.counts$Treatment)=="CO2, Absent"] <- "Elevated CO2, Invasives Absent"
levels(repeatedm.counts$Treatment)[levels(repeatedm.counts$Treatment)=="AIR, Present"] <- "Control, Invasives Present"
levels(repeatedm.counts$Treatment)[levels(repeatedm.counts$Treatment)=="CO2, Present"] <- "Elevated CO2, Invasives Present"
levels(repeatedm.counts$CO2.Treatment)[levels(repeatedm.counts$CO2.Treatment)=="AIR"] <- "Control"
levels(repeatedm.counts$CO2.Treatment)[levels(repeatedm.counts$CO2.Treatment)=="CO2"] <- "Elevated CO2"


head(repeatedm.counts)
repeatedm.counts<-repeatedm.counts[which(repeatedm.counts$Week>0), ]
repeatedm.counts$Mesocosm<-as.factor(repeatedm.counts$Mesocosm)
repeatedm.counts$Week<-as.factor(repeatedm.counts$Week)
repeatedm.counts$Tile.ID<-as.factor(repeatedm.counts$Tile.ID)
repeatedm.counts$occupied<-1-(repeatedm.counts$bare)



repeatedm.counts16<-repeatedm.counts[repeatedm.counts$Week==16, ]

repeatedm.counts14<-repeatedm.counts[repeatedm.counts$Week==14, ]

repeatedm.counts12<-repeatedm.counts[repeatedm.counts$Week==12, ]

repeatedm.counts10<-repeatedm.counts[repeatedm.counts$Week==10, ]

repeatedm.counts8<-repeatedm.counts[repeatedm.counts$Week==8, ]

repeatedm.counts6<-repeatedm.counts[repeatedm.counts$Week==6, ]

repeatedm.counts4<-repeatedm.counts[repeatedm.counts$Week==4, ]

repeatedm.counts2<-repeatedm.counts[repeatedm.counts$Week==2, ]

?glmmadmb

#Week 16:
poisson16<- fitdistr(repeatedm.counts16$serpulid, "Poisson")
qqp(repeatedm.counts16$serpulid, "pois", poisson16$estimate)
glmmadmb.serpulid.16<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts16, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.16)

#Week 14:
poisson14<- fitdistr(repeatedm.counts14$serpulid, "Poisson")
qqp(repeatedm.counts14$serpulid, "pois", poisson14$estimate)
glmmadmb.serpulid.14<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts14, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.14)


#Week 12:
poisson12<- fitdistr(repeatedm.counts12$serpulid, "Poisson")
qqp(repeatedm.counts12$serpulid, "pois", poisson12$estimate)
glmmadmb.serpulid.12<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts12, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.12)

?glmmadmb
#Week 10: ######3 should maybe try neg binomial??
poisson10<- fitdistr(repeatedm.counts10$serpulid, "Poisson")
qqp(repeatedm.counts10$serpulid, "pois", poisson10$estimate)

nbinom10 <- fitdistr(repeatedm.counts10$serpulid, "Negative Binomial")
qqp(repeatedm.counts10$serpulid, "nbinom", size = nbinom10$estimate[[1]], mu = nbinom10$estimate[[2]])


glmmadmb.serpulid.10<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts10, family = "nbinom1",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.10)


#Week 8:
poisson8<- fitdistr(repeatedm.counts8$serpulid, "Poisson")
qqp(repeatedm.counts8$serpulid, "pois", poisson8$estimate)

nbinom8 <- fitdistr(repeatedm.counts8$serpulid, "Negative Binomial")
qqp(repeatedm.counts8$serpulid, "nbinom", size = nbinom8$estimate[[1]], mu = nbinom8$estimate[[2]])

glmmadmb.serpulid.8<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts8, family = "nbinom1",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.8)


#Week 6:
poisson6<- fitdistr(repeatedm.counts6$serpulid, "Poisson")
qqp(repeatedm.counts6$serpulid, "pois", poisson6$estimate)
nbinom6 <- fitdistr(repeatedm.counts6$serpulid, "Negative Binomial")
qqp(repeatedm.counts6$serpulid, "nbinom", size = nbinom6$estimate[[1]], mu = nbinom6$estimate[[2]])

glmmadmb.serpulid.6<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts6, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.6)

#Week 4:
poisson4<- fitdistr(repeatedm.counts4$serpulid, "Poisson")
qqp(repeatedm.counts4$serpulid, "pois", poisson4$estimate)
glmmadmb.serpulid.4<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts4, family = "nbinom1",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.4)

#Week 2:
poisson2<- fitdistr(repeatedm.counts2$serpulid, "Poisson")
qqp(repeatedm.counts2$serpulid, "pois", poisson2$estimate)
glmmadmb.serpulid.2<-glmmadmb(formula = serpulid ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.counts2, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.serpulid.2)



library(glmmADMB)
#################### Native species richness ... just counts # of species therefore poissson
#Week 16:
poisson16<- fitdistr(repeatedm.stats16$num.species.no.bot, "Poisson")
qqp(repeatedm.stats16$num.species.no.bot, "pois", poisson16$estimate)
nbinom16 <- fitdistr(repeatedm.stats16$num.species.no.bot, "Negative Binomial")
qqp(repeatedm.stats16$num.species.no.bot, "nbinom", size = nbinom16$estimate[[1]], mu = nbinom16$estimate[[2]])


glmmadmb.num.species.no.bot.16<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats16, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.16)

#Week 14:
poisson14<- fitdistr(repeatedm.stats14$num.species.no.bot, "Poisson")
qqp(repeatedm.stats14$num.species.no.bot, "pois", poisson14$estimate)
glmmadmb.num.species.no.bot.14<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats14, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.14)


#Week 12:
poisson12<- fitdistr(repeatedm.stats12$num.species.no.bot, "Poisson")
qqp(repeatedm.stats12$num.species.no.bot, "pois", poisson12$estimate)

nbinom12 <- fitdistr(repeatedm.stats12$num.species.no.bot, "Negative Binomial")
qqp(repeatedm.stats12$num.species.no.bot, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])


glmmadmb.num.species.no.bot.12<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats12, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.12)

?glmmadmb
#Week 10: ######3 should maybe try neg binomial??
poisson10<- fitdistr(repeatedm.stats10$num.species.no.bot, "Poisson")
qqp(repeatedm.stats10$num.species.no.bot, "pois", poisson10$estimate)

nbinom10 <- fitdistr(repeatedm.stats10$num.species.no.bot, "Negative Binomial")
qqp(repeatedm.stats10$num.species.no.bot, "nbinom", size = nbinom10$estimate[[1]], mu = nbinom10$estimate[[2]])


glmmadmb.num.species.no.bot.10<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats10, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.10)


#Week 8:
poisson8<- fitdistr(repeatedm.stats8$num.species.no.bot, "Poisson")
qqp(repeatedm.stats8$num.species.no.bot, "pois", poisson8$estimate)

nbinom8 <- fitdistr(repeatedm.stats8$num.species.no.bot, "Negative Binomial")
qqp(repeatedm.stats8$num.species.no.bot, "nbinom", size = nbinom8$estimate[[1]], mu = nbinom8$estimate[[2]])

glmmadmb.num.species.no.bot.8<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats8, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.8)


#Week 6:
poisson6<- fitdistr(repeatedm.stats6$num.species.no.bot, "Poisson")
qqp(repeatedm.stats6$num.species.no.bot, "pois", poisson6$estimate)
nbinom6 <- fitdistr(repeatedm.stats6$num.species.no.bot, "Negative Binomial")
qqp(repeatedm.stats6$num.species.no.bot, "nbinom", size = nbinom6$estimate[[1]], mu = nbinom6$estimate[[2]])

glmmadmb.num.species.no.bot.6<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats6, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.6)

#Week 4:
poisson4<- fitdistr(repeatedm.stats4$num.species.no.bot, "Poisson")
qqp(repeatedm.stats4$num.species.no.bot, "pois", poisson4$estimate)
glmmadmb.num.species.no.bot.4<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats4, family = "nbinom1",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.4)

#Week 2:
poisson2<- fitdistr(repeatedm.stats2$num.species.no.bot, "Poisson")
qqp(repeatedm.stats2$num.species.no.bot, "pois", poisson2$estimate)
glmmadmb.num.species.no.bot.2<-glmmadmb(formula = num.species.no.bot ~ CO2.Treatment*Invasives + (1 | Mesocosm), data = repeatedm.stats2, family = "poisson",zeroInflation = TRUE, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.num.species.no.bot.2)






################### Membranipora eaten
######################################
#Week 16:
head(repeatedm.stats16)

glmmadmb.prop.mem.eaten.16.binom<-glmmadmb(formula = cbind(mem.eaten, membranipora) ~  CO2.Treatment*Invasives, data = repeatedm.stats16 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.16.binom)


#Week 14:
glmmadmb.prop.mem.eaten.14<-glmmadmb(formula = cbind((mem.eaten), (membranipora)) ~  CO2.Treatment*Invasives, data = repeatedm.stats14 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.14)

#Week 12:
glmmadmb.prop.mem.eaten.12<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats12 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.12)


#Week 10:
glmmadmb.prop.mem.eaten.10<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats10 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.10)



#Week 8:
glmmadmb.prop.mem.eaten.8<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats8 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.8)
plot(glmmadmb.prop.mem.eaten.8)

#Week 6:
glmmadmb.prop.mem.eaten.6<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats6 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.6)
plot(glmmadmb.prop.mem.eaten.6)

#Week 4:
glmmadmb.prop.mem.eaten.4<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats4 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.4)
plot(glmmadmb.prop.mem.eaten.4)

#Week 2:
glmmadmb.prop.mem.eaten.2<-glmmadmb(formula = cbind(mem.eaten, (membranipora)) ~  CO2.Treatment*Invasives  , data = repeatedm.stats2 , family = "binomial",zeroInflation = T, admb.opts = admbControl(shess = FALSE, noinit = FALSE)) 
summary(glmmadmb.prop.mem.eaten.2)
plot(glmmadmb.prop.mem.eaten.2)

