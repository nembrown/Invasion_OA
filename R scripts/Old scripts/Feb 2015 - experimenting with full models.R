repeatedm.stats<-read.csv(file.choose())
#File: All week R input .csv

#If doing continuous non zero data then load the gamma file
#If doing % cover data then load the binary file

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


repeatedm.stats$occupied<-100-(repeatedm.stats$bare)
repeatedm.stats$botryllid<-(repeatedm.stats$bot.eaten)/(repeatedm.stats$botryllid + repeatedm.stats$bot.eaten)
repeatedm.stats$botryllid[is.na(repeatedm.stats$botryllid)]<-0

repeatedm.stats$botryllid<-(repeatedm.stats$mem.eaten)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$botryllid[is.na(repeatedm.stats$botryllid)]<-0

repeatedm.stats$prop.mem.dead<-(repeatedm.stats$mem.dead)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$prop.mem.dead[is.na(repeatedm.stats$prop.mem.dead)]<-0

repeatedm.stats$total.bot<-(repeatedm.stats$botryllid + repeatedm.stats$bot.eaten)
repeatedm.stats$mussel<-(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$mussel<-(repeatedm.stats$mussel + repeatedm.stats$mem.dead)

####### LME 
library(nlme)
library(visreg)
library(lme4)

##########################################
# Split plot:
#lme: lme(y~A*B, random = 1|plot) or if I don't have weeks here then, : lme(y~Co2*Invasives, random=1|mesocosm)
#This is actually correct and will work for a single week. Haven't actually included the idea of all the weeks
#If I don't specify a week, it will average over all weeks, which isn't good. 

zlme1<-lme(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats, method="REML")
plot(zlme1)
summary(zlme1)

## add in varIdent sturcture ... can't do with Tile. ID for some reason ... ??? 
vf1<-varIdent(form=~1|Week)
zlme2<-lme(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm, data=repeatedm.stats, weights=vf1)
plot(zlme2)
summary(zlme2)
#residuals look better
#Becasue averaging over all week snot sure what the parameter values are telling me - intercept 1.55 
### this might work overall in the long run ... come back to this later 



#If I use repeatedm.stats (i.e. with 0), then the summary produces much less text --> is it because I specified as factor earlier? 

head(repeatedm.stats)
z6<-lme(botryllid  ~ pH.uptowk*Invasives*Week, random=~1|Mesocosm, data=repeatedm.stats)
plot(z6)
summary(z6)

### I think as soon as I treat week as a factor then the residuals are in lines.... don't look good/// form lines. 
### Can I make this a glmm??
# YEs. Either GLMM full model or GLMM with just week ... remember still need the residuals to look okay with that 
#model distribution (e.g. binomial), check if the model fits the data. 

#going to work with glmmpql b/c that's what zuur uses and is suggested in the 2008 Bolker TREE paper ...they also talk about lmer from the lme4 package and
# glmmML from the package of the same name


library(MASS)
botryllid.PQL<- glmmPQL(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm, family=poisson, data=repeatedm.stats)
summary(botryllid.PQL)
plot(botryllid.PQL)
#Residuals still have problems, because of week. Is this the right function to plot residuals? 
qqnorm(residuals(botryllid.PQL));qqline(residuals(botryllid.PQL), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(botryllid.PQL)~repeatedm.stats$Week)


###try negative binomial
# run first non mixed to get theta value:
botryllid.glm.nb<- glm.nb(botryllid ~ CO2.Treatment*Invasives*Week, data=repeatedm.stats)
summary(botryllid.glm.nb)

botryllid.PQL.nb2<- glmmPQL(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile, negative.binomial(theta =227040 , link = log), data=repeatedm.stats)
plot(botryllid.PQL.nb2)
summary(botryllid.PQL.nb2)

###### Looking up repeated measures online -> need to include week as random: 
###(Day | Subject),
botryllid.PQL<- glmmPQL(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Week/Mesocosm, family=poisson, data=repeatedm.stats)
plot(botryllid.PQL)
summary(botryllid.PQL)


repeatedm.stats$Week<-as.numeric(repeatedm.stats$Week)

#### Not sure how this works...
##Can also include Tile (i.e. individual/subject as a random factor ******** I think this is important)
botryllid.PQL2<- glmmPQL(botryllid ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile.ID, family=poisson, data=repeatedm.stats)
plot(botryllid.PQL2)
summary(botryllid.PQL2)

#Model validation: (from https://danieljhocking.wordpress.com/tag/glmmpql/)
# Check for residual pattern within groups and difference between groups     

#might be useful to plot by random factor
library(lattice)
xyplot(residuals(botryllid.PQL) ~ fitted(botryllid.PQL) | repeatedm.stats$Mesocosm)
xyplot(residuals(botryllid.PQL) ~ fitted(botryllid.PQL)) 
#residuals in general look bad. 
plot.new
plot(fitted(botryllid.PQL), residuals(botryllid.PQL),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(botryllid.PQL), residuals(botryllid.PQL)))




  #Without week, residuals look okay     
PQL.nowk<- glmmPQL(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm, family=poisson, data=repeatedm.stats)
plot(PQL.nowk)
plot(resid(PQL.nowk)~repeatedm.stats$Week)





head(repeatedm.stats12)
## Try lmer instead of glmmPQL - laplace approximation, wald's Z
library(lme4)
glmer1<- glmer((0.01*botryllid) ~ CO2.Treatment*Invasives + (1|Mesocosm), family=binomial, data=repeatedm.stats12)
summary(glmer1)
plot(lmer1)
qqnorm(residuals(lmer1));qqline(residuals(lmer1), col="red")

#When I include week not sure if that's okay? 
#What's the best for repeated measures? 
#What is the best way to evaluate the model, model assumptions etc ... how do I know if it's a good model for my data? 

lmer2<- glmer(botryllid ~ CO2.Treatment*Invasives + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer2)

anova(lmer1, lmer2, test = "Chisq")
#So obviously week is important

lmer3<- glmer(botryllid ~ CO2.Treatment*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer3)
anova(lmer1, lmer3, test = "Chisq")
# so obviously invasives are important

lmer4<- glmer(botryllid ~ Invasives*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer4)
plot(lmer4)
anova(lmer1, lmer4, test = "Chisq")
# so CO2 not important 

lmer5<- glmer(botryllid ~ CO2.Treatment:Invasives  + Invasives*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer5)
anova(lmer1, lmer5, test = "Chisq")

#am concerned about dropping terms with nesting.




## Try glmmML instead of glmmPQL
library(glmmML)
glmmML1<-glmmML(botryllid ~ CO2.Treatment*Invasives*Week, cluster=Mesocosm, family=poisson, data=repeatedm.stats)
summary(glmmML1)
#this is odd becasue all significant now??? 

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




#One sign of overdispersion is that the residual deviance is significantly higher than the residual degrees of freedom



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

plot(residuals(glmmadmb.botryllid["within"])~fitted(glmmadmb.botryllid["within"]))


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
glm.diag.plots(botryllid.PQL2, glm.diag(botryllid.PQL2))
# can't use for glmm, only glm I think



#From https://danieljhocking.wordpress.com/2011/07/18/model-validation-interpreting-residual-plots/#comments
#Took me a while to figure this out too. I suggest you generate negative binomial quantiles: nbquant<-rnbinom(n=length(count), scale=1.480, mu=mean(count))
#Then plot those quantiles against your residuals:qqnorm(nbquant, resid(glmmPQLnb1)I think you'll find that your points fall along a straight line
#, suggesting that you did choose the right distribution. I had this issue too, took me a while to figure out, and once I did, I got a nice 
#straight line.


nbquant<-rnbinom(n=length(repeatedm.stats$botryllid), mu=mean(repeatedm.stats$botryllid), size=5)
### not sure what size needs to be - it says number of succcesses?? 
qqplot(nbquant, resid(botryllid.PQL2))
#what is this plotting??? SHould it be a straight line if perfect fit?? 

#fitdirst ## avoid spurious accuracy

x <- rpois(24, 1)
fitdistr(x, "Poisson")
## now do this directly with more control.
fitdistr(x, dpois, list(shape = 1, rate = 0.1), lower = 0.001)
#### Don't understand what's going on here. 

#################ADDED FEB 2015
## Distributions
### Gamma probably shouldn't be used for % cover ... 


library(MASS)



### another technique from: http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
require(car)

require(MASS)
gamma <- fitdistr(repeatedm.stats16$botryllid, "Gamma")
qqp(repeatedm.stats16$botryllid, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#####"Check out the plots I've generated using qqp. The y axis represents the observations and the x axis represents 
###the quantiles modeled by the distribution. The solid red line represents a perfect distribution fit and the dashed 
#red lines are the confidence intervals of the perfect distribution fit. You want to pick the distribution for which 
#the largest number of observations falls between the dashed lines. In this case, that's the lognormal distribution, 
#in which only one observation falls outside the dashed lines. Now, armed with the knowledge of which probability 
#distribution fits best, I can try fitting a model."

### for fun try neg.binomial:
nbinom <- fitdistr(repeatedm.stats16$botryllid, "Negative Binomial")
qqp(repeatedm.stats16$botryllid, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])


#binomial??

#poisson
poisson <- fitdistr(repeatedm.stats16$botryllid, "Poisson")
qqp(repeatedm.stats16$botryllid, "pois", poisson$estimate)

#log normal:
qqp(repeatedm.stats10$botryllid, "lnorm")


?dbeta
?fitdistr


### Binomial distrbutoin is not the appropriate distribution, should use beta - according to this thread: 
#http://stats.stackexchange.com/questions/99425/distribution-for-percentage-data
#http://stats.stackexchange.com/questions/94539/remove-effect-of-a-factor-on-continuous-proportion-data-using-regression-in-r/94766#94766
library(betareg)
library(lmtest)
library(fitdistrplus)
beta<-fitdist((0.01*repeatedm.stats16$botryllid), "beta", start=NULL)
qqp(repeatedm.stats16$botryllid, "beta", shape1 = beta$estimate[[1]], shape2 = beta$estimate[[2]])

dx <- seq(0, 1, len = 200) 
plot(dx, dbeta(dx, shape1 = beta$estimate[[1]], shape2 = beta$estimate[[2]]), type = "l", col="red", ylim=c(0, 15)) 
d.AA <- density(0.01*repeatedm.stats16$botryllid)
lines(d.AA)





#######Beta doesn't seemto work with fitdist when I have lots of zeros --> see if it works forthe model
library(betareg)
?betareg
head(repeatedm.stats16)

botryllid.betareg<- betareg((0.01*botryllid) ~ pH, data = repeatedm.stats16)
plot(botryllid.betareg)
summary(botryllid.betareg)

plot(botryllid.betareg, which = 1:4, type = "pearson")
plot(botryllid.betareg, which = 5, type = "deviance", sub.caption="")
plot(botryllid.betareg, which = 1, type = "deviance", sub.caption="")



####################################################################################################################



##### TRY JUST ONE WEEK at a time: 
library(nlme)
head(repeatedm.stats)
head(repeatedm.stats12)
z4<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm.stats12)
plot(z4)
summary(z4)
EP<-resid(z4, type="pearson")
qqnorm(EP);qqline(EP, col="red")

plot(EP ~ fitted(z4))


#########################
EP<-resid(z, type="pearson")
ED<-resid(Muss.glm, type="deviance")
mu<-predict(Muss.glm, type="response")
E<-mydata$Number - mu
plot(x=mu, y=E, main = "Response residuals")
plot(x=mu, y=EP, main = "Pearson residuals")
plot(x=mu, y=ED, main = "Deviance residuals")
########################

#Residuals look bad so use a glm??? but glm doesn't really predict the values well at all 


z5<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats12, weights=varPower())
plot(z5)
summary(z5)
## included varPower - still bad ... try var ident: 
#vf<-varIdent(form=~1|Treatment)

vf<-varIdent(form=~1|CO2.Treatment*Invasives)
z6<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm.stats12, weights=vf)
plot(z6)
summary(z6)
## residuals look a lot better... is this okay though?? SAying that it's okay for variance to vary across treatment
plot(resid(z6, type="pearson")~fitted(z6))
## making sure its pearson that's plotted. 


### I think with the varIDent structure is the best so far... at least it's predicting the variables correctly and withvarIdent is predicting correctly
###############help(gnls, pack=nlme)

#Using Var combined
#varComb(varIdent(form=~1|Others.Type), varIdent(form=~1|Game),     
#varIdent(form=~1|Type), varIdent(form=~1|Male))

vf2<-varComb(varIdent(form=~1|CO2.Treatment), varIdent(form=~1|CO2.Treatment*Invasives))
z7<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats12, weights=vf2)
plot(z7)
summary(z7)

##check AIC to see which is better model
anova(z6,z7)
anova(z5,z6)
anova(z6,z4)
#### ns for the last one, so z6 is the best model I think. ********** Need to figure out if this is okay?? 


### Plot it out
ggplot(repeatedm.stats12, aes(x=Treatment, y=fitted(z4))) +geom_boxplot()+ scale_colour_brewer(palette="Set1")

ggplot(repeatedm.stats12, aes(x=Treatment, y=botryllid)) +geom_boxplot()+ scale_colour_brewer(palette="Set1")
#rename my plot to p, makes it easier to work with

#################################33
#Transformations
z8<-lme(log(botryllid+1) ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats12)
plot(z8)
summary(z8)
### doesn't look good... and remember: don't log transform count data, Bitches!!!! (Jared's blog post)




############################33

## GLMMS: 

############Poisson/ neg binomial

z.PQL<- glmmPQL((0.01*botryllid) ~ CO2.Treatment*Invasives, random=~1|Mesocosm, family=binomial, data=repeatedm.stats12)
plot(z.PQL)
qqnorm(residuals(z.PQL));qqline(residuals(z.PQL), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(z.PQL)~fitted(z.PQL))
#residuals look a lot better 
summary(z.PQL)

head(repeatedm.stats12)






## TRy negative binomial ... not better ... why is theta so huge? The plot function works better here... I think poisson is okay... 
z.nb<- glm.nb(botryllid ~ CO2.Treatment*Invasives,data=repeatedm.stats12)
summary(z.nb)
plot(z.nb)

z.nb2<- glmmPQL(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, negative.binomial(theta =21218   , link = log), data=repeatedm.stats12)
plot(z.nb2)
summary(z.nb2)



###### Problem that I'm getting each time is that the residuals vs fitted values are heteroskedastic.... 
###estimates seem really off ... not even close - say it's 1.76 ... but really like 5 ?? 

ggplot(repeatedm.stats12, aes(x=Treatment, y=fitted(z.nb2))) +geom_boxplot()+ scale_colour_brewer(palette="Set1")

ggplot(repeatedm.stats12, aes(x=Treatment, y=botryllid)) +geom_boxplot()+ scale_colour_brewer(palette="Set1")





############## GAMMA
head(repeatedm.stats12)
z.gamma<- glmmPQL(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, family=Gamma, data=repeatedm.stats12)
plot(z.gamma)
qqnorm(residuals(z.gamma));qqline(residuals(z.gamma), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(z.gamma)~fitted(z.gamma))
#residuals look a lot better 
summary(z.gamma)


############ BINOMIAL
head(repeatedm.stats12)
z.binomial<- glmmPQL((0.01*botryllid) ~ CO2.Treatment*Invasives, random=~1|Mesocosm, family=binomial, data=repeatedm.stats12)
plot(z.binomial)
qqnorm(residuals(z.binomial));qqline(residuals(z.binomial), col="red")

plot(resid(z.binomial)~fitted(z.binomial))

summary(z.binomial)






#lmer: lmer(y~A*B + (1|plot)) ... or: lmer(y~CO2* Invasives + (1|Mesocosm))
library(lme4)
z5<-lmer(botryllid ~ CO2.Treatment*Invasives + (1|Mesocosm), REML=TRUE, data=repeatedm.stats12)
plot(z5)
summary(z5)
#The output of lmer is harder to intepret ---> doesn't present any pvalues

### Not sure this is what I want to test 
#I need a model where Mesocosm is random, CO2, Invasives and Week are fixed, and Invasives is nested within CO2 but that facto is random



# Using only linear trend of year
#shrubs$yearL <- as.numeric(as.character(shrubs$year))-2000  # center
# at year 2000
> print(model2<-lmer(response~yearL*snow*warm  + (1|plot), method="ML",
                     > data=shrubs), cor=FALSE)
>
  > # dropping 3-factor interactions, clearly suggests two significant
  > interactions with year
> print(model3<-lmer(response~(yearL+snow+warm)^2  + (1|plot),
                     > method="ML", data=shrubs), cor=FALSE)
> anova(model1, model2, model3)
>
  > # plot two significant interactions
  > shrubs.rs <- melt(shrubs, id=c("plot", "yearL", "warm", "snow"),
                      > measure="response")
>

z<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm.stats)
plot(z)
summary(z)
visreg(z, xvar = "CO2.Treatment", by = "Invasives", scales=list(rot = 90))





z2<-lme(botryllid ~ pH.uptowk*Invasives, random=~1|Mesocosm, data=repeatedm.stats)
plot(z2)
summary(z2)
anova(z2)
#plot looks better for by pH
visreg(z2, xvar = "pH.uptowk", by = "Invasives", scales=list(rot = 90))

anova(z)
anova(z2)
## do I use summary or anova - give different results - need to analyze differently










