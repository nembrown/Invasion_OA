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
repeatedm.stats<-repeatedm.stats[repeatedm.stats$Invasives=="Present",]

repeatedm.stats<-repeatedm.stats[which(repeatedm.stats$Week>0), ]
repeatedm.stats$Mesocosm<-as.factor(repeatedm.stats$Mesocosm)
repeatedm.stats$Week<-as.factor(repeatedm.stats$Week)
repeatedm.stats$Tile.ID<-as.factor(repeatedm.stats$Tile.ID)
repeatedm.stats$prop.mem.eaten2<-1-(repeatedm.stats$bare)


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
repeatedm.stats$prop.prop.mem.eaten2<-(repeatedm.stats$prop.mem.eaten2)/(repeatedm.stats$botryllid + repeatedm.stats$prop.mem.eaten2)
repeatedm.stats$prop.prop.mem.eaten2[is.na(repeatedm.stats$prop.prop.mem.eaten2)]<-0

repeatedm.stats$prop.mem.eaten2<-(repeatedm.stats$mem.eaten)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$prop.mem.eaten2[is.na(repeatedm.stats$prop.mem.eaten2)]<-0

repeatedm.stats$mem.dead<-(repeatedm.stats$mem.dead)/(repeatedm.stats$mussel + repeatedm.stats$mem.eaten + repeatedm.stats$mem.dead)
repeatedm.stats$mem.dead2[is.na(repeatedm.stats$prop.mem.dead2)]<-0

repeatedm.stats$total.bot<-(repeatedm.stats$botryllid + repeatedm.stats$prop.mem.eaten2)
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

zlme1<-lme(prop.mem.eaten2 ~ CO2.Treatment*Week, random=~1|Tile.ID, data=repeatedm.stats, method="REML")
plot(zlme1)
summary(zlme1)



#If I use repeatedm.stats (i.e. with 0), then the summary produces much less text --> is it because I specified as factor earlier? 

repeatedm.stats$Week.factor<-as.factor(repeatedm.stats$Week)

head(repeatedm.stats)
z6<-lme(prop.mem.eaten2  ~ av.pH.all*Week, random=~1|Tile.ID, data=repeatedm.stats)
plot(z6)
summary(z6)

### I think as soon as I treat week as a factor then the residuals are in lines.... don't look good/// form lines. 
### Can I make this a glmm??
# YEs. Either GLMM full model or GLMM with just week ... remember still need the residuals to look okay with that 
#model distribution (e.g. binomial), check if the model fits the data. 

#going to work with glmmpql b/c that's what zuur uses and is suggested in the 2008 Bolker TREE paper ...they also talk about lmer from the lme4 package and
# glmmML from the package of the same name


library(MASS)
prop.mem.eaten2.PQL<- glmmPQL(prop.mem.eaten2 ~ pH.uptowk, random=~1|Tile.ID, family=binomial, data=repeatedm.stats14)
summary(prop.mem.eaten2.PQL)
plot(prop.mem.eaten2.PQL)




#Residuals still have problems, because of week. Is this the right function to plot residuals? 
qqnorm(residuals(prop.mem.eaten2.PQL));qqline(residuals(prop.mem.eaten2.PQL), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(prop.mem.eaten2.PQL)~repeatedm.stats$Week)


###try negative binomial
# run first non mixed to get theta value:
prop.mem.eaten2.glm.nb<- glm.nb(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, data=repeatedm.stats)
summary(prop.mem.eaten2.glm.nb)

prop.mem.eaten2.PQL.nb2<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile, negative.binomial(theta =227040 , link = log), data=repeatedm.stats)
plot(prop.mem.eaten2.PQL.nb2)
summary(prop.mem.eaten2.PQL.nb2)

###### Looking up repeated measures online -> need to include week as random: 
###(Day | Subject),
prop.mem.eaten2.PQL<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1|Week/Mesocosm, family=poisson, data=repeatedm.stats)
plot(prop.mem.eaten2.PQL)
summary(prop.mem.eaten2.PQL)


repeatedm.stats$Week<-as.numeric(repeatedm.stats$Week)

#### Not sure how this works...
##Can also include Tile (i.e. individual/subject as a random factor ******** I think this is important)
prop.mem.eaten2.PQL2<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile.ID, family=poisson, data=repeatedm.stats)
plot(prop.mem.eaten2.PQL2)
summary(prop.mem.eaten2.PQL2)

#Model validation: (from https://danieljhocking.wordpress.com/tag/glmmpql/)
# Check for residual pattern within groups and difference between groups     

#might be useful to plot by random factor
library(lattice)
xyplot(residuals(prop.mem.eaten2.PQL) ~ fitted(prop.mem.eaten2.PQL) | repeatedm.stats$Mesocosm)
xyplot(residuals(prop.mem.eaten2.PQL) ~ fitted(prop.mem.eaten2.PQL)) 
#residuals in general look bad. 
plot.new
plot(fitted(prop.mem.eaten2.PQL), residuals(prop.mem.eaten2.PQL),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(prop.mem.eaten2.PQL), residuals(prop.mem.eaten2.PQL)))




  #Without week, residuals look okay     
PQL.nowk<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm, family=poisson, data=repeatedm.stats)
plot(PQL.nowk)
plot(resid(PQL.nowk)~repeatedm.stats$Week)






## Try lmer instead of glmmPQL - laplace approximation, wald's Z
library(lme4)
lmer1<- glmer(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer1)
plot(lmer1)
qqnorm(residuals(lmer1));qqline(residuals(lmer1), col="red")


lmer2<- glmer(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week + (1|Mesocosm/Tile), family=poisson, data=repeatedm.stats)
summary(lmer1)
plot(lmer1)
#When I include week not sure if that's okay? 
#What's the best for repeated measures? 
#What is the best way to evaluate the model, model assumptions etc ... how do I know if it's a good model for my data? 

lmer2<- glmer(prop.mem.eaten2 ~ CO2.Treatment*Invasives + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer2)

anova(lmer1, lmer2, test = "Chisq")
#So obviously week is important

lmer3<- glmer(prop.mem.eaten2 ~ CO2.Treatment*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer3)
anova(lmer1, lmer3, test = "Chisq")
# so obviously invasives are important

lmer4<- glmer(prop.mem.eaten2 ~ Invasives*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer4)
plot(lmer4)
anova(lmer1, lmer4, test = "Chisq")
# so CO2 not important 

lmer5<- glmer(prop.mem.eaten2 ~ CO2.Treatment:Invasives  + Invasives*Week + (1|Mesocosm), family=poisson, data=repeatedm.stats)
summary(lmer5)
anova(lmer1, lmer5, test = "Chisq")

#am concerned about dropping terms with nesting.




## Try glmmML instead of glmmPQL
library(glmmML)
glmmML1<-glmmML(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, cluster=Mesocosm, family=poisson, data=repeatedm.stats)
summary(glmmML1)
#this is odd becasue all significant now??? 




### Try glmm ADMB
head(repeatedm.stats)
repeatedm.stats$Mesocosm<-as.factor(repeatedm.stats$Mesocosm)
repeatedm.stats$Tile<-as.factor(repeatedm.stats$Tile)

install.packages("glmmADMB", repos="http://r-forge.r-project.org", type="source")
install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)
#glmmadmb(formula = NCalls ~ (FoodTreatment + ArrivalTime) * SexParent + offset(logBroodSize) + (1 | Nest), data = Owls, family = "poisson",zeroInflation = TRUE) 

glmmadmb<-glmmadmb(formula = prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week + (1 | Mesocosm/Tile.ID), data = repeatedm.stats, family = "poisson",zeroInflation = FALSE) 
plot(glmmadmb)
summary(glmmadmb)

is.factor(repeatedm.stats$Week)


glmmadmb<-glmmadmb(formula = prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week + (1 | Mesocosm/Tile.ID), data = repeatedm.stats, family = "poisson",zeroInflation = TRUE) 
summary(glmmadmb)





###########################################################################
###Model validation - Graphing distributions
 
#Residuals still have problems, because of week. Is this the right function to plot residuals? 
qqnorm(residuals(glmmadmb));qqline(residuals(glmmadmb), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(glmmadmb)~repeatedm.stats$Week)
plot(resid(glmmadmb)~fitted(glmmadmb))
#Makes a pattern

library(boot)
glm.diag.plots(prop.mem.eaten2.PQL2, glm.diag(prop.mem.eaten2.PQL2))
# can't use for glmm, only glm I think



#From https://danieljhocking.wordpress.com/2011/07/18/model-validation-interpreting-residual-plots/#comments
#Took me a while to figure this out too. I suggest you generate negative binomial quantiles: nbquant<-rnbinom(n=length(count), scale=1.480, mu=mean(count))
#Then plot those quantiles against your residuals:qqnorm(nbquant, resid(glmmPQLnb1)I think you'll find that your points fall along a straight line
#, suggesting that you did choose the right distribution. I had this issue too, took me a while to figure out, and once I did, I got a nice 
#straight line.


nbquant<-rnbinom(n=length(repeatedm.stats$prop.mem.eaten2), mu=mean(repeatedm.stats$prop.mem.eaten2), size=5)
### not sure what size needs to be - it says number of succcesses?? 
qqplot(nbquant, resid(prop.mem.eaten2.PQL2))
#what is this plotting??? SHould it be a straight line if perfect fit?? 

#fitdirst ## avoid spurious accuracy

x <- rpois(24, 1)
fitdistr(x, "Poisson")
## now do this directly with more control.
fitdistr(x, dpois, list(shape = 1, rate = 0.1), lower = 0.001)
#### Don't understand what's going on here. 



####################################################################################################################


##### TRY JUST ONE WEEK at a time: 

head(repeatedm.stats)
head(repeatedm.stats14)
z4<-lme(prop.mem.eaten2 ~ CO2.Treatment, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats14)
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


z5<-lme(prop.mem.eaten2 ~ CO2.Treatment, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats14, weights=varPower())
plot(z5)
## included varPower - still bad ... try var ident: 
#vf<-varIdent(form=~1|Treatment)

vf<-varIdent(form=~1|CO2.Treatment)
z6<-lme(prop.mem.eaten2 ~ CO2.Treatment, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats14, weights=vf)
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
z7<-lme(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats14, weights=vf2)
plot(z7)
summary(z7)

##check AIC to see which is better model
anova(z6,z7)
anova(z5,z6)
anova(z6,z4)
#### ns for the last one, so z6 is the best model I think. ********** Need to figure out if this is okay?? 


### Plot it out
ggplot(repeatedm.stats14, aes(x=Treatment, y=fitted(z4))) +geom_boxplot()+ scale_colour_brewer(palette="Set1")

ggplot(repeatedm.stats14, aes(x=Treatment, y=prop.mem.eaten2)) +geom_boxplot()+ scale_colour_brewer(palette="Set1")
#rename my plot to p, makes it easier to work with

#################################33
#Transformations
z8<-lme(log(prop.mem.eaten2+1) ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, data=repeatedm.stats14)
plot(z8)
summary(z8)
### doesn't look good... and remember: don't log transform count data, Bitches!!!! (Jared's blog post)




############################33

## GLMMS: 

############Poisson/ neg binomial

z.PQL<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, family=poisson, data=repeatedm.stats14)
plot(z.PQL)
qqnorm(residuals(z.PQL));qqline(residuals(z.PQL), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(z.PQL)~fitted(z.PQL))
#residuals look a lot better 
summary(z.PQL)

head(repeatedm.stats14)

## TRy negative binomial ... not better ... why is theta so huge? The plot function works better here... I think poisson is okay... 
z.nb<- glm.nb(prop.mem.eaten2 ~ CO2.Treatment*Invasives,data=repeatedm.stats14)
summary(z.nb)
plot(z.nb)

z.nb2<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, negative.binomial(theta =21218   , link = log), data=repeatedm.stats14)
plot(z.nb2)
summary(z.nb2)



###### Problem that I'm getting each time is that the residuals vs fitted values are heteroskedastic.... 
###estimates seem really off ... not even close - say it's 1.76 ... but really like 5 ?? 

ggplot(repeatedm.stats14, aes(x=Treatment, y=fitted(z.nb2))) +geom_boxplot()+ scale_colour_brewer(palette="Set1")

ggplot(repeatedm.stats14, aes(x=Treatment, y=prop.mem.eaten2)) +geom_boxplot()+ scale_colour_brewer(palette="Set1")





############## GAMMA
head(repeatedm.stats14)
z.gamma<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Tile.ID, family=Gamma, data=repeatedm.stats14)
plot(z.gamma)
qqnorm(residuals(z.gamma));qqline(residuals(z.gamma), col="red")
#Jocelyn's checking graphs ---> but is the normality of the residuals really what I want to know?? 
# For GLM, you still need to have residuals that make sense - but they don't have to be normal I dont think. 
plot(resid(z.gamma)~fitted(z.gamma))
#residuals look a lot better 
summary(z.gamma)


############ BINOMIAL
head(repeatedm.stats14)
z.binomial<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment, random=~1|Mesocosm/Tile.ID, family=binomial, data=repeatedm.stats14)
plot(z.binomial)
qqnorm(residuals(z.binomial));qqline(residuals(z.binomial), col="red")

plot(resid(z.binomial)~fitted(z.binomial))

summary(z.binomial)






#lmer: lmer(y~A*B + (1|plot)) ... or: lmer(y~CO2* Invasives + (1|Mesocosm))

z5<-lmer(prop.mem.eaten2 ~ CO2.Treatment*Invasives + (1|Mesocosm), REML=TRUE, data=repeatedm.stats)
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

z<-lme(prop.mem.eaten2 ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm.stats)
plot(z)
summary(z)
visreg(z, xvar = "CO2.Treatment", by = "Invasives", scales=list(rot = 90))





z2<-lme(prop.mem.eaten2 ~ pH.uptowk*Invasives, random=~1|Mesocosm, data=repeatedm.stats)
plot(z2)
summary(z2)
anova(z2)
#plot looks better for by pH
visreg(z2, xvar = "pH.uptowk", by = "Invasives", scales=list(rot = 90))

anova(z)
anova(z2)
## do I use summary or anova - give different results - need to analyze differently


############################

### Quantile regression without random factors ... not useful really unless don't have Invasives as a factor
zbot<-lm(repeatedm.stats$prop.mem.eaten2 ~ repeatedm.stats$av.pH.all)
plot(zbot)
summary(zbot)

library(quantreg)

zbot.rq<-rq(prop.mem.eaten2 ~ av.pH.all, data=repeatedm.stats, tau=0.5)
r1.av.pH.all<-resid(zbot.rq)

plot(x=repeatedm.stats$av.pH.all, y=repeatedm.stats$prop.mem.eaten2, cex=.25, type="n", xlab="av.pH.all", ylab="prop.mem.eaten2 of prop.mem.eaten2s")
points(x=repeatedm.stats$av.pH.all, y=repeatedm.stats$prop.mem.eaten2, cex=.75, col="blue")
abline(rq(prop.mem.eaten2~av.pH.all,data=repeatedm.stats, tau=.5),col="blue")
abline(lm(prop.mem.eaten2~av.pH.all,data=repeatedm.stats),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(prop.mem.eaten2~av.pH.all,data=repeatedm.stats,tau=taus[i]),col="gray")}


summary(zbot.rq)
summary(zbot.rq, se="nid")


### So quantile regression will probably only be feasible if I DON'T have invasives in there ... i.e. in things that involve prop.mem.eaten2s
#Becaus ethen I can discount the prop.mem.eaten2s and then I can do just by CO2 or JUST by pHuptowk. 
#Also will only be feasible if I do ONE WEEK 
#Theremore might be useful for prop.mem.eaten2 percent cover, bot. eaten at certain weeks. 

library(quantreg)
head(repeatedm.stats14)

Rm.rq<-rq(prop.mem.eaten2 ~ pH.uptowk, data=repeatedm.stats14, tau=0.50)
r1.pH<-resid(Rm.rq)
plot(Rm.rq)

hist(r1.pH)
c1.pH<-coef(Rm.rq)
summary(Rm.rq)
summary(Rm.rq, se="nid")
plot(x=repeatedm.stats14$pH.uptowk, y=repeatedm.stats14$prop.mem.eaten2, cex=.25, type="n", xlab="pH", ylab="bot eaten")
points(x=repeatedm.stats14$pH.uptowk, y=repeatedm.stats14$prop.mem.eaten2, cex=.75, col="blue")
abline(rq(prop.mem.eaten2 ~ pH.uptowk,data=repeatedm.stats14, tau=.5),col="blue")
abline(lm(prop.mem.eaten2 ~ pH.uptowk,data=repeatedm.stats14),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats14,tau=taus[i]),col="gray")}









