repeatedm<-read.csv(file.choose())
#File: All week R input .csv

levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="AIR, Absent"] <- "Control, Invasives Absent"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="CO2, Absent"] <- "Elevated CO2, Invasives Absent"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="AIR, Present"] <- "Control, Invasives Present"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="CO2, Present"] <- "Elevated CO2, Invasives Present"
levels(repeatedm$CO2.Treatment)[levels(repeatedm$CO2.Treatment)=="AIR"] <- "Control"
levels(repeatedm$CO2.Treatment)[levels(repeatedm$CO2.Treatment)=="CO2"] <- "Elevated CO2"


head(repeatedm)



library(plyr)

library(ggplot2)

library(doBy)


repeatedm.nowk0<-repeatedm[which(repeatedm$Week>0), ]
head(repeatedm.nowk0)

windows()

# shannon.diversity.no.bot # sub in any name
cdata <- summaryBy(shannon.diversity.no.bot ~ pH.uptowk + Treatment + Invasives + CO2.Treatment, data=repeatedm.nowk0, FUN=c(length,mean,sd))
head(cdata)
#summarize the shannon.diversity.no.bot data by Treatment and week in terms of mean and standard dev

names(cdata)[names(cdata)=="shannon.diversity.no.bot.length"]<-"N"
head(cdata)
#changed column heading to N

cdata$shannon.diversity.no.bot.se<-cdata$shannon.diversity.no.bot.sd/sqrt(cdata$N)
head(cdata)
#created standard error from SD



p<-ggplot(cdata, aes(x=Week, y=shannon.diversity.no.bot.mean, linetype=Treatment, colour=pH.uptowk )) + geom_errorbar(aes(ymin=shannon.diversity.no.bot.mean-shannon.diversity.no.bot.se, ymax=shannon.diversity.no.bot.mean+shannon.diversity.no.bot.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "Time(weeks)", y="Mean % cover\nC. inflata")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=14)) + theme(axis.title.y = element_text(angle=0))

p + scale_linetype_manual(values=c("dashed","solid","dashed", "solid" ))+ theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

OR:
windows() # THis is a different thing coloured ---> this is one I use in my powerpoint

p<-ggplot(cdata, aes(x=pH.uptowk, y=shannon.diversity.no.bot.mean, colour=Invasives)) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "pH.uptowk", y="Shannon diversity\nindex (H)\nnative species")+ theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=0))

p<- p + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

p+ geom_smooth(method=lm)



####### LME 
library(nlme)
library(visreg)
library(lme4)

head(repeatedm.nowk0)
repeatedm.nowk0$Week<-as.factor(repeatedm.nowk0$Week)
repeatedm.nowk0$Mesocosm<-as.factor(repeatedm.nowk0$Mesocosm)


##### Trying ANCOVA type design.....
# eed to figure ou the proper syntax. 

### ANCOVA uses pH


zancova<- lme(num.species.no.bot ~ av.pH.all*Invasives*Week, random=~1|Mesocosm/Invasives, data=repeatedm.nowk0)
summary(zancova)





head(repeatedm.nowk0)
z4<-lme(shannon.diversity.no.bot ~ Week*CO2.Treatment*Invasives, random=~1|Mesocosm/CO2.Treatment, data=repeatedm.nowk0)
summary(z4)
plot(z4)

# lmer equivalen, with polynomial contrasts for year and treatment
#  contrast for snow (reference=A)
# print(model1<-lmer(response~ year*snow*warm + (1|plot), method="ML",
  #                   > data=shrubs), cor=FALSE)


znowk0<-lmer(num.species.no.bot ~ Week*CO2.Treatment*Invasives + (1|Mesocosm), REML=TRUE, data=repeatedm.nowk0)
summary(znowk0)

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

z<-lme(shannon.diversity.no.bot ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm.nowk0)
plot(z)
summary(z)
visreg(z, xvar = "CO2.Treatment", by = "Invasives", scales=list(rot = 90))





z2<-lme(shannon.diversity.no.bot ~ pH.uptowk*Invasives, random=~1|Mesocosm, data=repeatedm.nowk0)
plot(z2)
summary(z2)
anova(z2)
#plot looks better for by pH
visreg(z2, xvar = "pH.uptowk", by = "Invasives", scales=list(rot = 90))

anova(z)
anova(z2)
## do I use summary or anova - give different results - need to analyze differently








@#####Yes there is quantile regression using lqmm for formulas with random effects but I'm finding difficult o visiualize residuals etc ... 
  #So might not be the best ... maybe just try a glm ... try getting Week in there. 
### residuals not great -- maybe can try quantile regression?? lqmm
library(lqmm)
#lqmm(fixed, random, group, covariance = "pdDiag", tau = 0.5,
    # nK = 7, type = "normal", rule = 1, data = sys.frame(sys.parent()),
   #  subset, weights, na.action = na.fail, control = list(),
  #   contrasts = NULL, fit = TRUE)

zlqmm<-lqmm(shannon.diversity.no.bot ~ pH.uptowk*Invasives, random=~1, group= Mesocosm, tau = 0.5,
     nK = 7, type = "normal", rule = 1, data = repeatedm.nowk0, na.action = na.fail, control = list(),
     contrasts = NULL, fit = TRUE)
summary(zlqmm)

plot(zlqmm)

### can I make this work for lqmm??

plot(x=repeatedm$pH.uptowk, y=repeatedm$shannon.diversity.no.bot, cex=.25, type="n", xlab="pH", ylab="shannon.diversity.no.bot of mussels")
points(x=repeatedm$pH.uptowk, y=repeatedm$shannon.diversity.no.bot, cex=.75, col="blue")
abline(rq(shannon.diversity.no.bot~pH.uptowk,data=repeatedm, tau=.5),col="blue")
abline(lm(shannon.diversity.no.bot~pH.uptowk,data=repeatedm),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(shannon.diversity.no.bot~pH,data=repeatedm,tau=taus[i]),col="gray")}




############# Just regular plotting with pH.uptowk

p<- ggplot(repeatedm.nowk0, aes(x=repeatedm.nowk0$pH.uptowk, y=repeatedm.nowk0$shannon.diversity.no.bot, colour=repeatedm.nowk0$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH.uptowk", y="shannon.diversity.no.bot percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p



#WIth lm model
p<- ggplot(repeatedm.nowk0, aes(x=repeatedm.nowk0$pH.uptowk, y=repeatedm.nowk0$shannon.diversity.no.bot)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH.uptowk", y="shannon.diversity.no.bot percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p





####LM
zbot<-lm(repeatedm.nowk0$shannon.diversity.no.bot ~ repeatedm.nowk0$pH.uptowk)
plot(zbot)
summary(zbot)

library(quantreg)

zbot.rq<-rq(shannon.diversity.no.bot ~ pH.uptowk, data=repeatedm.nowk0, tau=0.5)
r1.pH.uptowk<-resid(zbot.rq)

plot(x=repeatedm.nowk0$pH.uptowk, y=repeatedm.nowk0$shannon.diversity.no.bot, cex=.25, type="n", xlab="pH.uptowk", ylab="shannon.diversity.no.bot of shannon.diversity.no.bots")
points(x=repeatedm.nowk0$pH.uptowk, y=repeatedm.nowk0$shannon.diversity.no.bot, cex=.75, col="blue")
abline(rq(shannon.diversity.no.bot~pH.uptowk,data=repeatedm.nowk0, tau=.5),col="blue")
abline(lm(shannon.diversity.no.bot~pH.uptowk,data=repeatedm.nowk0),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(shannon.diversity.no.bot~pH.uptowk,data=repeatedm.nowk0,tau=taus[i]),col="gray")}


summary(zbot.rq)
summary(zbot.rq, se="nid")
