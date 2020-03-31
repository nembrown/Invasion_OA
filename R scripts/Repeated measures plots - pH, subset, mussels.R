repeatedm<-read.csv(file.choose())
#File: All week R input .csv

levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="AIR, Absent"] <- "Control, Invasives Absent"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="CO2, Absent"] <- "Elevated CO2, Invasives Absent"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="AIR, Present"] <- "Control, Invasives Present"
levels(repeatedm$Treatment)[levels(repeatedm$Treatment)=="CO2, Present"] <- "Elevated CO2, Invasives Present"
levels(repeatedm$CO2.Treatment)[levels(repeatedm$CO2.Treatment)=="AIR"] <- "Control"
levels(repeatedm$CO2.Treatment)[levels(repeatedm$CO2.Treatment)=="CO2"] <- "Elevated CO2"




head(repeatedm)
head(repeatedm)
repeatedm<-repeatedm[repeatedm$Invasives=="Present",]

repeatedm$Mesocosm<-as.factor(repeatedm$Mesocosm)
repeatedm$Week<-as.factor(repeatedm$Week)

repeatedm16<-repeatedm[repeatedm$Week==16, ]

repeatedm14<-repeatedm[repeatedm$Week==14, ]

repeatedm12<-repeatedm[repeatedm$Week==12, ]

repeatedm10<-repeatedm[repeatedm$Week==10, ]

repeatedm8<-repeatedm[repeatedm$Week==8, ]

repeatedm6<-repeatedm[repeatedm$Week==6, ]

repeatedm4<-repeatedm[repeatedm$Week==4, ]

repeatedm2<-repeatedm[repeatedm$Week==2, ]



library(plyr)

library(ggplot2)

library(doBy)


repeatedm<-repeatedm[which(repeatedm$Week>0), ]
head(repeatedm)

windows()

# botryllid # sub in any name
cdata <- summaryBy(botryllid ~ pH.uptowk + Treatment + Invasives + CO2.Treatment, data=repeatedm, FUN=c(length,mean,sd))
head(cdata)
#summarize the botryllid data by Treatment and week in terms of mean and standard dev

names(cdata)[names(cdata)=="botryllid.length"]<-"N"
head(cdata)
#changed column heading to N

cdata$botryllid.se<-cdata$botryllid.sd/sqrt(cdata$N)
head(cdata)
#created standard error from SD



p<-ggplot(cdata, aes(x=Week, y=botryllid.mean, linetype=Treatment, colour=pH.uptowk )) + geom_errorbar(aes(ymin=botryllid.mean-botryllid.se, ymax=botryllid.mean+botryllid.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "Time(weeks)", y="Mean % cover\nC. inflata")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=14)) + theme(axis.title.y = element_text(angle=0))

p + scale_linetype_manual(values=c("dashed","solid","dashed", "solid" ))+ theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

OR:
windows() # THis is a different thing coloured ---> this is one I use in my powerpoint

p<-ggplot(cdata, aes(x=pH.uptowk, y=botryllid.mean, colour=Invasives)) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "pH.uptowk", y="Shannon diversity\nindex (H)\nnative species")+ theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=0))

p<- p + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

p+ geom_smooth(method=lm)



####### LME 
library(nlme)
library(visreg)




z<-lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm, data=repeatedm)
plot(z)
summary(z)
visreg(z, xvar = "CO2.Treatment", by = "Invasives", scales=list(rot = 90))

z2<-lme(botryllid ~ pH.uptowk*Invasives, random=~1|Mesocosm, data=repeatedm)
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

zlqmm<-lqmm(botryllid ~ pH.uptowk*Invasives, random=~1, group= Mesocosm, tau = 0.5,
     nK = 7, type = "normal", rule = 1, data = repeatedm, na.action = na.fail, control = list(),
     contrasts = NULL, fit = TRUE)
summary(zlqmm)

plot(zlqmm)

### can I make this work for lqmm??

plot(x=repeatedm$pH.uptowk, y=repeatedm$botryllid, cex=.25, type="n", xlab="pH", ylab="botryllid of mussels")
points(x=repeatedm$pH.uptowk, y=repeatedm$botryllid, cex=.75, col="blue")
abline(rq(botryllid~pH.uptowk,data=repeatedm, tau=.5),col="blue")
abline(lm(botryllid~pH.uptowk,data=repeatedm),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(botryllid~pH,data=repeatedm,tau=taus[i]),col="gray")}




############# Just regular plotting with pH.uptowk

p<- ggplot(repeatedm10, aes(x=repeatedm$pH.uptowk, y=repeatedm$botryllid, colour=repeatedm$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH.uptowk", y="botryllid percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p



#WIth lm model
p<- ggplot(repeatedm10, aes(x=repeatedm$pH.uptowk, y=repeatedm$botryllid)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH.uptowk", y="botryllid percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p





####LM
zbot<-lm(repeatedm$botryllid ~ repeatedm$pH.uptowk)
plot(zbot)
summary(zbot)

library(quantreg)

zbot.rq<-rq(botryllid ~ pH.uptowk, data=repeatedm, tau=0.5)
r1.pH.uptowk<-resid(zbot.rq)

plot(x=repeatedm$pH.uptowk, y=repeatedm$botryllid, cex=.25, type="n", xlab="pH.uptowk", ylab="botryllid of botryllids")
points(x=repeatedm$pH.uptowk, y=repeatedm$botryllid, cex=.75, col="blue")
abline(rq(botryllid~pH.uptowk,data=repeatedm, tau=.5),col="blue")
abline(lm(botryllid~pH.uptowk,data=repeatedm),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(botryllid~pH.uptowk,data=repeatedm,tau=taus[i]),col="gray")}


summary(zbot.rq)
summary(zbot.rq, se="nid")
