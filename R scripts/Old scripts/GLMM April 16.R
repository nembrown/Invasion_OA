mydata<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
#repeated measures file

head(mydata)
library(nlme)
library(lme4)

zbotryllid<- lme(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, data=mydata10)
# makes a funnel, so doesn't work 
summary(zbotryllid)
plot(zbotryllid)



mydata$botryllid.trans<-log(mydata$botryllid + 1)

zbotryllid.trans<- lme(botryllid.trans ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, data=mydata10)
plot(zbotryllid.trans)
summary(zbotryllid.trans)

fix(mydata10)

mydata16<-mydata[mydata$Week==16, ]

mydata14<-mydata[mydata$Week==14, ]

mydata12<-mydata[mydata$Week==12, ]

mydata10<-mydata[mydata$Week==10, ]

mydata8<-mydata[mydata$Week==8, ]

mydata6<-mydata[mydata$Week==6, ]

mydata4<-mydata[mydata$Week==4, ]

mydata2<-mydata[mydata$Week==2, ]


#Try including week

zbot<-glmer(botryllid ~ CO2.Treatment*Invasives + Week + (1|Mesocosm/Invasives), family=poisson, data=mydata)
plot(zbot)

#drop the CO2 treatment x invasives term

zbot2<-glmer(botryllid ~ CO2.Treatment + Invasives + Week + (1|Mesocosm/Invasives), family=poisson, data=mydata)
plot(zbot2)

anova(zbot,zbot2)








library(MASS)
zbotryllid<-glmmPQL(botryllid ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=binomial, data=mydata10)
plot(zbotryllid)

windows(
summary(zbotryllid)
anova(zbotryllid)
qqnorm(zbotryllid)


qqnorm(zbotryllus)


znumspecgauss<-glmmPQL(num.species ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=gaussian ,data=mydata2)
plot(znumspecgauss)
summary(znumspecgauss)


znumspecgamma<-glmmPQL(num.species ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=Gamma ,data=mydata14)
windows()
plot(znumspecgamma)
summary(znumspecgamma)


zshannon.diversity<-glmmPQL(shannon.diversity ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=Gamma ,data=mydata4)
plot(zshannon.diversity)
summary(zshannon.diversity)


#-------------------
#count data
zshannon.diversity<-glmmPQL(shannon.diversity ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=Gamma ,data=mydata8)
plot(zshannon.diversity)
summary(zshannon.diversity)





