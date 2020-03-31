mydata<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
#repeated measures file


zmem.eaten2<- lme(mem.eaten ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, data=mydata10)
# makes a funnel, so doesn't work 
summary(zmem.eaten2)
plot(zmem.eaten2)


head(mydata)
library(nlme)

mydata$mem.eaten.trans<-log(mydata$mem.eaten + 1)

zmem.eaten.trans<- lme(mem.eaten.trans ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, data=mydata10)
plot(zmem.eaten.trans)
summary(zmem.eaten.trans)

fix(mydata10)

mydata16<-mydata[mydata$Week==16, ]

mydata14<-mydata[mydata$Week==14, ]

mydata12<-mydata[mydata$Week==12, ]

mydata10<-mydata[mydata$Week==10, ]

mydata8<-mydata[mydata$Week==8, ]

mydata6<-mydata[mydata$Week==6, ]

mydata4<-mydata[mydata$Week==4, ]

mydata2<-mydata[mydata$Week==2, ]


library(MASS)
zmem.eaten<-glmmPQL(mem.eaten ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, family=binomial, data=mydata10)
plot(zmem.eaten)

windows(
summary(zmem.eaten)
anova(zmem.eaten)
qqnorm(zmem.eaten)


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





