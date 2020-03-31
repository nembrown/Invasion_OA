prop<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
names(prop)

hist(mydata$mem.eaten)
sort(mydata$mem.eaten)
fix(mydata)
names(mydata)
#install.packages("gamlss")
library(gamlss)
library(lme4)

is.factor(mydata$Tile)

mydata$Tile<-as.factor(mydata$Tile)	
mydata$Invasives<-as.factor(mydata$Invasives)	
mydata$CO2.Treatment<-as.factor(mydata$CO2.Treatment)	
mydata$Treatment<-as.factor(mydata$Treatment)



mydata$mem.eaten.trans<-log(mydata$mem.eaten + 1)
mydata$mem.eaten.cube<-mydata$mem.eaten^(1/3)

hist(mydata$mem.eaten[mydata$Week!=0])
hist(mydata$mem.eaten.cube)


zmem.eaten.trans<- lme(mem.eaten.trans ~ CO2.Treatment*Invasives, random=~1|Mesocosm/Invasives, data=mydata10)



fm01<-lmer(mem.eaten ~ CO2.Treatment*Invasives + 1|Mesocosm, data=mydata)
qqnorm(residuals(fm01));qqline(residuals(fm01), col="red")  #check for normality in the residuals


fm02<-lmer(mem.eaten.cube ~ CO2.Treatment*Invasives + 1|Mesocosm, data=mydata)
qqnorm(residuals(fm02));qqline(residuals(fm02), col="red")  #check for normality in the residuals

fm03<-lmer(mem.eaten.trans ~ CO2.Treatment*Invasives + 1|Mesocosm, data=mydata)
qqnorm(residuals(fm03));qqline(residuals(fm03), col="red")  #check for normality in the residuals


fm04<-lmer(mem.eaten.trans ~ CO2.Treatment*Invasives + Invasives|Mesocosm, data=mydata)
qqnorm(residuals(fm04));qqline(residuals(fm04), col="red")  #check for normality in the residuals

fm05<-lmer(mem.eaten.trans ~ CO2.Treatment*Invasives + Mesocosm|Invasives, data=mydata)
qqnorm(residuals(fm05));qqline(residuals(fm05), col="red")  #check for normality in the residuals
plot(fm05)






fm06<-lmer(mem.eaten ~ CO2.Treatment+ 1|Mesocosm, data=mydata)
qqnorm(residuals(fm06));qqline(residuals(fm06), col="red")  #check for normality in the residuals
summary(fm06)
