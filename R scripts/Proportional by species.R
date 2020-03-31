propbyspecies<-read.csv(file.choose())
#File: pH input .csv

head(propbyspecies)

library(plyr)

library(ggplot2)

library(doBy)

cdata <- summaryBy(Proportional.difference ~ Treatment + Species, data=propbyspecies, FUN=c(length,mean,sd))
head(cdata)
#summarize the mussel data by treatment and week in terms of mean and standard dev

names(cdata)[names(cdata)=="Proportional.difference.length"]<-"N"
head(cdata)
#changed column heading to N

cdata$Proportional.difference.se<-cdata$Proportional.difference.sd/sqrt(cdata$N)
head(cdata)
#created standard error from SD


p<-ggplot(cdata, aes(x=Species, y=Proportional.difference.mean, colour=Treatment)) + geom_errorbar(aes(ymin=Proportional.difference.mean-Proportional.difference.se, ymax=Proportional.difference.mean+Proportional.difference.se), colour="black", width=.1) + geom_line() + geom_point()
#rename my plot to p, makes it easier to work with

p + labs(x = "Species", y="Proportional.difference (NO Invasives - Invasives)")+theme_bw()