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

repeatedm.present<-repeatedm[repeatedm$Invasives=="Present", ]

repeatedm.week14.present<-repeatedm.present[repeatedm.present$Week=="14", ]
head(repeatedm.week14.present)

windows()

# Botryllid # sub in any name
cdata <- summaryBy(botryllid ~ pH.uptowk + Treatment + CO2.Treatment, data=repeatedm.week14.present, FUN=c(length,mean,sd))
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

p<-ggplot(cdata, aes(x=pH.uptowk, y=botryllid.mean, colour=Treatment)) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "pH.uptowk", y="Botryllid percent cover")+ theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())






############# Just regular plotting with pH.uptowk

p<- ggplot(repeatedm.week14.present, aes(x=repeatedm.week14.present$pH.uptowk, y=repeatedm.week14.present$botryllid, colour=repeatedm.week14.present$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH.uptowk", y="botryllid percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p
#WIth lm model
p<- ggplot(repeatedm.week14.present, aes(x=repeatedm.week14.present$pH.uptowk, y=repeatedm.week14.present$botryllid)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH.uptowk", y="botryllid percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p





####LM
zbot<-lm(repeatedm.week14.present$botryllid ~ repeatedm.week14.present$pH.uptowk)
plot(zbot)
summary(zbot)

library(quantreg)

zbot.rq<-rq(botryllid ~ pH.uptowk, data=repeatedm.week14.present, tau=0.5)
r1.pH.uptowk<-resid(zbot.rq)

plot(x=repeatedm.week14.present$pH.uptowk, y=repeatedm.week14.present$botryllid, cex=.25, type="n", xlab="pH.uptowk week14", ylab="botryllid percent cover")
points(x=repeatedm.week14.present$pH.uptowk, y=repeatedm.week14.present$botryllid, cex=.75, col="blue")
abline(rq(botryllid~pH.uptowk,data=repeatedm.week14.present, tau=.5),col="blue")
abline(lm(botryllid~pH.uptowk,data=repeatedm.week14.present),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(botryllid~pH.uptowk,data=repeatedm.week14.present,tau=taus[i]),col="gray")}


summary(zbot.rq)
summary(zbot.rq, se="nid")
