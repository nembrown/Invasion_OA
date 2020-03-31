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


repeatedm8<-repeatedm[repeatedm$Week==8, ]

windows()

# Botryllid # sub in any name
cdata <- summaryBy(membranipora ~ CO2.Treatment + botryllid + Week, data=repeatedm8, FUN=c(length,mean,sd))
head(cdata)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata)[names(cdata)=="membranipora.length"]<-"N"
head(cdata)
#changed column heading to N

cdata$membranipora.se<-cdata$membranipora.sd/sqrt(cdata$N)
head(cdata)
#created standard error from SD



p<-ggplot(cdata, aes(x=botryllid, y=membranipora.mean, colour=CO2.Treatment )) + geom_errorbar(aes(ymin=membranipora.mean-membranipora.se, ymax=membranipora.mean+membranipora.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with


###############
repeatedm16<-repeatedm[repeatedm$Week==16, ]
repeatedm14<-repeatedm[repeatedm$Week==14, ]
repeatedm12<-repeatedm[repeatedm$Week==12, ]
repeatedm10<-repeatedm[repeatedm$Week==10, ]
repeatedm8<-repeatedm[repeatedm$Week==8, ]
repeatedm6<-repeatedm[repeatedm$Week==6, ]
repeatedm4<-repeatedm[repeatedm$Week==4, ]




ggplot(repeatedm8, aes(x=botryllid, y=membranipora)) + geom_point(shape=1) + geom_smooth(method=lm) 


OR

ggplot(repeatedm10, aes(x=bot        ryllid, y=mussel, color=CO2.Treatment)) + geom_point(shape=1) + scale_colour_hue(l=50) + geom_smooth(method=lm,se=FALSE) 






##########################

p<- p + labs(x = "Percent cover Botryllid", y="Mean % cover\nMembranipora")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=14)) + theme(axis.title.y = element_text(angle=0))

p + scale_linetype_manual(values=c("dashed","solid","dashed", "solid" ))+ theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

OR:
windows()


ggplot(repeatedm8, aes(x=botryllid, y=membranipora)) + geom_point(shape=1) + geom_smooth(method=lm) 







p<-ggplot(cdata, aes(x=Week, y=membranipora.mean, colour=Treatment, linetype=CO2.Treatment )) + geom_errorbar(aes(ymin=membranipora.mean-membranipora.se, ymax=membranipora.mean+membranipora.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=3.5)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "Time(weeks)", y="Mean Species Richness")+ theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=0))

p + scale_colour_manual(values=c("#FF7256","#00CDCD","#FF7256","#00CDCD"))+ theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())




 + scale_shape_manual(breaks=c("AIR","CO2"), values=c(16,17))

+  scale_linetype_manual(values=c("dashed", "dashed", "solid", solid"))






#-----------------------------------------------------------------------
#Bryozoan

cdata1 <- summaryBy(av.bryo ~ Treatment + Week, data=repeatedm, FUN=c(length,mean,sd))
head(cdata1)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata1)[names(cdata1)=="av.bryo.length"]<-"N"
head(cdata1)
#changed column heading to N

cdata1$av.bryo.se<-cdata1$av.bryo.sd/sqrt(cdata1$N)
head(cdata1)
#created standard error from SD

p1<-ggplot(cdata1, aes(x=Week, y=av.bryo.mean, colour=Treatment)) + geom_errorbar(aes(ymin=av.bryo.mean-av.bryo.se, ymax=av.bryo.mean+av.bryo.se), colour="black", width=.1) + geom_line(size=0.5) + geom_point(size=2.5)
#rename my plot to p, makes it easier to work with

p1<- p1 + labs(x = "Time(weeks)", y="Mean % cover\nM. membranacea")+theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=12)) + theme(axis.title.y = element_text(angle=0))

p1 + theme(legend.text = element_text(colour="black", size = 14))+theme(legend.title = element_text(colour="black", size=16))





#------------------------------------------------------------------------
#Bare space
cdata2 <- summaryBy(av.membranipora ~ Treatment + Week, data=repeatedm, FUN=c(length,mean,sd))
head(cdata2)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata2)[names(cdata2)=="av.membranipora.length"]<-"N"
head(cdata2)
#changed column heading to N

cdata2$av.membranipora.se<-cdata2$av.membranipora.sd/sqrt(cdata2$N)
head(cdata2)
#created standard error from SD

p2<-ggplot(cdata2, aes(x=Week, y=av.membranipora.mean, colour=Treatment)) + geom_errorbar(aes(ymin=av.membranipora.mean-av.membranipora.se, ymax=av.membranipora.mean+av.membranipora.se), colour="black", width=.1) + geom_line() + geom_point()
#rename my plot to p, makes it easier to work with

p2 + labs(x = "Time(weeks)", y="Mean percent cover of free space")+theme_bw() + ylim(0,100)

#------------------------------------------------------------------------
#Hydroid
cdata3 <- summaryBy(av.hydriod ~ Treatment + Week, data=repeatedm, FUN=c(length,mean,sd))
head(cdata3)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata3)[names(cdata3)=="av.hydriod.length"]<-"N"
head(cdata3)
#changed column heading to N

cdata3$av.hydriod.se<-cdata3$av.hydriod.sd/sqrt(cdata3$N)
head(cdata3)
#created standard error from SD

p3<-ggplot(cdata3, aes(x=Week, y=av.hydriod.mean, colour=Treatment)) + geom_errorbar(aes(ymin=av.hydriod.mean-av.hydriod.se, ymax=av.hydriod.mean+av.hydriod.se), colour="black", width=.1) + geom_line(size=0.5) + geom_point(size=2.5)
#rename my plot to p, makes it easier to work with

windows()
p3<-p3 + labs(x = "Time(weeks)", y="Mean % cover\nO. dichotoma")+theme_bw()+ theme(text = element_text(size=14), axis.text = element_text(size=12)) + theme(axis.title.y = element_text(angle=0))

p3 + theme(legend.text = element_text(colour="black", size = 14))+theme(legend.title = element_text(colour="black", size=16))




#
#----------------------------------------------------------------------------
Corella
cdata5 <- summaryBy(av.membranipora ~ Treatment + Week, data=repeatedm, FUN=c(length,mean,sd))
head(cdata5)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata5)[names(cdata5)=="av.membranipora.length"]<-"N"
head(cdata5)
#changed column heading to N

cdata5$av.membranipora.se<-cdata5$av.membranipora.sd/sqrt(cdata5$N)
head(cdata5)
#created standard error from SD

p5<-ggplot(cdata5, aes(x=Week, y=av.membranipora.mean, colour=Treatment)) + geom_errorbar(aes(ymin=av.membranipora.mean-av.membranipora.se, ymax=av.membranipora.mean+av.membranipora.se), colour="black", width=.1) + geom_line() + geom_point()
#rename my plot to p, makes it easier to work with

windows()
p5 + labs(x = "Time(weeks)", y="Mean percent cover of C. inflata")+theme_bw()

#----------------------------------------------------------------------------
Botryllus
cdata6 <- summaryBy(av.bottrylus ~ Treatment + Week, data=repeatedm, FUN=c(length,mean,sd))
head(cdata6)
#summarize the membranipora data by treatment and week in terms of mean and standard dev

names(cdata6)[names(cdata6)=="av.bottrylus.length"]<-"N"
head(cdata6)
#changed column heading to N

cdata6$av.bottrylus.se<-cdata6$av.bottrylus.sd/sqrt(cdata6$N)
head(cdata6)
#created standard error from SD

p6<-ggplot(cdata6, aes(x=Week, y=av.bottrylus.mean, colour=Treatment)) + geom_errorbar(aes(ymin=av.bottrylus.mean-av.bottrylus.se, ymax=av.bottrylus.mean+av.bottrylus.se), colour="black", width=.1) + geom_line() + geom_point()
#rename my plot to p, makes it easier to work with

windows()
p6 + labs(x = "Time(weeks)", y="Mean percent cover of B. schlosseri")+theme_bw()
