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


windows()

 
head(repeatedm)

repeatedm$occupied<-100-(repeatedm$bare)
repeatedm$botryllid<-(repeatedm$occupied - repeatedm$botryllid)
repeatedm$botryllid<-repeatedm$botryllid/(repeatedm$occupied)
repeatedm$botryllid[is.na(repeatedm$botryllid)]<-0


repeatedm$botryllid<-(repeatedm$bot.eaten)/(repeatedm$botryllid + repeatedm$bot.eaten)
repeatedm$prop.bot.eaten[is.na(repeatedm$prop.bot.eaten)]<-0

repeatedm$botryllid<-(repeatedm$mem.eaten)/(repeatedm$membranipora + repeatedm$mem.eaten + repeatedm$mem.dead)
repeatedm$botryllid[is.na(repeatedm$botryllid)]<-0

repeatedm$prop.mem.dead<-(repeatedm$mem.dead)/(repeatedm$membranipora + repeatedm$mem.dead)
repeatedm$mem.dead[is.na(repeatedm$prop.mem.dead2)]<-0

repeatedm$total.bot<-(repeatedm$botryllid + repeatedm$bot.eaten)

head(repeatedm)

repeatedm$propcover.mem<-(repeatedm$botryllid)/(repeatedm$occupied)
repeatedm$propcover.bot[is.na(repeatedm$botryllid)]<-0

repeatedm$propcover.mem<-(repeatedm$membranipora)/(repeatedm$occupied)
repeatedm$botryllid[is.na(repeatedm$botryllid)]<-0

repeatedm$propcover.proto<-(repeatedm$protozoa)/(repeatedm$occupied)
repeatedm$botryllid[is.na(repeatedm$botryllid)]<-0

repeatedm$propcover.mussel<-(repeatedm$mussel)/(repeatedm$occupied)
repeatedm$botryllid[is.na(repeatedm$botryllid)]<-0


repeatedm$native.occupied<-(repeatedm$occupied)-(repeatedm$botryllid)
repeatedm$prop.native.occupied<-(repeatedm$native.occupied)/(repeatedm$occupied)


repeatedm$nat.propcover.mem<-(repeatedm$membranipora)/(repeatedm$native.occupied)
repeatedm$nat.propcover.proto<-(repeatedm$protozoa)/(repeatedm$native.occupied)
repeatedm$nat.propcover.mussel<-(repeatedm$mussel)/(repeatedm$native.occupied)
repeatedm$barn<-(repeatedm$barn)/(repeatedm$native.occupied)



## handling Nas:
length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}
#sub in: na.rm=TRUE and length2 in cdata (remember serpulid.length2)

# botryllid # sub in any name
cdata <- summaryBy(botryllid ~ CO2.Treatment + Invasives + Treatment + Week , data=repeatedm, FUN=c(length2,mean,sd), na.rm=TRUE)
head(cdata)
#summarize the botryllid data by treatment and week in terms of mean and standard dev

names(cdata)[names(cdata)=="botryllid.length2"]<-"N"
head(cdata)
#changed column heading to N

cdata$botryllid.se<-cdata$botryllid.sd/sqrt(cdata$N)
head(cdata)
#created standard error from SD


cbbPalette3 <- c("#000000", "#FF6600")

[cdata$Invasives=="Present",]
[cdata$Invasives=="Absent",]
[cdata$CO2.Treatment=="Control",]
[cdata$CO2.Treatment=="Elevated CO2",]

windows()
p<-ggplot(cdata[cdata$Invasives=="Present",], aes(x=Week, y=botryllid.mean, linetype=CO2.Treatment, colour=Invasives, shape=CO2.Treatment)) + geom_errorbar(aes(ymin=botryllid.mean-botryllid.se, ymax=botryllid.mean+botryllid.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=4)+ scale_colour_manual(values=cbbPalette3)

p<- p + labs(x = "Time (weeks)", y="Botryllus % cover")+ theme_bw()+ theme(text = element_text(size=20), axis.text = element_text(size=20)) + theme(axis.title.y = element_text(angle=90))

p + scale_linetype_manual(values=c("solid", "dashed"))+ theme(legend.text = element_text(colour="black", size = 20))+theme(legend.title = element_text(colour="black", size=20))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())



### Making lines disappear
cbbPalette4 <- c("#000000", "#FF6600","NA", "NA" )
cbbPalette4 <- c("#000000", "#FF6600","#000000", "#FF6600")

p<-ggplot(cdata, aes(x=Week, y=botryllid.mean, linetype=Treatment, colour=Treatment, shape=CO2.Treatment)) + geom_errorbar(aes(ymin=botryllid.mean-botryllid.se, ymax=botryllid.mean+botryllid.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=4)+ scale_colour_manual(values=cbbPalette4)

p<- p + labs(x = "Time (weeks)", y="Native species richness")+ theme_bw()+ theme(text = element_text(size=20), axis.text = element_text(size=20)) + theme(axis.title.y = element_text(angle=90))

p + scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed"))+ theme(legend.text = element_text(colour="black", size = 20))+theme(legend.title = element_text(colour="black", size=20))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())








cbbPalette4 <- c("#000000", "#FF6600", "#000000", "#FF6600")
p<-ggplot(cdata, aes(x=Week, y=botryllid.mean, linetype=Treatment, colour=Treatment, shape=Treatment)) + geom_errorbar(aes(ymin=botryllid.mean-botryllid.se, ymax=botryllid.mean+botryllid.se), colour="black", width=.3) +geom_line(size=1) + geom_point(size=4)+ scale_colour_manual(values=cbbPalette4)
#rename my plot to p, makes it easier to work with

p<- p + labs(x = "Time (weeks)", y="Proportion botryllid eaten by nudibranchs")+ theme_bw()+ theme(text = element_text(size=20), axis.text = element_text(size=20)) + theme(axis.title.y = element_text(angle=90))

p + scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed"))+ theme(legend.text = element_text(colour="black", size = 20))+theme(legend.title = element_text(colour="black", size=20))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())






