
setwd("/Users/Norah Brown/Documents/Projects/Summer 2015 food availability experiments")

invasion.exp.data<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
#either read in the percent cover file or the with presence file
head(invasion.exp.data)


fix(invasion.exp.data)



invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==16, ]
invasion.exp.data.14<-invasion.exp.data[invasion.exp.data$Week==14, ]
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]
invasion.exp.data.10<-invasion.exp.data[invasion.exp.data$Week==10, ]
invasion.exp.data.8<-invasion.exp.data[invasion.exp.data$Week==8, ]
invasion.exp.data.6<-invasion.exp.data[invasion.exp.data$Week==6, ]
invasion.exp.data.4<-invasion.exp.data[invasion.exp.data$Week==4, ]


fix(invasion.exp.data.16)
library(doBy)
library(plyr)
library(ggplot2) 
library(doBy)
library(grid)


cbbPalette.all.purple<- c( "#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")

cbbPalette.all.purple.purple<-c( "#666666","#CC99CC", "#666666","#CC99CC")

######### 

cdata.botryllid <- summaryBy(botryllid ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.botryllid)[names(cdata.botryllid)=="botryllid.length"]<-"N"
cdata.botryllid$botryllid.se<-cdata.botryllid$botryllid.sd/sqrt(cdata.botryllid$N)
head(cdata.botryllid)

plot.botryllid<-ggplot(cdata.botryllid, aes(x=Week, y=botryllid.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=botryllid.mean-botryllid.se, ymax=botryllid.mean+botryllid.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.botryllid<- plot.botryllid + labs(x = "Time (weeks)", y="botryllid")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.botryllid<-plot.botryllid + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.botryllid <- plot.botryllid + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.botryllid 

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.botryllid.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$botryllid*0.01), colour=Invasives)) + scale_colour_discrete(name = "Invasives")+geom_point(size=5) + guides(fill=FALSE) + scale_colour_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", size=1.5, method="lm")
plot.botryllid.16<- plot.botryllid.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species botryllid"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.botryllid.16<- plot.botryllid.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.botryllid.16<- plot.botryllid.16+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.botryllid.16

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.botryllid.14<- ggplot(invasion.exp.data.14, aes(x=invasion.exp.data.14$min.10.pH, y=(invasion.exp.data.14$botryllid*0.01), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill = "grey80", size=1.5)
plot.botryllid.14<- plot.botryllid.14 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species botryllid"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.botryllid.14<- plot.botryllid.14 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.botryllid.14<- plot.botryllid.14+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.botryllid.14








#bare space --> occupied space

invasion.exp.data$occupied.space<-100-(invasion.exp.data$bare)


cdata.occupied.space <- summaryBy(occupied.space ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.occupied.space)[names(cdata.occupied.space)=="occupied.space.length"]<-"N"
cdata.occupied.space$occupied.space.se<-cdata.occupied.space$occupied.space.sd/sqrt(cdata.occupied.space$N)
head(cdata.occupied.space)

plot.occupied.space<-ggplot(cdata.occupied.space, aes(x=Week, y=occupied.space.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=occupied.space.mean-occupied.space.se, ymax=occupied.space.mean+occupied.space.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.occupied.space<- plot.occupied.space + labs(x = "Time (weeks)", y="Occupied space")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.occupied.space<-plot.occupied.space + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.occupied.space <- plot.occupied.space + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.occupied.space 

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.occupied.space.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$occupied.space*0.01), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill = "grey80", size=1.5)
plot.occupied.space.16<- plot.occupied.space.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species occupied.space"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.occupied.space.16<- plot.occupied.space.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.occupied.space.16<- plot.occupied.space.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.occupied.space.16


##### biomass

plot.tile.brick.biomass.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$tile.brick.biomass), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'lm', size=1.5)
plot.tile.brick.biomass.16<- plot.tile.brick.biomass.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Biomass (g) per brick"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.tile.brick.biomass.16<- plot.tile.brick.biomass.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.tile.brick.biomass.16<- plot.tile.brick.biomass.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.tile.brick.biomass.16


##### biomass
# per unit space ... remember this is wet weight. 
invasion.exp.data.16$stand.biomass<-(invasion.exp.data.16$biomass)/(invasion.exp.data.16$occupied.space)

plot.stand.biomass.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$stand.biomass), colour=Invasives)) + scale_colour_discrete(name = "Invasives") + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80",method = 'glm', method.args = list(family = "Gamma"), size=1.5)
plot.stand.biomass.16<- plot.stand.biomass.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Biomass (g) per 1% cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stand.biomass.16<- plot.stand.biomass.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stand.biomass.16<- plot.stand.biomass.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stand.biomass.16


##### num.species.no.bot redone
 

cdata.num.species.no.bot <- summaryBy(num.species.no.bot ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.num.species.no.bot)[names(cdata.num.species.no.bot)=="num.species.no.bot.length"]<-"N"
cdata.num.species.no.bot$num.species.no.bot.se<-cdata.num.species.no.bot$num.species.no.bot.sd/sqrt(cdata.num.species.no.bot$N)
head(cdata.num.species.no.bot)

plot.num.species.no.bot<-ggplot(cdata.num.species.no.bot, aes(x=Week, y=num.species.no.bot.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.species.no.bot.mean-num.species.no.bot.se, ymax=num.species.no.bot.mean+num.species.no.bot.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.species.no.bot<- plot.num.species.no.bot + labs(x = "Time (weeks)", y="Species num.species.no.bot")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.species.no.bot<-plot.num.species.no.bot + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.species.no.bot <- plot.num.species.no.bot + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.species.no.bot 


plot.num.species.no.bot.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$num.species.no.bot), colour=Invasives))  + geom_point(size=5) + guides(fill=FALSE) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)+ scale_fill_manual(values=cbbPalette.all.purple)+ scale_colour_discrete(name = "Invasives")
plot.num.species.no.bot.16<- plot.num.species.no.bot.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species num.species.no.bot"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.species.no.bot.16<- plot.num.species.no.bot.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.species.no.bot.16<- plot.num.species.no.bot.16+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.species.no.bot.16


plot.num.species.no.bot.12<- ggplot(invasion.exp.data.12, aes(x=invasion.exp.data.12$min.10.pH, y=(invasion.exp.data.12$num.species.no.bot), colour=Invasives)) + scale_colour_discrete(name = "Invasives")+ geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.species.no.bot.12<- plot.num.species.no.bot.12 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species num.species.no.bot"))  + theme(text = element_text(size=12), axis.text = element_text(size=12))+theme(axis.title.y = element_text(angle=90))
plot.num.species.no.bot.12<- plot.num.species.no.bot.12 + theme(legend.text = element_text(colour="black", size = 12))+ theme(legend.title = element_text(colour="black", size=12))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.species.no.bot.12<- plot.num.species.no.bot.12+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.species.no.bot.12



##### shannon.diversity redone

cdata.shannon.diversity <- summaryBy(shannon.diversity ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.shannon.diversity)[names(cdata.shannon.diversity)=="shannon.diversity.length"]<-"N"
cdata.shannon.diversity$shannon.diversity.se<-cdata.shannon.diversity$shannon.diversity.sd/sqrt(cdata.shannon.diversity$N)
head(cdata.shannon.diversity)

plot.shannon.diversity<-ggplot(cdata.shannon.diversity, aes(x=Week, y=shannon.diversity.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=shannon.diversity.mean-shannon.diversity.se, ymax=shannon.diversity.mean+shannon.diversity.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.shannon.diversity<- plot.shannon.diversity + labs(x = "Time (weeks)", y="Shannon diversity (H')")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.shannon.diversity<-plot.shannon.diversity + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.shannon.diversity <- plot.shannon.diversity + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.shannon.diversity 


plot.shannon.diversity.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$shannon.diversity), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "Gamma"), size=1.5)
plot.shannon.diversity.16<- plot.shannon.diversity.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Shannon diversity"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.shannon.diversity.16<- plot.shannon.diversity.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.shannon.diversity.16<- plot.shannon.diversity.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.shannon.diversity.16

plot.shannon.diversity.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$shannon.diversity), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'lm', size=1.5)
plot.shannon.diversity.16<- plot.shannon.diversity.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Shannon diversity"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.shannon.diversity.16<- plot.shannon.diversity.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.shannon.diversity.16<- plot.shannon.diversity.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.shannon.diversity.16

##### biomass mussel
head(invasion.exp.data.16)
plot.Mussel.wet.weight.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=log(invasion.exp.data.16$Mussel.wet.weight+0.1), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'lm', size=1.5)
plot.Mussel.wet.weight.16<- plot.Mussel.wet.weight.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Log Mussel.wet.weight"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.Mussel.wet.weight.16<- plot.Mussel.wet.weight.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.Mussel.wet.weight.16<- plot.Mussel.wet.weight.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.Mussel.wet.weight.16

gamma.16<-fitdistr(invasion.exp.data.16$Mussel.wet.weight + 0.01, "gamma")
qqp(invasion.exp.data.16$Mussel.wet.weight, "gamma", shape = gamma.16$estimate[[1]], rate = gamma.16$estimate[[2]])

qqp(invasion.exp.data.16$Mussel.wet.weight, "norm")
qqp(log(invasion.exp.data.16$Mussel.wet.weight+0.1), "norm")


lm.Mussel.wet.weight<-lm(formula = log(Mussel.wet.weight+0.1) ~ min.10.pH*Invasives, data = invasion.exp.data.16)
summary(lm.Mussel.wet.weight)
plot(lm.Mussel.wet.weight)
anova(lm.Mussel.wet.weight, test = "F")

lm.Mussel.wet.weight<-lm(formula = log(Mussel.wet.weight+0.1) ~ CO2.Treatment*Invasives, data = invasion.exp.data.16)
summary(lm.Mussel.wet.weight)
plot(lm.Mussel.wet.weight)
anova(lm.Mussel.wet.weight, test = "F")

invasion.exp.data.16$Weight.per.ind<-(invasion.exp.data.16$Weight.per.ind+0.01)

### weight per individuals
plot.Weight.per.ind.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=log(invasion.exp.data.16$Weight.per.ind+0.01), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'lm', size=1.5)
plot.Weight.per.ind.16<- plot.Weight.per.ind.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Weight.per.ind"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.Weight.per.ind.16<- plot.Weight.per.ind.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.Weight.per.ind.16<- plot.Weight.per.ind.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.Weight.per.ind.16

gamma.16<-fitdistr(invasion.exp.data.16$Weight.per.ind + 0.01, "gamma")
qqp(invasion.exp.data.16$Weight.per.ind, "gamma", shape = gamma.16$estimate[[1]], rate = gamma.16$estimate[[2]])

qqp(log(invasion.exp.data.16$Weight.per.ind+0.01), "norm")

lm.Weight.per.ind<-lm(formula = log(Weight.per.ind+0.01) ~ min.10.pH*Invasives, data = invasion.exp.data.16)
summary(lm.Weight.per.ind)
plot(lm.Weight.per.ind)
anova(lm.Weight.per.ind, test = "F")

### need to get the actual numbers!! 
plot.nudis.per.tank.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$nudis.per.tank), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.nudis.per.tank.16<- plot.nudis.per.tank.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("# nudibranchs per tank"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.nudis.per.tank.16<- plot.nudis.per.tank.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.nudis.per.tank.16<- plot.nudis.per.tank.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.nudis.per.tank.16



poisson.16<-fitdistr(invasion.exp.data.16$nudis.per.tank, "Poisson")
qqp(invasion.exp.data.16$nudis.per.tank, "pois", poisson.16$estimate)

glm.poisson.nudis.per.tank.16<-glm(formula = (nudis.per.tank) ~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "poisson")
summary(glm.poisson.nudis.per.tank.16)
plot(glm.poisson.nudis.per.tank.16)
anova(glm.poisson.nudis.per.tank.16, test = "Chi")

##### Didemnum
cdata.didemnum <- summaryBy(didemnum ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.didemnum)[names(cdata.didemnum)=="didemnum.length"]<-"N"
cdata.didemnum$didemnum.se<-cdata.didemnum$didemnum.sd/sqrt(cdata.didemnum$N)
head(cdata.didemnum)

plot.didemnum<-ggplot(cdata.didemnum, aes(x=Week, y=didemnum.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=didemnum.mean-didemnum.se, ymax=didemnum.mean+didemnum.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.didemnum<- plot.didemnum + labs(x = "Time (weeks)", y="Species didemnum")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.didemnum<-plot.didemnum + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.didemnum <- plot.didemnum + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.didemnum 

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.didemnum.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$didemnum*0.01), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill = "grey80", size=1.5)
plot.didemnum.16<- plot.didemnum.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species didemnum"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.didemnum.16<- plot.didemnum.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.didemnum.16<- plot.didemnum.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.didemnum.16


#Didemnum #
cdata.num.didemnum <- summaryBy(num.didemnum ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.num.didemnum)[names(cdata.num.didemnum)=="num.didemnum.length"]<-"N"
cdata.num.didemnum$num.didemnum.se<-cdata.num.didemnum$num.didemnum.sd/sqrt(cdata.num.didemnum$N)
head(cdata.num.didemnum)

plot.num.didemnum<-ggplot(cdata.num.didemnum, aes(x=Week, y=num.didemnum.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.didemnum.mean-num.didemnum.se, ymax=num.didemnum.mean+num.didemnum.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.didemnum<- plot.num.didemnum + labs(x = "Time (weeks)", y="Species num.didemnum")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.didemnum<-plot.num.didemnum + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.didemnum <- plot.num.didemnum + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.didemnum 


plot.num.didemnum.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$num.didemnum), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.didemnum.16<- plot.num.didemnum.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species num.didemnum"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.didemnum.16<- plot.num.didemnum.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.didemnum.16<- plot.num.didemnum.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.didemnum.16


#### nudi eggs
head(invasion.exp.data.16)
invasion.exp.data.16$nudi.eggs.per.tank.total<-(invasion.exp.data.16$nudi.eggs.per.tank)+(invasion.exp.data.16$num.nudi.eggs)


plot.nudi.eggs.per.tank.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$nudi.eggs.per.tank), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.nudi.eggs.per.tank.16<- plot.nudi.eggs.per.tank.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("# nudibranch eggs laid per tank"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.nudi.eggs.per.tank.16<- plot.nudi.eggs.per.tank.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.nudi.eggs.per.tank.16<- plot.nudi.eggs.per.tank.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.nudi.eggs.per.tank.16


poisson.16<-fitdistr(invasion.exp.data.16$nudi.eggs.per.tank, "Poisson")
qqp(invasion.exp.data.16$nudi.eggs.per.tank, "pois", poisson.16$estimate)

glm.poisson.nudi.eggs.per.tank.16<-glm(formula = (nudi.eggs.per.tank) ~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "poisson")
summary(glm.poisson.nudi.eggs.per.tank.16)
plot(glm.poisson.nudi.eggs.per.tank.16)
anova(glm.poisson.nudi.eggs.per.tank.16, test = "Chi")


plot.num.nudi.eggs.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$num.nudi.eggs), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.nudi.eggs.16<- plot.num.nudi.eggs.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("# nudibranch eggs laid  per tank"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.nudi.eggs.16<- plot.num.nudi.eggs.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.nudi.eggs.16<- plot.num.nudi.eggs.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.nudi.eggs.16

plot.num.nudi.eggs<-ggplot(cdata.num.nudi.eggs, aes(x=Week, y=num.nudi.eggs.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.nudi.eggs.mean-num.nudi.eggs.se, ymax=num.nudi.eggs.mean+num.nudi.eggs.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.nudi.eggs<- plot.num.nudi.eggs + labs(x = "Time (weeks)", y="Species num.nudi.eggs")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.nudi.eggs<-plot.num.nudi.eggs + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.nudi.eggs <- plot.num.nudi.eggs + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.nudi.eggs 

head(invasion.exp.data.16)

### Disporella total 

cdata.disporella <- summaryBy(disporella ~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.disporella)[names(cdata.disporella)=="disporella.length"]<-"N"
cdata.disporella$disporella.se<-cdata.disporella$disporella.sd/sqrt(cdata.disporella$N)
head(cdata.disporella)

plot.disporella<-ggplot(cdata.disporella, aes(x=Week, y=disporella.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.mean-disporella.se, ymax=disporella.mean+disporella.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella<- plot.disporella + labs(x = "Time (weeks)", y="disporella")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella<-plot.disporella + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella <- plot.disporella + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella 


binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=(invasion.exp.data.16$disporella*0.01), colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.16<- plot.disporella.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))+ylim(0,0.25)
plot.disporella.16<- plot.disporella.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.16<- plot.disporella.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.16



#### Disporella flat 
cdata.disporella.flat <- summaryBy(disporella.flat~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.disporella.flat)[names(cdata.disporella.flat)=="disporella.flat.length"]<-"N"
cdata.disporella.flat$disporella.flat.se<-cdata.disporella.flat$disporella.flat.sd/sqrt(cdata.disporella.flat$N)
head(cdata.disporella.flat)

plot.disporella.flat<-ggplot(cdata.disporella.flat, aes(x=Week, y=disporella.flat.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.flat.mean-disporella.flat.se, ymax=disporella.flat.mean+disporella.flat.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella.flat<- plot.disporella.flat + labs(x = "Time (weeks)", y="disporella.flat")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella.flat<-plot.disporella.flat + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella.flat <- plot.disporella.flat + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.flat

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.flat.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$disporella.flat, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.flat.16<- plot.disporella.flat.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.flat.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.disporella.flat.16<- plot.disporella.flat.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.flat.16<- plot.disporella.flat.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))

plot.disporella.flat.16



###Disporella erect.all

cdata.disporella.erect.all <- summaryBy(disporella.erect.all~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.disporella.erect.all)[names(cdata.disporella.erect.all)=="disporella.erect.all.length"]<-"N"
cdata.disporella.erect.all$disporella.erect.all.se<-cdata.disporella.erect.all$disporella.erect.all.sd/sqrt(cdata.disporella.erect.all$N)
head(cdata.disporella.erect.all)

plot.disporella.erect.all<-ggplot(cdata.disporella.erect.all, aes(x=Week, y=disporella.erect.all.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.erect.all.mean-disporella.erect.all.se, ymax=disporella.erect.all.mean+disporella.erect.all.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella.erect.all<- plot.disporella.erect.all + labs(x = "Time (weeks)", y="disporella.erect.all")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.all<-plot.disporella.erect.all + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella.erect.all <- plot.disporella.erect.all + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.all

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.erect.all.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$disporella.erect.all, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.erect.all.16<- plot.disporella.erect.all.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.erect.all.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.all.16<- plot.disporella.erect.all.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.erect.all.16<- plot.disporella.erect.all.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.all.16


#### Disporella with knobs
cdata.disporella.erect.knobs <- summaryBy(disporella.erect.knobs~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.disporella.erect.knobs)[names(cdata.disporella.erect.knobs)=="disporella.erect.knobs.length"]<-"N"
cdata.disporella.erect.knobs$disporella.erect.knobs.se<-cdata.disporella.erect.knobs$disporella.erect.knobs.sd/sqrt(cdata.disporella.erect.knobs$N)
head(cdata.disporella.erect.knobs)

plot.disporella.erect.knobs<-ggplot(cdata.disporella.erect.knobs, aes(x=Week, y=disporella.erect.knobs.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.erect.knobs.mean-disporella.erect.knobs.se, ymax=disporella.erect.knobs.mean+disporella.erect.knobs.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella.erect.knobs<- plot.disporella.erect.knobs + labs(x = "Time (weeks)", y="disporella.erect.knobs")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.knobs<-plot.disporella.erect.knobs + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella.erect.knobs <- plot.disporella.erect.knobs + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.knobs

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.erect.knobs.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$disporella.erect.knobs, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.erect.knobs.16<- plot.disporella.erect.knobs.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.erect.knobs.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.knobs.16<- plot.disporella.erect.knobs.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.erect.knobs.16<- plot.disporella.erect.knobs.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.knobs.16


#### Disporella with fan
cdata.disporella.erect.fan <- summaryBy(disporella.erect.fan~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.disporella.erect.fan)[names(cdata.disporella.erect.fan)=="disporella.erect.fan.length"]<-"N"
cdata.disporella.erect.fan$disporella.erect.fan.se<-cdata.disporella.erect.fan$disporella.erect.fan.sd/sqrt(cdata.disporella.erect.fan$N)
head(cdata.disporella.erect.fan)

plot.disporella.erect.fan<-ggplot(cdata.disporella.erect.fan, aes(x=Week, y=disporella.erect.fan.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.erect.fan.mean-disporella.erect.fan.se, ymax=disporella.erect.fan.mean+disporella.erect.fan.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella.erect.fan<- plot.disporella.erect.fan + labs(x = "Time (weeks)", y="disporella.erect.fan")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.fan<-plot.disporella.erect.fan + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella.erect.fan <- plot.disporella.erect.fan + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.fan

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.erect.fan.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$disporella.erect.fan, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.erect.fan.16<- plot.disporella.erect.fan.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.erect.fan.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.fan.16<- plot.disporella.erect.fan.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.erect.fan.16<- plot.disporella.erect.fan.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.fan.16

#### Disporella ruffles
cdata.disporella.erect.ruffles <- summaryBy(disporella.erect.ruffles~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.disporella.erect.ruffles)[names(cdata.disporella.erect.ruffles)=="disporella.erect.ruffles.length"]<-"N"
cdata.disporella.erect.ruffles$disporella.erect.ruffles.se<-cdata.disporella.erect.ruffles$disporella.erect.ruffles.sd/sqrt(cdata.disporella.erect.ruffles$N)
head(cdata.disporella.erect.ruffles)

plot.disporella.erect.ruffles<-ggplot(cdata.disporella.erect.ruffles, aes(x=Week, y=disporella.erect.ruffles.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=disporella.erect.ruffles.mean-disporella.erect.ruffles.se, ymax=disporella.erect.ruffles.mean+disporella.erect.ruffles.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.disporella.erect.ruffles<- plot.disporella.erect.ruffles + labs(x = "Time (weeks)", y="disporella.erect.ruffles")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.ruffles<-plot.disporella.erect.ruffles + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.disporella.erect.ruffles <- plot.disporella.erect.ruffles + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.ruffles

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.disporella.erect.ruffles.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$disporella.erect.ruffles, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.disporella.erect.ruffles.16<- plot.disporella.erect.ruffles.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("disporella.erect.ruffles.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.disporella.erect.ruffles.16<- plot.disporella.erect.ruffles.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.disporella.erect.ruffles.16<- plot.disporella.erect.ruffles.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disporella.erect.ruffles.16



######## COUNTS

library(dplyr)
install.packages("tidyr")
library(tidyr)
install.packages("reshape2")
library(reshape2)
invasion.exp.data.16$Treatment<-as.factor(invasion.exp.data.16$Treatment)

is.factor(invasion.exp.data.16$Treatment)
levels(invasion.exp.data.16$Treatment)

which( colnames(invasion.exp.data.16)=="num.disporella" )

long.morph<- gather(invasion.exp.data.16, morphology, count,num.disporella.flat, num.disporella.erect.knobs, num.disporella.erect.fan, num.disporella.erect.ruffles)
head(long.morph)

ggplot(data = long.morph, aes(x = Treatment, y =count , fill = morphology)) + 
  geom_bar(stat="identity")+ scale_fill_brewer(palette = "Reds")

long.morph$morphology<- factor(long.morph$morphology, as.character(long.morph$morphology))

lvl4<-c("Flat", "Erect with knobs", "Erect with fan", "Erect with ruffles")
ce = ddply(long.morph, "Treatment", mutate, percent_total = count/sum(count) * 100)
ggplot(ce, aes(x=Treatment, y=percent_total, fill=morphology)) + 
  geom_bar(stat='identity') + scale_fill_brewer(palette="Reds", labels=lvl4)



#### Disporella counts all 
cdata.num.disporella <- summaryBy(num.disporella~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella)[names(cdata.num.disporella)=="num.disporella.length"]<-"N"
cdata.num.disporella$num.disporella.se<-cdata.num.disporella$num.disporella.sd/sqrt(cdata.num.disporella$N)
head(cdata.num.disporella)

plot.num.disporella<-ggplot(cdata.num.disporella, aes(x=Week, y=num.disporella.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.mean-num.disporella.se, ymax=num.disporella.mean+num.disporella.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella<- plot.num.disporella + labs(x = "Time (weeks)", y="num.disporella")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella<-plot.num.disporella + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella <- plot.num.disporella + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella

plot.num.disporella.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.16<- plot.num.disporella.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.16<- plot.num.disporella.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.16<- plot.num.disporella.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.16



###### Number disporella flat
cdata.num.disporella.flat <- summaryBy(num.disporella.flat~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella.flat)[names(cdata.num.disporella.flat)=="num.disporella.flat.length"]<-"N"
cdata.num.disporella.flat$num.disporella.flat.se<-cdata.num.disporella.flat$num.disporella.flat.sd/sqrt(cdata.num.disporella.flat$N)
head(cdata.num.disporella.flat)

plot.num.disporella.flat<-ggplot(cdata.num.disporella.flat, aes(x=Week, y=num.disporella.flat.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.flat.mean-num.disporella.flat.se, ymax=num.disporella.flat.mean+num.disporella.flat.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella.flat<- plot.num.disporella.flat + labs(x = "Time (weeks)", y="num.disporella.flat")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella.flat<-plot.num.disporella.flat + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella.flat <- plot.num.disporella.flat + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.flat

plot.num.disporella.flat.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella.flat, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.flat.16<- plot.num.disporella.flat.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella.flat"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.flat.16<- plot.num.disporella.flat.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.flat.16<- plot.num.disporella.flat.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.flat.16

#### Disporella erect.all
cdata.num.disporella.erect.all <- summaryBy(num.disporella.erect.all~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella.erect.all)[names(cdata.num.disporella.erect.all)=="num.disporella.erect.all.length"]<-"N"
cdata.num.disporella.erect.all$num.disporella.erect.all.se<-cdata.num.disporella.erect.all$num.disporella.erect.all.sd/sqrt(cdata.num.disporella.erect.all$N)
head(cdata.num.disporella.erect.all)

plot.num.disporella.erect.all<-ggplot(cdata.num.disporella.erect.all, aes(x=Week, y=num.disporella.erect.all.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.erect.all.mean-num.disporella.erect.all.se, ymax=num.disporella.erect.all.mean+num.disporella.erect.all.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella.erect.all<- plot.num.disporella.erect.all + labs(x = "Time (weeks)", y="num.disporella.erect.all")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.all<-plot.num.disporella.erect.all + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella.erect.all <- plot.num.disporella.erect.all + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.all

plot.num.disporella.erect.all.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella.erect.all, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.erect.all.16<- plot.num.disporella.erect.all.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella.erect.all"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.all.16<- plot.num.disporella.erect.all.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.erect.all.16<- plot.num.disporella.erect.all.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.all.16


#### Disporella erect knobs 
cdata.num.disporella.erect.knobs <- summaryBy(num.disporella.erect.knobs~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella.erect.knobs)[names(cdata.num.disporella.erect.knobs)=="num.disporella.erect.knobs.length"]<-"N"
cdata.num.disporella.erect.knobs$num.disporella.erect.knobs.se<-cdata.num.disporella.erect.knobs$num.disporella.erect.knobs.sd/sqrt(cdata.num.disporella.erect.knobs$N)
head(cdata.num.disporella.erect.knobs)

plot.num.disporella.erect.knobs<-ggplot(cdata.num.disporella.erect.knobs, aes(x=Week, y=num.disporella.erect.knobs.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.erect.knobs.mean-num.disporella.erect.knobs.se, ymax=num.disporella.erect.knobs.mean+num.disporella.erect.knobs.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella.erect.knobs<- plot.num.disporella.erect.knobs + labs(x = "Time (weeks)", y="num.disporella.erect.knobs")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.knobs<-plot.num.disporella.erect.knobs + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella.erect.knobs <- plot.num.disporella.erect.knobs + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.knobs

plot.num.disporella.erect.knobs.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella.erect.knobs, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.erect.knobs.16<- plot.num.disporella.erect.knobs.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella.erect.knobs"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.knobs.16<- plot.num.disporella.erect.knobs.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.erect.knobs.16<- plot.num.disporella.erect.knobs.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.knobs.16



#### Number disporella Fan
cdata.num.disporella.erect.fan <- summaryBy(num.disporella.erect.fan~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella.erect.fan)[names(cdata.num.disporella.erect.fan)=="num.disporella.erect.fan.length"]<-"N"
cdata.num.disporella.erect.fan$num.disporella.erect.fan.se<-cdata.num.disporella.erect.fan$num.disporella.erect.fan.sd/sqrt(cdata.num.disporella.erect.fan$N)
head(cdata.num.disporella.erect.fan)

plot.num.disporella.erect.fan<-ggplot(cdata.num.disporella.erect.fan, aes(x=Week, y=num.disporella.erect.fan.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.erect.fan.mean-num.disporella.erect.fan.se, ymax=num.disporella.erect.fan.mean+num.disporella.erect.fan.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella.erect.fan<- plot.num.disporella.erect.fan + labs(x = "Time (weeks)", y="num.disporella.erect.fan")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.fan<-plot.num.disporella.erect.fan + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella.erect.fan <- plot.num.disporella.erect.fan + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.fan

plot.num.disporella.erect.fan.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella.erect.fan, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.erect.fan.16<- plot.num.disporella.erect.fan.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella.erect.fan"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.fan.16<- plot.num.disporella.erect.fan.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.erect.fan.16<- plot.num.disporella.erect.fan.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.fan.16



#### Number Disporella ruffles
cdata.num.disporella.erect.ruffles <- summaryBy(num.disporella.erect.ruffles~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) })
names(cdata.num.disporella.erect.ruffles)[names(cdata.num.disporella.erect.ruffles)=="num.disporella.erect.ruffles.length"]<-"N"
cdata.num.disporella.erect.ruffles$num.disporella.erect.ruffles.se<-cdata.num.disporella.erect.ruffles$num.disporella.erect.ruffles.sd/sqrt(cdata.num.disporella.erect.ruffles$N)
head(cdata.num.disporella.erect.ruffles)

plot.num.disporella.erect.ruffles<-ggplot(cdata.num.disporella.erect.ruffles, aes(x=Week, y=num.disporella.erect.ruffles.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.disporella.erect.ruffles.mean-num.disporella.erect.ruffles.se, ymax=num.disporella.erect.ruffles.mean+num.disporella.erect.ruffles.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.disporella.erect.ruffles<- plot.num.disporella.erect.ruffles + labs(x = "Time (weeks)", y="num.disporella.erect.ruffles")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.ruffles<-plot.num.disporella.erect.ruffles + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.disporella.erect.ruffles <- plot.num.disporella.erect.ruffles + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.ruffles

plot.num.disporella.erect.ruffles.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.disporella.erect.ruffles, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', method.args = list(family = "poisson"), size=1.5)
plot.num.disporella.erect.ruffles.16<- plot.num.disporella.erect.ruffles.16 + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.disporella.erect.ruffles"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.disporella.erect.ruffles.16<- plot.num.disporella.erect.ruffles.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.disporella.erect.ruffles.16<- plot.num.disporella.erect.ruffles.16 + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.disporella.erect.ruffles.16


#########
#Proportion disporella

#### Proportion of the disporella counts that are erect fan
invasion.exp.data$num.prop.disporella.erect.fan<-(invasion.exp.data$num.disporella.erect.fan)/(invasion.exp.data$num.disporella)
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]

length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.num.prop.disporella.erect.fan <- summaryBy(num.prop.disporella.erect.fan~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=c(length2,mean,sd), na.rm=TRUE)
names(cdata.num.prop.disporella.erect.fan)[names(cdata.num.prop.disporella.erect.fan)=="num.prop.disporella.erect.fan.length2"]<-"N"
cdata.num.prop.disporella.erect.fan$num.prop.disporella.erect.fan.se<-cdata.num.prop.disporella.erect.fan$num.prop.disporella.erect.fan.sd/sqrt(cdata.num.prop.disporella.erect.fan$N)
head(cdata.num.prop.disporella.erect.fan)

plot.num.prop.disporella.erect.fan<-ggplot(cdata.num.prop.disporella.erect.fan, aes(x=Week, y=num.prop.disporella.erect.fan.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.prop.disporella.erect.fan.mean-num.prop.disporella.erect.fan.se, ymax=num.prop.disporella.erect.fan.mean+num.prop.disporella.erect.fan.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.prop.disporella.erect.fan<- plot.num.prop.disporella.erect.fan + labs(x = "Time (weeks)", y="num.prop.disporella.erect.fan")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.fan<-plot.num.prop.disporella.erect.fan + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.prop.disporella.erect.fan <- plot.num.prop.disporella.erect.fan + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.fan

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.num.prop.disporella.erect.fan.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect.fan, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.num.prop.disporella.erect.fan.16<- plot.num.prop.disporella.erect.fan.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.erect.fan.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.fan.16<- plot.num.prop.disporella.erect.fan.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect.fan.16<- plot.num.prop.disporella.erect.fan.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.fan.16



#### Proportion of the disporella counts that are erect in total ### can't totally figure out na.rm
?summaryBy
invasion.exp.data$num.prop.disporella.erect<-(invasion.exp.data$num.disporella.erect.all)/(invasion.exp.data$num.disporella)
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]
 
cdata.num.prop.disporella.erect <- summaryBy(num.prop.disporella.erect~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data,FUN=function(x) { c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE), length=length2(x))})
names(cdata.num.prop.disporella.erect)[names(cdata.num.prop.disporella.erect)=="num.prop.disporella.erect.length"]<-"N"
cdata.num.prop.disporella.erect$num.prop.disporella.erect.se<-cdata.num.prop.disporella.erect$num.prop.disporella.erect.sd/sqrt(cdata.num.prop.disporella.erect$N)
head(cdata.num.prop.disporella.erect)

plot.num.prop.disporella.erect<-ggplot(cdata.num.prop.disporella.erect, aes(x=Week, y=num.prop.disporella.erect.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.prop.disporella.erect.mean-num.prop.disporella.erect.se, ymax=num.prop.disporella.erect.mean+num.prop.disporella.erect.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + labs(x = "Time (weeks)", y="num.prop.disporella.erect")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect<-plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.prop.disporella.erect <- plot.num.prop.disporella.erect + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.num.prop.disporella.erect.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.num.prop.disporella.erect.16<- plot.num.prop.disporella.erect.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.erect.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.16<- plot.num.prop.disporella.erect.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect.16<- plot.num.prop.disporella.erect.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.16

invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]


#### Proportion of the disporella counts that are knobs shaped
invasion.exp.data$num.prop.disporella.erect.knobs<-(invasion.exp.data$num.disporella.erect.knobs)/(invasion.exp.data$num.disporella)
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]

length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.num.prop.disporella.erect.knobs <- summaryBy(num.prop.disporella.erect.knobs~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=c(length2,mean,sd), na.rm=TRUE)
names(cdata.num.prop.disporella.erect.knobs)[names(cdata.num.prop.disporella.erect.knobs)=="num.prop.disporella.erect.knobs.length2"]<-"N"
cdata.num.prop.disporella.erect.knobs$num.prop.disporella.erect.knobs.se<-cdata.num.prop.disporella.erect.knobs$num.prop.disporella.erect.knobs.sd/sqrt(cdata.num.prop.disporella.erect.knobs$N)
head(cdata.num.prop.disporella.erect.knobs)

plot.num.prop.disporella.erect.knobs<-ggplot(cdata.num.prop.disporella.erect.knobs, aes(x=Week, y=num.prop.disporella.erect.knobs.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.prop.disporella.erect.knobs.mean-num.prop.disporella.erect.knobs.se, ymax=num.prop.disporella.erect.knobs.mean+num.prop.disporella.erect.knobs.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.prop.disporella.erect.knobs<- plot.num.prop.disporella.erect.knobs + labs(x = "Time (weeks)", y="num.prop.disporella.erect.knobs")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.knobs<-plot.num.prop.disporella.erect.knobs + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.prop.disporella.erect.knobs <- plot.num.prop.disporella.erect.knobs + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.knobs

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.num.prop.disporella.erect.knobs.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect.knobs, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.num.prop.disporella.erect.knobs.16<- plot.num.prop.disporella.erect.knobs.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.erect.knobs.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.knobs.16<- plot.num.prop.disporella.erect.knobs.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect.knobs.16<- plot.num.prop.disporella.erect.knobs.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.knobs.16


#### Proportion of the disporella counts that are erect ruffles
invasion.exp.data$num.prop.disporella.erect.ruffles<-(invasion.exp.data$num.disporella.erect.ruffles)/(invasion.exp.data$num.disporella)
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]

length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.num.prop.disporella.erect.ruffles <- summaryBy(num.prop.disporella.erect.ruffles~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=c(length2,mean,sd), na.rm=TRUE)
names(cdata.num.prop.disporella.erect.ruffles)[names(cdata.num.prop.disporella.erect.ruffles)=="num.prop.disporella.erect.ruffles.length2"]<-"N"
cdata.num.prop.disporella.erect.ruffles$num.prop.disporella.erect.ruffles.se<-cdata.num.prop.disporella.erect.ruffles$num.prop.disporella.erect.ruffles.sd/sqrt(cdata.num.prop.disporella.erect.ruffles$N)
head(cdata.num.prop.disporella.erect.ruffles)

plot.num.prop.disporella.erect.ruffles<-ggplot(cdata.num.prop.disporella.erect.ruffles, aes(x=Week, y=num.prop.disporella.erect.ruffles.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.prop.disporella.erect.ruffles.mean-num.prop.disporella.erect.ruffles.se, ymax=num.prop.disporella.erect.ruffles.mean+num.prop.disporella.erect.ruffles.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.prop.disporella.erect.ruffles<- plot.num.prop.disporella.erect.ruffles + labs(x = "Time (weeks)", y="num.prop.disporella.erect.ruffles")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.ruffles<-plot.num.prop.disporella.erect.ruffles + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.prop.disporella.erect.ruffles <- plot.num.prop.disporella.erect.ruffles + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.ruffles

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.num.prop.disporella.erect.ruffles.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect.ruffles, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.num.prop.disporella.erect.ruffles.16<- plot.num.prop.disporella.erect.ruffles.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.erect.ruffles.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.erect.ruffles.16<- plot.num.prop.disporella.erect.ruffles.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect.ruffles.16<- plot.num.prop.disporella.erect.ruffles.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect.ruffles.16



#### Proportion of the disporella counts that are flat
invasion.exp.data$num.prop.disporella.flat<-(invasion.exp.data$num.disporella.flat)/(invasion.exp.data$num.disporella)
invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==12, ]

length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.num.prop.disporella.flat <- summaryBy(num.prop.disporella.flat~ Treatment + Invasives + CO2.Treatment + Week, data=invasion.exp.data, FUN=c(length2,mean,sd), na.rm=TRUE)
names(cdata.num.prop.disporella.flat)[names(cdata.num.prop.disporella.flat)=="num.prop.disporella.flat.length2"]<-"N"
cdata.num.prop.disporella.flat$num.prop.disporella.flat.se<-cdata.num.prop.disporella.flat$num.prop.disporella.flat.sd/sqrt(cdata.num.prop.disporella.flat$N)
head(cdata.num.prop.disporella.flat)

plot.num.prop.disporella.flat<-ggplot(cdata.num.prop.disporella.flat, aes(x=Week, y=num.prop.disporella.flat.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=num.prop.disporella.flat.mean-num.prop.disporella.flat.se, ymax=num.prop.disporella.flat.mean+num.prop.disporella.flat.se), width=0.25) + geom_line(size=1) + geom_point(size=5) + scale_colour_manual(values=cbbPalette.all.purple)
plot.num.prop.disporella.flat<- plot.num.prop.disporella.flat + labs(x = "Time (weeks)", y="num.prop.disporella.flat")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.flat<-plot.num.prop.disporella.flat + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num.prop.disporella.flat <- plot.num.prop.disporella.flat + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.flat

binomial_smooth <- function(...) {geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}
plot.num.prop.disporella.flat.16<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.flat, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)
plot.num.prop.disporella.flat.16<- plot.num.prop.disporella.flat.16 + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.flat.16"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num.prop.disporella.flat.16<- plot.num.prop.disporella.flat.16 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.flat.16<- plot.num.prop.disporella.flat.16+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.flat.16

#binomial


####binomial week 8
plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.8, aes(x=invasion.exp.data.8$min.10.pH, y=invasion.exp.data.8$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ binomial_smooth(fill="grey80", size=1.5)

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Proportion dead barnacles"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect



#for counts

plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', family = 'poisson', size=1.5)

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("num.prop.disporella.erect"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect

#for counts week 8

plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.8, aes(x=invasion.exp.data.8$min.10.pH, y=invasion.exp.data.8$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', family = 'poisson', size=1.5)

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("# eaten membranipora"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect

#for counts week 10

plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.10, aes(x=invasion.exp.data.10$min.10.pH, y=invasion.exp.data.10$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', family = 'poisson', size=1.5)

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("# dead barnacles"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect



### week 8
### Take out an outlier:
invasion.exp.data.8.no.outlier<-invasion.exp.data.8[invasion.exp.data.8$num.prop.disporella.erect<20, ]

head(invasion.exp.data.8.no.outlier)

plot.num.prop.disporella.erect.no.outlier<- ggplot(invasion.exp.data.8.no.outlier, aes(x=invasion.exp.data.8.no.outlier$min.10.pH, y=invasion.exp.data.8.no.outlier$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple) + geom_smooth(fill = "grey80", method = 'glm', family = 'poisson', size=1.5)

plot.num.prop.disporella.erect.no.outlier<- plot.num.prop.disporella.erect.no.outlier + theme_bw() + xlab(bquote('pH 10th percentile')) + ylab(expression("# live barnacles"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect.no.outlier<- plot.num.prop.disporella.erect.no.outlier + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect.no.outlier


###for normal
plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill="grey80", size=1.5, method='lm')

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Shannon diversity"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect



### can do glm gamma
plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect, colour=Invasives)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all.purple)+ geom_smooth(fill="grey80", size=1.5, method='glm', family='Gamma')

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Shannon diversity"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect




#####
#plot for presentations

plot.num.prop.disporella.erect<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$min.10.pH, y=invasion.exp.data.16$num.prop.disporella.erect, colour=factor(invasion.exp.data.16$Invasives, levels = c("High", "Low", "high")))) + geom_point(size=4, alpha = 6.5/8) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette6)+ geom_smooth(method = 'glm', family = 'Gamma')

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("num.prop.disporella.erect"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num.prop.disporella.erect<- plot.num.prop.disporella.erect+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num.prop.disporella.erect





################
#Stats March 2016

#### do I need glmm if I don't have random effects



library(glmmADMB)
require(car)
require(MASS)
library(betareg)
library(lmtest)
library(fitdistrplus)

head(invasion.exp.data)
invasion.exp.data.16$Invasives<-as.factor(invasion.exp.data$Invasives)
invasion.exp.data.16$CO2.Treatment<-as.factor(invasion.exp.data$CO2.Treatment)

invasion.exp.data.16$Invasives <- factor(invasion.exp.data.16$Invasives, levels = c("None", "Low", "High"))


invasion.exp.data.8$Invasives <- factor(invasion.exp.data.8$Invasives, levels = c("None", "Low", "High"))
invasion.exp.data.16$Invasives <- factor(invasion.exp.data.16$Invasives, levels = c("High", "Low", "None"))


#########################
invasion.exp.data<-invasion.exp.data[,-5]
# Week 12
head(invasion.exp.data)
fix(invasion.exp.data)

beta.16<-fitdist((0.01*invasion.exp.data.16$num.prop.disporella.erect), "beta", start=NULL)
qqp(invasion.exp.data$num.prop.disporella.erect, "beta", shape1 = beta.16$estimate[[1]], shape2 = beta.16$estimate[[2]])

binom.16<-fitdist((0.01*invasion.exp.data.16$num.prop.disporella.erect), "binom", start=0.1)
qqp(rbinom(invasion.exp.data$num.prop.disporella.erect), size=0.1)




###Don't need a mixed model b/c don't have random effects (i.e. one rep per mesocosm)

?fitdist
?qqp
?glm

library(betareg)
?betareg
full <- c(invasion.exp.data$min.10.pH, invasion.exp.data$Invasives)


###I think betareg can only be used for numeric predictors
#betareg.num.prop.disporella.erect.16<- betareg(formula = 0.01*num.prop.disporella.erect ~ full$min.10.pH, data = invasion.exp.data)
#can't do a glm with a beta distributions.... 
#glm.num.prop.disporella.erect.16<- glm(formula = 0.01*num.prop.disporella.erect ~ min.10.pH*Invasives, data = invasion.exp.data, family = "beta")

invasion.exp.data.16$num.prop.disporella.erect<-(invasion.exp.data.16$num.prop.disporella.erect)+0.0001
glm.binomial.num.prop.disporella.erect.16<- glm(formula = cbind(num.prop.disporella.erect, 1.0001-num.prop.disporella.erect)~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.16)
summary(glm.binomial.num.prop.disporella.erect.16)
anova(glm.binomial.num.prop.disporella.erect.16, test = "Chi")


head(invasion.exp.data.16.na.omit)
invasion.exp.data.16.na.omit<-na.omit(invasion.exp.data.16)
glm.binomial.num.prop.disporella.erect.16.na.omit<- glm(formula = cbind(num.prop.disporella.erect, 1.0001-num.prop.disporella.erect)~ min.10.pH*Invasives, data = invasion.exp.data.16.na.omit, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.16.na.omit)
summary(glm.binomial.num.prop.disporella.erect.16.na.omit)
anova(glm.binomial.num.prop.disporella.erect.16.na.omit, test = "Chi")


#### week 8 binomial
invasion.exp.data.8$num.prop.disporella.erect<-(invasion.exp.data.8$num.prop.disporella.erect)+0.0001
glm.binomial.num.prop.disporella.erect.8<- glm(formula = cbind(num.prop.disporella.erect, 1.0001-num.prop.disporella.erect)~ min.10.pH*Invasives, data = invasion.exp.data.8, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.8)
summary(glm.binomial.num.prop.disporella.erect.8)
anova(glm.binomial.num.prop.disporella.erect.8, test = "Chi")



head(invasion.exp.data.16)
##need to change from 100 to 1 for proportion dead. 
#, contrasts = c("contr.sum", "contr.poly")
#library(car)
#Anova(glm.bionomial.num.prop.disporella.erect.16, type = 3)


###visreg - visualize the fit of the model to the data
#You can use visreg to visualize the model fit on 
#the transformed scale (the function uses predict(z) to generate the result). 
#Replace "x" in the command below with the name of your x-variable (in quotes). 
#The glm method fits a linear model on the transformed scale, and this is what you 
#will visualize. The dots are not the transformed data, however. They are 
#"working values" obtained by transforming the residuals of the fitted model 
#on the original scale. glm repeatedly recalculates the working values and the 
#fitted model as it converges on the maximum likelihood estimate. visreg shows you 
#the results from the final iteration.

library(visreg)
visreg(glm.binomial.num.prop.disporella.erect.8, xvar = "min.10.pH", by= "Invasives", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))




##############################################################
## poisson for counts and num species:
library(glmmADMB)
require(car)
require(MASS)
library(poissonreg)
library(lmtest)
library(fitdistrplus)

head(invasion.exp.data.16)

?glm
###PS it's not the link function but the family that determines whether you should use Chi-square or F 
#(specifically, whether the scale parameter is fixed [Poisson, binomial] or 
#estimated [Gaussian, Gamma, quasi-likelihood fits] 
#	(1) Yes, chi-square for fixed and F for estimated.
head(invasion.exp.data.16)

nbinom12 <- fitdistr(invasion.exp.data.16$num.disporella.erect.all, "Negative Binomial")
qqp(invasion.exp.data.16$num.disporella.erect.all, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])


poisson.16<-fitdistr(invasion.exp.data.16$num.disporella.erect.all, "Poisson")
qqp(invasion.exp.data.16$num.disporella.erect.all, "pois", poisson.16$estimate)

glm.poisson.num.disporella.erect.all.16<-glm(formula = (num.disporella.erect.all) ~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "poisson")
summary(glm.poisson.num.disporella.erect.all.16)
plot(glm.poisson.num.disporella.erect.all.16)
anova(glm.poisson.num.disporella.erect.all.16, test = "Chi")

library(visreg)
visreg(glm.poisson.num.prop.disporella.erect.16, xvar = "min.10.pH", by= "Invasives", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))



##### Week 8
glm.poisson.num.prop.disporella.erect.8<-glm(formula = (num.prop.disporella.erect) ~ min.10.pH*Invasives, data = invasion.exp.data.8, family = "poisson")
summary(glm.poisson.num.prop.disporella.erect.8)
plot(glm.poisson.num.prop.disporella.erect.8)
anova(glm.poisson.num.prop.disporella.erect.8, test = "Chi")


##### Week 10
glm.poisson.num.prop.disporella.erect.10<-glm(formula = (num.prop.disporella.erect) ~ min.10.pH*Invasives, data = invasion.exp.data.10, family = "poisson")
summary(glm.poisson.num.prop.disporella.erect.10)
plot(glm.poisson.num.prop.disporella.erect.10)
anova(glm.poisson.num.prop.disporella.erect.10, test = "Chi")



###nbinom

glm.nbinom.num.prop.disporella.erect.16<-glm.nb(formula = (num.prop.disporella.erect) ~ min.10.pH*Invasives, data = invasion.exp.data.16)
summary(glm.nbinom.num.prop.disporella.erect.16)
plot(glm.nbinom.num.prop.disporella.erect.16)
anova(glm.nbinom.num.prop.disporella.erect.16, test = "Chi")

visreg(glm.nbinom.num.prop.disporella.erect.16, xvar = "min.10.pH", by= "Invasives", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))

?glm.nb

####################### GAMMA For species div
library(glmmADMB)
require(car)
require(MASS)
library(poissonreg)
library(lmtest)
library(fitdistrplus)

#####################

# Week 12
head(invasion.exp.data.16)
?dgamma


gamma.16<-fitdistr(invasion.exp.data.16$num.prop.disporella.erect + 0.01, "gamma")
qqp(invasion.exp.data.16$num.prop.disporella.erect, "gamma", shape = gamma.16$estimate[[1]], rate = gamma.16$estimate[[2]])


glm.gamma.num.prop.disporella.erect.16<- glm(formula = num.prop.disporella.erect ~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "Gamma")
plot(glm.gamma.num.prop.disporella.erect.16)

summary(glm.gamma.num.prop.disporella.erect.16)
anova(glm.gamma.num.prop.disporella.erect.16, test = "F")

visreg(glm.gamma.num.prop.disporella.erect.16, xvar = "min.10.pH", by= "Invasives", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))






####################### NORMAL
library(lme4)


qqp(log(invasion.exp.data.16$Mussel.wet.weight+0.1), "norm")

qqp(invasion.exp.data.16$num.prop.disporella.erect, "lnorm")

qqp(sqrt(invasion.exp.data.16$Mussel.wet.weight), "norm")

fix(invasion.exp.data.16)


lm.Mussel.wet.weight<-lm(formula = log(Mussel.wet.weight+0.1) ~ min.10.pH*Invasives, data = invasion.exp.data.16)
summary(lm.Mussel.wet.weight)
plot(lm.Mussel.wet.weight)
anova(lm.Mussel.wet.weight, test = "F")


lognormal.16<-fitdistr(invasion.exp.data.16$num.prop.disporella.erect+0.01, "lognormal")
qqp(invasion.exp.data.16$num.prop.disporella.erect, "lnorm", lognormal.16$estimate)


glm.gamma.num.prop.disporella.erect.16<- glm(formula = num.prop.disporella.erect ~ min.10.pH*Invasives, data = invasion.exp.data.16, family = "lognormal")


qqp(invasion.exp.data.16$num.prop.disporella.erect, "lnorm")
?glm
?qqp

visreg(lm.num.prop.disporella.erect, xvar = "min.10.pH", by= "Invasives", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))




###might need gamlss pacage to do exponential functions
