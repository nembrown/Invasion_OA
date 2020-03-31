repeatedm.JC<-read.csv(file.choose())
#File: All week R input .csv

head(repeatedm.JC)

levels(repeatedm.JC$Treatment)[levels(repeatedm.JC$Treatment)=="AIR"] <- "Ambient CO2"
levels(repeatedm.JC$Treatment)[levels(repeatedm.JC$Treatment)=="CO2"] <- "Elevated CO2"

repeatedm.JC$hydrogen.concentration<-10^(-repeatedm.JC$min.10.pH)
repeatedm.JC$total<-100

repeatedm.JC.14<-repeatedm.JC[repeatedm.JC$Week==14, ]

repeatedm.JC.12<-repeatedm.JC[repeatedm.JC$Week==12, ]

repeatedm.JC.12to14<-rbind(repeatedm.JC[repeatedm.JC$Week==14, ], repeatedm.JC[repeatedm.JC$Week==12,])

repeatedm.JC.12to14.new <-aggregate(repeatedm.JC.12to14, by=list(repeatedm.JC.12to14$Tile.ID, repeatedm.JC.12to14$Mesocosm), 
                                    FUN=mean, na.rm=TRUE)

length(repeatedm.JC.12$Invasives)
head(repeatedm.JC.12to14.new)


repeatedm.JC.12to14.new$Invasives<-repeatedm.JC.14$Invasives
repeatedm.JC.12to14.new$CO2.Treatment<-repeatedm.JC.14$CO2.Treatment
repeatedm.JC.12to14.new$Treatment<-repeatedm.JC.14$Treatment
repeatedm.JC.12to14.new$Tile<-repeatedm.JC.14$Tile
repeatedm.JC.12to14.new$Tile.ID<-repeatedm.JC.14$Tile.ID
repeatedm.JC.12to14.new$Mesocosm<-repeatedm.JC.14$Mesocosm



library(plyr)

library(ggplot2)

library(doBy)



cdata.num_bot <- summaryBy(num_bot ~ Treatment + Week, data=repeatedm.JC, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.num_bot)[names(cdata.num_bot)=="num_bot.length"]<-"N"
cdata.num_bot$num_bot.se<-cdata.num_bot$num_bot.sd/sqrt(cdata.num_bot$N)
head(cdata.num_bot)


plot.num_bot<-ggplot(cdata.num_bot, aes(x=Week, y=num_bot.mean, linetype=Treatment, shape=Treatment)) + geom_errorbar(aes(ymin=num_bot.mean-num_bot.se, ymax=num_bot.mean+num_bot.se), width=0.25) + geom_line(size=1) + geom_point(aes(shape=Treatment),size=5, col="#666666") + scale_colour_manual(values=cbbPalette.all.purple)+scale_shape_manual(values=c(1, 2))
plot.num_bot<- plot.num_bot + labs(x = "Time (weeks)", y="Botryllus new recruitment\n(# colonies)")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+xlim(0,16)
plot.num_bot<-plot.num_bot + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.50,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.num_bot <- plot.num_bot+ theme(legend.position="none") + theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num_bot 


cbbPalette3 <- c("#000000", "#FF6600")


head(repeatedm.JC.14)

plot.num_bot.12<- ggplot(repeatedm.JC.12, aes(x=repeatedm.JC.14$hydrogen.concentration, y=(repeatedm.JC.12$num_bot)))  + geom_point(aes(shape=Treatment),size=5, col="#666666") + guides(fill=FALSE) + geom_smooth(colour="#666666", alpha=0.15, size=1.5,  method = 'glm', method.args = list(family = "poisson"))+scale_shape_manual(values=c(1,2))
plot.num_bot.12<- plot.num_bot.12 + theme_bw() + labs(x = "min.10.pH", y="Botryllus new recruitment\n(# colonies) 2nd peak")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num_bot.12<- plot.num_bot.12 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num_bot.12<- plot.num_bot.12+ theme(legend.position="none")+theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num_bot.12 <-plot.num_bot.12 + scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7, 10^-7.8), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7, 7.8))
plot.num_bot.12


plot.num_bot.14<- ggplot(repeatedm.JC.14, aes(x=repeatedm.JC.14$hydrogen.concentration, y=(repeatedm.JC.14$num_bot)))  + geom_point(aes(shape=Treatment),size=5, col="#666666") + guides(fill=FALSE) + geom_smooth(colour="#666666", alpha=0.15, size=1.5,  method = 'glm', method.args = list(family = "poisson"))+scale_shape_manual(values=c(1,2))
plot.num_bot.14<- plot.num_bot.14 + theme_bw() + labs(x = "min.10.pH", y="Botryllus new recruitment\n(# colonies) 2nd peak")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num_bot.14<- plot.num_bot.14 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num_bot.14<- plot.num_bot.14+ theme(legend.position="none")+theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num_bot.14 <-plot.num_bot.14 + scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7, 10^-7.8), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7, 7.8))
plot.num_bot.14

plot.num_bot.12to14<- ggplot(repeatedm.JC.12to14.new, aes(x=repeatedm.JC.14$hydrogen.concentration, y=(repeatedm.JC.12to14.new$num_bot)))  + geom_point(aes(shape=Treatment),size=5, col="#666666") + guides(fill=FALSE) + geom_smooth(colour="#666666", alpha=0.15, size=1.5,  method = 'glm', method.args = list(family = "poisson"))+scale_shape_manual(values=c(1,2))
plot.num_bot.12to14<- plot.num_bot.12to14 + theme_bw() + labs(x = "min.10.pH", y="Botryllus new recruitment\n(# colonies) 2nd peak")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.num_bot.12to14<- plot.num_bot.12to14 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.num_bot.12to14<- plot.num_bot.12to14+ theme(legend.position="none")+theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.num_bot.12to14 <-plot.num_bot.12to14 + scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7, 10^-7.8), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7, 7.8))
plot.num_bot.12to14

warnings()

plot.botryllid.16.invasion<- ggplot(invasion.exp.data.16, aes(x=invasion.exp.data.16$hydrogen.concentration, y=(invasion.exp.data.16$botryllid/invasion.exp.data.16$total), col="#666666")) + geom_point(size=5,aes(shape=Treatment), col="#666666") + guides(fill=FALSE) + geom_smooth(aes(lty=Invasives, fill=Invasives, weight=total), colour="#666666", alpha=0.15,size=1.5, method="glm", method.args = list(family = "binomial")) +scale_shape_manual(values=c(1,19,2,17))+scale_linetype_manual(values=c("dashed", "solid"))+ scale_fill_manual(values=c( "#FFFFFF","#666666"))
plot.botryllid.16.invasion<- plot.botryllid.16.invasion + theme_bw() + xlab(bquote('Minimum 10th percentile pH')) + ylab("Proportion cover Botryllus")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.botryllid.16.invasion<- plot.botryllid.16.invasion + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.botryllid.16.invasion<- plot.botryllid.16.invasion+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.botryllid.16.invasion <-plot.botryllid.16.invasion +  scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7))
plot.botryllid.16.invasion 



library(cowplot)
plot_grid(plot.num_bot, plot.num_bot.12to14 , ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))


plot_grid(plot.botryllid, plot.botryllid.16.invasion ,plot.num_bot, plot.num_bot.12to14 , ncol=2, align='h', labels=c('(a)', '(b)', '(c)','(d)', label_size=12))

##############Stats
library(glmmADMB)
require(car)
require(MASS)
library(poissonreg)
library(lmtest)
library(fitdistrplus)

repeatedm.JC.12<-repeatedm.JC[repeatedm.JC$Week==12, ]
head(repeatedm.JC.12)


poisson.12<-fitdistr(repeatedm.JC.12$num_bot, "Poisson")
qqp(repeatedm.JC.12$num_bot, "pois", poisson.12$estimate)

poisson<-fitdistr(repeatedm.JC$num_bot, "Poisson")
qqp(repeatedm.JC$num_bot, "pois", poisson$estimate)

nbinom <- fitdistr(repeatedm.JC$num_bot, "Negative Binomial")
qqp(repeatedm.JC$num_bot, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])


glm.poisson.num_bot.12<-glm(formula = (num_bot) ~ Treatment, data = repeatedm.JC.12, family = "poisson")
summary(glm.poisson.num_bot.12)
plot(glm.poisson.num_bot.12)
anova(glm.poisson.num_bot.12, test = "Chi")

repeatedm.JC$Week<-as.factor(repeatedm.JC$Week)
repeatedm.JC$Treatment<-as.factor(repeatedm.JC$Treatment)
repeatedm.JC$Mesocosm<-as.factor(repeatedm.JC$Mesocosm)
###### Do I have Week as continous or not??? 
### Or do I have 1|Meso and not Week|Meso
### Have to use glmm.... 
#### Could I treat week as a continuous variable?? 

repeatedm.JC$Week<-as.numeric(repeatedm.JC$Week)
?glmmadmb

library(glmmADMB)
install.packages("reshape2")
library(reshape2)

is.factor(repeatedm.JC$Week)

### Negative binomial worked

glm.new<- glmmadmb(num_bot ~ Treatment*Week + (Week|Mesocosm), data=repeatedm.JC, family="nbinom", zeroInflation =F, admb.opts = admbControl(shess = FALSE, noinit = FALSE, impSamp=100,maxfn=1000,imaxfn=500,maxph=5))
summary(glm.new)
plot(resid(glm.new)~fitted(glm.new))





###Week 8

repeatedm.JC.8<-repeatedm.JC[repeatedm.JC$Week==8, ]
poisson.8<-fitdistr(repeatedm.JC.8$num_bot, "Poisson")
qqp(repeatedm.JC.8$num_bot, "pois", poisson.8$estimate)

glm.poisson.num_bot.8<-glm(formula = (num_bot) ~ Treatment, data = repeatedm.JC.8, family = "poisson")
summary(glm.poisson.num_bot.8)
plot(glm.poisson.num_bot.8)
anova(glm.poisson.num_bot.8, test = "Chi")
