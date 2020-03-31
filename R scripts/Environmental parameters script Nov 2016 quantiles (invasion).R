
setwd("/Users/Norah Brown/Documents/Projects/Summer 2015 food availability experiments")

repeat.env.para.invasion<-read.csv(file.choose())
#File: salinity input .csv

head(repeat.env.para.invasion)


repeat.env.para.invasion$min.10.pH<-quantile(repeat.env.para.invasion$pH,0.10)


library(grid)

library(plyr)

library(ggplot2)

library(doBy)



library(plyr)
do.call("rbind", tapply(repeat.env.para.invasion$pH, repeat.env.para.invasion$Tile.ID, quantile))

repeat.env.para.invasion.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion$pH, repeat.env.para.invasion$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.pH.quantiles)

write.csv(repeat.env.para.invasion.pH.quantiles, "repeat.env.para.invasion.16.pH.quantiles.csv")

##Week 10

repeat.env.para.invasion.14<- repeat.env.para.invasion[repeat.env.para.invasion$Week != 16, ]

repeat.env.para.invasion.14.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.14$pH, repeat.env.para.invasion.14$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.14.pH.quantiles)


write.csv(repeat.env.para.invasion.14.pH.quantiles, "repeat.env.para.invasion.14.pH.quantiles.csv")


##Week 12

repeat.env.para.invasion.12<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=12, ]

repeat.env.para.invasion.12.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.12$pH, repeat.env.para.invasion.12$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.12.pH.quantiles)


write.csv(repeat.env.para.invasion.12.pH.quantiles, "repeat.env.para.invasion.12.pH.quantiles.csv")



##Week 10

repeat.env.para.invasion.10<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=10, ]

repeat.env.para.invasion.10.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.10$pH, repeat.env.para.invasion.10$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.10.pH.quantiles)


write.csv(repeat.env.para.invasion.10.pH.quantiles, "repeat.env.para.invasion.10.pH.quantiles.csv")


##Week 8

repeat.env.para.invasion.8<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=8, ]

repeat.env.para.invasion.8.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.8$pH, repeat.env.para.invasion.8$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.8.pH.quantiles)


write.csv(repeat.env.para.invasion.8.pH.quantiles, "repeat.env.para.invasion.8.pH.quantiles.csv")


##Week 6

repeat.env.para.invasion.6<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=6, ]

repeat.env.para.invasion.6.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.6$pH, repeat.env.para.invasion.6$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.6.pH.quantiles)


write.csv(repeat.env.para.invasion.6.pH.quantiles, "repeat.env.para.invasion.6.pH.quantiles.csv")

##Week 4

repeat.env.para.invasion.4<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=4, ]

repeat.env.para.invasion.4.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.4$pH, repeat.env.para.invasion.4$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.4.pH.quantiles)


write.csv(repeat.env.para.invasion.4.pH.quantiles, "repeat.env.para.invasion.4.pH.quantiles.csv")


##Week 2

repeat.env.para.invasion.2<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=2, ]

repeat.env.para.invasion.2.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.2$pH, repeat.env.para.invasion.2$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.2.pH.quantiles)


write.csv(repeat.env.para.invasion.2.pH.quantiles, "repeat.env.para.invasion.2.pH.quantiles.csv")


##Week 0

repeat.env.para.invasion.0<- repeat.env.para.invasion[repeat.env.para.invasion$Week <=0, ]

repeat.env.para.invasion.0.pH.quantiles<-do.call("rbind", tapply(repeat.env.para.invasion.0$pH, repeat.env.para.invasion.0$Tile.ID, quantile, c(0.05, 0.1, 0.2, .5, .75, .90)))

head(repeat.env.para.invasion.0.pH.quantiles)


write.csv(repeat.env.para.invasion.0.pH.quantiles, "repeat.env.para.invasion.0.pH.quantiles.csv")









###psalinity = salinity
cdatasalinity <- summaryBy(salinity ~ Combined.Treatment + CO2 + Food.quality + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
head(cdatasalinity)
#summarize the mussel data by Combined.Treatment and week in terms of mean and standard dev

names(cdatasalinity)[names(cdatasalinity)=="salinity.length"]<-"N"
head(cdatasalinity)
#changed column heading to N

cdatasalinity$salinity.se<-cdatasalinity$salinity.sd/sqrt(cdatasalinity$N)
head(cdatasalinity)
#created standard error from SD

#cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")

psalinity<-ggplot(cdatasalinity, aes(x=Week, y=salinity.mean, color = Food.quality, shape=CO2, linetype=CO2)) + geom_errorbar(aes(ymin=salinity.mean-salinity.se, ymax=salinity.mean+salinity.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)
#rename my plot to p, makes it easier to work with

psalinity<- psalinity + labs(x = "Time (weeks)", y="Mean salinity (ppt)")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

psalinity<- psalinity + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.81,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

psalinity <- psalinity + theme(legend.position="none")

psalinity




###pTemperature = Temperature
cdatatemp <- summaryBy(temp ~ Combined.Treatment + CO2 + Food.quality+ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatatemp)[names(cdatatemp)=="temp.length"]<-"N"
cdatatemp$temp.se<-cdatatemp$temp.sd/sqrt(cdatatemp$N)
pTemperature<-ggplot(cdatatemp, aes(x=Week, y=temp.mean, color = Food.quality, shape=CO2, linetype=CO2)) + geom_errorbar(aes(ymin=temp.mean-temp.se, ymax=temp.mean+temp.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

pTemperature<- pTemperature + labs(x = "Time (weeks)", y="Mean temperature (°C)")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

pTemperature<- pTemperature + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pTemperature



####ppH = pH

head(repeat.env.para.invasion)

cdatapH <- summaryBy(pH ~ CO2 + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatapH)[names(cdatapH)=="pH.length"]<-"N"
cdatapH$pH.se<-cdatapH$pH.sd/sqrt(cdatapH$N)
head(cdatapH)

cbbPalette7 <- c( "#636363", "#bdbdbd")


ppH<-ggplot(cdatapH, aes(x=Week, y=pH.mean, colour = CO2)) + geom_errorbar(aes(ymin=pH.mean-pH.se, ymax=pH.mean+pH.se), colour="black", width=0.25) + geom_line(size=0.5) + geom_point(size=3) + scale_colour_manual(values=cbbPalette7)

ppH<- ppH + labs(x = "Time (Weeks)", y="Mean pH")+ theme_bw()+ theme(text = element_text(size=11), axis.text = element_text(size=11)) + theme(axis.title.y = element_text(angle=90))

ppH<- ppH + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=11), legend.position=c(.90,.90))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

ppH


#####
cdatapH.Combined <- summaryBy(pH ~ Combined.Treatment + CO2 + Food.quality + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatapH.Combined)[names(cdatapH.Combined)=="pH.length"]<-"N"
cdatapH.Combined$pH.se<-cdatapH.Combined$pH.sd/sqrt(cdatapH.Combined$N)
head(cdatapH.Combined)

#cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")


ppH.Combined<-ggplot(cdatapH.Combined, aes(x=Week, y=pH.mean, color = Food.quality, shape=CO2, linetype=CO2)) + geom_errorbar(aes(ymin=pH.mean-pH.se, ymax=pH.mean+pH.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

ppH.Combined<- ppH.Combined + labs(x = "Time (Weeks)", y="Mean pH")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

ppH.Combined<- ppH.Combined + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.90,.90))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

ppH.Combined <- ppH.Combined + theme(legend.position="none")
ppH.Combined


multiplot(ppH.Combined, pTemperature, psalinity, cols = 1)




multiplot(ppH.Combined, psalinity,pTemperature, cols=1)

require(cowplot)
theme_set(theme_classic())

plot_grid(ppH, pTemperature, psalinity, ncol=1, align='v', labels=c('(a)', '(b)', '(c)'), label_size=12)


plot_grid(p5, p6, ncol=2, align='h', labels=c('(a)', '(b)'), label_size=12)



##### Weight minus brick
head(repeat.env.para.invasion)

cdatabiomass.tile <- summaryBy(biomass.tile ~ Combined.Treatment +CO2 + Food.quality+ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatabiomass.tile)[names(cdatabiomass.tile)=="biomass.tile.length"]<-"N"
cdatabiomass.tile$biomass.tile.se<-cdatabiomass.tile$biomass.tile.sd/sqrt(cdatabiomass.tile$N)
head(cdatabiomass.tile)

#cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")


pbiomass.tile<-ggplot(cdatabiomass.tile, aes(x=Week, y=biomass.tile.mean, colour = Food.quality,  shape=CO2, linetype=CO2)) + geom_errorbar(aes(ymin=biomass.tile.mean-biomass.tile.se, ymax=biomass.tile.mean+biomass.tile.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

pbiomass.tile<- pbiomass.tile + labs(x = "Time (Weeks)", y="Mean biomass (g)")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

pbiomass.tile<- pbiomass.tile + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.25,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pbiomass.tile+ theme(legend.position="none")


##### WEight over time ... with brick
head(repeat.env.para.invasion)

cdatawet.weight..g. <- summaryBy(wet.weight..g. ~ Combined.Treatment + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatawet.weight..g.)[names(cdatawet.weight..g.)=="wet.weight..g..length"]<-"N"
cdatawet.weight..g.$wet.weight..g..se<-cdatawet.weight..g.$wet.weight..g..sd/sqrt(cdatawet.weight..g.$N)
head(cdatawet.weight..g.)

cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")


pwet.weight..g.<-ggplot(cdatawet.weight..g., aes(x=Week, y=wet.weight..g..mean, colour = Combined.Treatment)) + geom_errorbar(aes(ymin=wet.weight..g..mean-wet.weight..g..se, ymax=wet.weight..g..mean+wet.weight..g..se), colour="black", width=0.25) + geom_line(size=0.5) + geom_point(size=3) + scale_colour_manual(values=cbbPalette6)

pwet.weight..g.<- pwet.weight..g. + labs(x = "Time (Weeks)", y="Mean wet.weight..g.")+ theme_bw()+ theme(text = element_text(size=11), axis.text = element_text(size=11)) + theme(axis.title.y = element_text(angle=90))

pwet.weight..g.<- pwet.weight..g. + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=11), legend.position=c(.90,.90))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pwet.weight..g.



##### Weight minus brick
head(repeat.env.para.invasion)

cdatabiomass.tile <- summaryBy(biomass.tile ~ Combined.Treatment + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatabiomass.tile)[names(cdatabiomass.tile)=="biomass.tile.length"]<-"N"
cdatabiomass.tile$biomass.tile.se<-cdatabiomass.tile$biomass.tile.sd/sqrt(cdatabiomass.tile$N)
head(cdatabiomass.tile)

cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")


pbiomass.tile<-ggplot(cdatabiomass.tile, aes(x=Week, y=biomass.tile.mean, colour = Combined.Treatment)) + geom_errorbar(aes(ymin=biomass.tile.mean-biomass.tile.se, ymax=biomass.tile.mean+biomass.tile.se), colour="black", width=0.25) + geom_line(size=0.5) + geom_point(size=3) + scale_colour_manual(values=cbbPalette6)

pbiomass.tile<- pbiomass.tile + labs(x = "Time (Weeks)", y="Mean biomass.tile")+ theme_bw()+ theme(text = element_text(size=11), axis.text = element_text(size=11)) + theme(axis.title.y = element_text(angle=90))

pbiomass.tile<- pbiomass.tile + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=11), legend.position=c(.25,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pbiomass.tile


##### Weight with just CO2
head(repeat.env.para.invasion)

cdatabiomass.tile.CO2 <- summaryBy(biomass.tile ~ CO2 + Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatabiomass.tile.CO2)[names(cdatabiomass.tile.CO2)=="biomass.tile.length"]<-"N"
cdatabiomass.tile.CO2$biomass.tile.se<-cdatabiomass.tile.CO2$biomass.tile.sd/sqrt(cdatabiomass.tile.CO2$N)
head(cdatabiomass.tile.CO2)


pbiomass.tile.CO2<-ggplot(cdatabiomass.tile.CO2, aes(x=Week, y=biomass.tile.mean, colour = CO2)) + geom_errorbar(aes(ymin=biomass.tile.mean-biomass.tile.se, ymax=biomass.tile.mean+biomass.tile.se), colour="black", width=0.25) + geom_line(size=0.5) + geom_point(size=3) + scale_colour_manual(values=cbbPalette7)

pbiomass.tile.CO2<- pbiomass.tile.CO2 + labs(x = "Time (Weeks)", y="Mean biomass.tile.CO2")+ theme_bw()+ theme(text = element_text(size=11), axis.text = element_text(size=11)) + theme(axis.title.y = element_text(angle=90))

pbiomass.tile.CO2<- pbiomass.tile.CO2 + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=11), legend.position=c(.25,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pbiomass.tile.CO2













#Multiplot how to
# Make a list from the ... arguments and plotlist
plots <- c(list(...), plotlist)

numPlots = length(plots)

# If layout is NULL, then use 'cols' to determine layout
if (is.null(layout)) {
  # Make the panel
  # ncol: Number of columns of plots
  # nrow: Number of rows needed, calculated from # of cols
  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                   ncol = cols, nrow = ceiling(numPlots/cols))
}

if (numPlots==1) {
  print(plots[[1]])
  
} else {
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    # Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    
    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                    layout.pos.col = matchidx$col))
  }
}
}






#### MARINA ONLY 
head(repeat.env.para.invasion)
###pSalinity..psu = Salinity..psu
cdataSalinity..psu <- summaryBy(Salinity..psu. ~  Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
head(cdataSalinity..psu)
#summarize the mussel data by Combined.Treatment and week in terms of mean and standard dev

names(cdataSalinity..psu)[names(cdataSalinity..psu)=="Salinity..psu..length"]<-"N"
head(cdataSalinity..psu)
#changed column heading to N

cdataSalinity..psu$Salinity..psu.se<-cdataSalinity..psu$Salinity..psu..sd/sqrt(cdataSalinity..psu$N)
head(cdataSalinity..psu)
#created standard error from SD

#cbbPalette6 <- c( "#de2d26", "#fc9272",  "#31a354","#a1d996","#636363", "#bdbdbd")

pSalinity..psu<-ggplot(cdataSalinity..psu, aes(x=Week, y=Salinity..psu..mean)) + geom_errorbar(aes(ymin=Salinity..psu..mean-Salinity..psu.se, ymax=Salinity..psu..mean+Salinity..psu.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)
#rename my plot to p, makes it easier to work with

pSalinity..psu<- pSalinity..psu + labs(x = "Time (weeks)", y="Mean Salinity marina (ppt))")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

pSalinity..psu<- pSalinity..psu + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.81,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pSalinity..psu <- pSalinity..psu + theme(legend.position="none")

pSalinity..psu


###ptemperature.Cerature = temperature.C
cdatatemperature.C. <- summaryBy(temperature.C. ~ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatatemperature.C.)[names(cdatatemperature.C.)=="temperature.C..length"]<-"N"
cdatatemperature.C.$temperature.C..se<-cdatatemperature.C.$temperature.C..sd/sqrt(cdatatemperature.C.$N)
ptemperature.C.erature<-ggplot(cdatatemperature.C., aes(x=Week, y=temperature.C..mean)) + geom_errorbar(aes(ymin=temperature.C..mean-temperature.C..se, ymax=temperature.C..mean+temperature.C..se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

ptemperature.C.erature<- ptemperature.C.erature + labs(x = "Time (weeks)", y="Mean temperature marina (°C)")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

ptemperature.C.erature<- ptemperature.C.erature + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

ptemperature.C.erature


###ptemperature.Cerature = temperature.C
cdatapH <- summaryBy(pH ~ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatapH)[names(cdatapH)=="pH.length"]<-"N"
cdatapH$pH.se<-cdatapH$pH.sd/sqrt(cdatapH$N)
ppHerature<-ggplot(cdatapH, aes(x=Week, y=pH.mean)) + geom_errorbar(aes(ymin=pH.mean-pH.se, ymax=pH.mean+pH.se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

ppHerature<- ppHerature + labs(x = "Time (weeks)", y="Mean pH marina")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

ppHerature<- ppHerature + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

ppHerature



###ptemperature.Cerature = temperature.C
cdatapCO2..microatm. <- summaryBy(pCO2..microatm. ~ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdatapCO2..microatm.)[names(cdatapCO2..microatm.)=="pCO2..microatm..length"]<-"N"
cdatapCO2..microatm.$pCO2..microatm..se<-cdatapCO2..microatm.$pCO2..microatm..sd/sqrt(cdatapCO2..microatm.$N)
ppCO2..microatm.<-ggplot(cdatapCO2..microatm., aes(x=Week, y=pCO2..microatm..mean)) + geom_errorbar(aes(ymin=pCO2..microatm..mean-pCO2..microatm..se, ymax=pCO2..microatm..mean+pCO2..microatm..se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

ppCO2..microatm.<- ppCO2..microatm. + labs(x = "Time (weeks)", y="Mean pCO2..microatm. marina")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

ppCO2..microatm.<- ppCO2..microatm. + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

ppCO2..microatm.


###
cdataTCO2..micromol.kgSW. <- summaryBy(TCO2..micromol.kgSW. ~ Week, data=repeat.env.para.invasion, FUN=c(length,mean,sd))
names(cdataTCO2..micromol.kgSW.)[names(cdataTCO2..micromol.kgSW.)=="TCO2..micromol.kgSW..length"]<-"N"
cdataTCO2..micromol.kgSW.$TCO2..micromol.kgSW..se<-cdataTCO2..micromol.kgSW.$TCO2..micromol.kgSW..sd/sqrt(cdataTCO2..micromol.kgSW.$N)
pTCO2..micromol.kgSW.<-ggplot(cdataTCO2..micromol.kgSW., aes(x=Week, y=TCO2..micromol.kgSW..mean)) + geom_errorbar(aes(ymin=TCO2..micromol.kgSW..mean-TCO2..micromol.kgSW..se, ymax=TCO2..micromol.kgSW..mean+TCO2..micromol.kgSW..se), colour="black", width=0.25) + geom_line(size=1) + geom_point(size=5) #+ scale_colour_manual(values=cbbPalette6)

pTCO2..micromol.kgSW.<- pTCO2..micromol.kgSW. + labs(x = "Time (weeks)", y="Mean TCO2..micromol.kgSW. marina")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))

pTCO2..micromol.kgSW.<- pTCO2..micromol.kgSW. + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.5, "cm"))

pTCO2..micromol.kgSW.


