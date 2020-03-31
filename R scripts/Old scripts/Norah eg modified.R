##Sept. 27 2013##

library(vegan)
library(sciplot)


setwd("/Users/Norah Brown/Documents/Projects/Summer 2015 food availability experiments")

par(family="serif",las=1)

head(species)
species <- read.csv(file.choose())
row.names(species) <- species$Mesocosm
species <- species[,-1]

just.species<- species
head(just.species)


##species richness##
species$num.species <- specnumber(just.species)


##species accumulation##

sp2 <- specaccum(just.species,"random")
plot(sp2)
summary(sp2)
coef(sp2)
plot(sp2, ci.type="poly", xlab="Number of Sites",col="blue", lwd=2, ci.lty=0, ci.col="lightblue",ylab="Number of species")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

##shnnon diversity##

species$shannon.diversity <- diversity(just.species,"shannon")


write.csv(species, "Oct.2015.meso.inventory.species.csv")

##now mds and stuff##
##read in environment dataset##

compiled.data <- read.csv(file.choose())
head(compiled.data)
row.names(compiled.data) <- compiled.data$ID
compiled.data <- compiled.data[,-1]


transformed <- decostand(just.species,method='hellinger')
ord.mine <- metaMDS(transformed,distance="euclidean",trymax=100)

plot(ord.mine, type = "n",cex=.25,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))

##this is the specific sites/tiles represented in the ordination##
site.scaling <- as.data.frame(ord.mine$points)
points(site.scaling,pch=16)


ef2 <- envfit(ord.mine, compiled.data)
ef2
plot(ef2, p.max = 0.05,cex=.55)

site.scaling$ID <- row.names(site.scaling)
compiled.data$ID <- row.names(compiled.data)

new.compiled <- merge(site.scaling,compiled.data,by=c("ID"),all.x=F,all.y=F)

##now replot ordination, with sites colour-coded##

plot(ord.mine, type = "n",cex=.25,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
points(new.compiled$MDS1,new.compiled$MDS2,col=new.compiled$Factor,pch=16)

##now look at stats##

#### READ in new data files that have all the weeks and WITHOUT presence data and with

species.all.with <- read.csv(file.choose())
species.all.without<-read.csv(file.choose())
compiled.all<-read.csv(file.choose())


row.names(compiled.all) <- compiled.all$ID
compiled.all <- compiled.all[,-1]
head(compiled.all)

#### ALL.WITH
row.names(species.all.with) <- species.all.with$ID
species.all.with <- species.all.with[,-1]

just.species.all.with <- species.all.with

##species richness##
species.all.with$num.species <- specnumber(just.species.all.with)


##shnnon diversity##

species.all.with$shannon.diversity <- diversity(just.species.all.with,"shannon")

##### ALL.WITHOUT
row.names(species.all.without) <- species.all.without$ID
species.all.without <- species.all.without[,-1]

just.species.all.without <- species.all.without

##species richness##
species.all.without$num.species <- specnumber(just.species.all.without)


##shnnon diversity##

species.all.without$shannon.diversity <- diversity(just.species.all.without,"shannon")

#######################

shannon.div <- as.data.frame(species.all.with$shannon.diversity)
shannon.div$ID <- row.names(species.all.with)
head(shannon.div)


compiled.all$ID <- row.names(compiled.all)
new.compiled.all <- merge(shannon.div,compiled.all,by=c("ID"),all.x=F,all.y=F)

##### GRAPHING
head(new.compiled.all)
windows()
lineplot.CI(Week,species.all.with$shannon.diversity,group=Factor,data=new.compiled.all,xlab="Week",ylab="Shannon.diversity",trace.label="Factor")

my.aov <- glm(species.all.with$shannon.diversity~as.factor(Treatment)*Invasives,data=new.compiled.all,family="gaussian")
summary(my.aov)
anova(my.aov,test="F")

head(new.compiled.all)
head(just.species)


######
Try a glm:

my.glm2 <- aov(species.all.with$shannon.diversity ~ Treatment*as.ordered(Week)*Invasives + Error(Treatment/Mesocosm/Week/Invasives),family="Gamma",data=new.compiled.all)

summary(my.glm2)

####
species$ID <- row.names(site.scaling)
compiled.data$ID <- row.names(compiled.data)
species.compiled<-merge(species,compiled.data,by=c("ID"),all.x=F,all.y=F)

head(species.compiled)
my.glm <- aov(shannon.diversity ~ Treatment*as.ordered(Week)*Invasives + Error(Treatment/Mesocosm/Week/Invasives),family="Gamma",data=new.compiled.all)
summary(my.glm)

hist(species.all.with$shannon.diversity)


##### TRY with different species
head(species.all.without)
species.all.without$ID <- row.names(species.all.without)
head(species.all.without)

compiled.all.without<-merge(species.all.without,compiled.all, by=c("ID"),all.x=F,all.y=F)
head(compiled.all.without)

lineplot.CI(Week,botryllid,group=Factor,data=compiled.all.without,xlab="Week",ylab="Botryllid cover",trace.label="Factor")

hist(compiled.all.without$botryllid)
### NOT NORMAL AT ALL ... therefor try a tranformation first

compiled.all.without$botryllid.arcsin<-asin(sqrt(compiled.all.without$botryllid))
hist(compiled.all.without$botryllid.arcsin)
##### DID NOT HELP AT ALL

glm.botryllid <- aov(botryllid ~ Treatment*Invasives*as.ordered(Week) + Error(Treatment/Mesocosm/Week/Invasives),family="Gamma",data=compiled.all.without)
summary(glm.botryllid)

TRY AN LM WITH TRANSFORMED DATA
aov(transformed.percent.cover ~ CO2*as.ordered(WEEK)*TUNI + Error(CO2/MESO/WEEK/TUNI))

my.lm <- lm(botryllid.arcsin~Treatment*as.ordered(Week)*Invasives ,data=compiled.all.without)
plot(my.lm)

#### PROBLEM with KAties code is that you can't do Error in a linear model, you need to do an lme --> linear mixed effects model. 


### TRY LMER PACKAGE

library(lme4)


#### TRY not nesting week & treatment to make it simpler
glm.botryllid <- aov(botryllid ~ Treatment*Invasives*as.ordered(Week) + Error(Mesocosm/Invasives),data=compiled.all.without)
summary(glm.botryllid)

glm.botryllid2 <- aov(botryllid ~ as.ordered(Week)*Treatment*Invasives + Error(Mesocosm/Invasives),data=compiled.all.without)
summary(glm.botryllid2)
plot(glm.botryllid2)
residuals(glm.botryllid2)




######################
Plotting other things

lineplot.CI(Week,protozoa,group=Factor,data=compiled.all.without,xlab="Week",ylab="Botryllid cover",trace.label="Factor")


fix(species.all.with)

write.csv(species.all.with, "test.csv")
