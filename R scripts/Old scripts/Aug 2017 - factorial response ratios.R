library(metafor)

invasion.exp.data<- read.csv(file.choose())

invasion.exp.data.16<-invasion.exp.data[invasion.exp.data$Week==16, ]
head(invasion.exp.data.16)


library(tidyr)

library(plyr)
library(dplyr)
library(ggplot2)
library(doBy)
library(grid)
library(cowplot)



invasion.exp.data.16<-rename(invasion.exp.data.16, c("alive.bot"="botryllus", "schizo"="schizoporella", "alive.mem"="membranipora", "alive.barn"="balanus"))
head(invasion.exp.data.16)
length(invasion.exp.data.16)
data_long_invasion <- gather(invasion.exp.data.16, species, cover, c(12:47), factor_key=TRUE)
head(data_long_invasion)



length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.cover <- summaryBy(cover ~ species + Treatment + Week, data=data_long_invasion, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length2(x), na.rm=TRUE) } )
names(cdata.cover)[names(cdata.cover)=="cover.length"]<-"N"
cdata.cover$cover.se<-cdata.cover$cover.sd/sqrt(cdata.cover$N)
head(cdata.cover)
fix(cdata.cover)

data_wide_invasion <- reshape(cdata.cover, 
             timevar = "Treatment",
             idvar = c("species", "Week"),
             direction = "wide")
head(data_wide_invasion)



head(data_wide_invasion)
data_wide_invasion$cover.mean.AIRAbsent<-data_wide_invasion$cover.mean.AIRAbsent+0.001
data_wide_invasion$cover.mean.AIRPresent<-data_wide_invasion$cover.mean.AIRPresent+0.001
data_wide_invasion$cover.mean.CO2Absent<-data_wide_invasion$cover.mean.CO2Absent+0.001
data_wide_invasion$cover.mean.CO2Present<-data_wide_invasion$cover.mean.CO2Present+0.001


#### Calculate effect sizes

C<-data_wide_invasion$cover.mean.AIRAbsent
A<-data_wide_invasion$cover.mean.AIRPresent
B<-data_wide_invasion$cover.mean.CO2Absent
AB<-data_wide_invasion$cover.mean.CO2Present
sC<-data_wide_invasion$cover.sd.AIRAbsent
sA<-data_wide_invasion$cover.sd.AIRPresent
sB<-data_wide_invasion$cover.sd.CO2Absent
sAB<-data_wide_invasion$cover.sd.CO2Present
nC<-data_wide_invasion$N.AIRAbsent
nA<-data_wide_invasion$N.AIRPresent
nB<-data_wide_invasion$N.CO2Absent
nAB<-data_wide_invasion$N.CO2Present


data_wide_invasion$LnRR.invasives.overall<-log(A+AB)-log(C+B)
  
data_wide_invasion$LnRR.CO2.overall<-log(B+AB)-log(C+A)
  
data_wide_invasion$LnRR.interaction<-log(AB)-log(A)-log(B)+log(C)

data_wide_invasion$s.invasives.overall<-(((1/(A+AB))^2)*(((sA^2)/(nA))+((sAB^2)/nAB)))+(((1/(C+B))^2)*(((sC^2)/(nC))+((sB^2)/nB)))

data_wide_invasion$s.CO2.overall<-(((1/(B+AB))^2)*(((sB^2)/(nB))+((sAB^2)/nAB)))+(((1/(C+A))^2)*(((sC^2)/(nC))+((sA^2)/nA)))


data_wide_invasion$s.interaction<-((sA^2)/((A^2)*(nA))) + ((sB^2)/((B^2)*(nB))) + ((sAB^2)/((AB^2)*(nAB))) + ((sC^2)/((C^2)*(nC)))



head(data_wide_invasion)



library(metafor)

pd <- position_dodge(width = 0.4)
cbbPalette.purple.blue.red<- c( "#F8766D", "#0000FF", "#CC99CC")

theme_set(theme_classic())


length(data_wide_invasion)
length(data_wide_invasion_s)
length(data_wide_invasion_LnRR)
head(data_wide_invasion_s)

data_wide_invasion_LnRR <- gather(data_wide_invasion, LnRR.type, LnRR, c(23:25), factor_key=TRUE)
data_wide_invasion_s <- gather(data_wide_invasion, s.LnRR.type, s, c(26:28), factor_key=TRUE)
data_wide_invasion_bind <-as.data.frame(cbind(data_wide_invasion_LnRR, data_wide_invasion_s[,26:27]))

head(data_wide_invasion_bind)
data_wide_invasion_bind_CO2 <- data_wide_invasion_bind[data_wide_invasion_bind$LnRR.type=="LnRR.interaction",]
xvals_ <- data_wide_invasion_bind_CO2[with(data_wide_invasion_bind_CO2, order(LnRR)), ]$species
data_wide_invasion_bind$species2<-factor(data_wide_invasion_bind$species, levels=xvals_)


plot.invasion<-ggplot(data_wide_invasion_bind, aes(x=species,y=LnRR,col=LnRR.type))+ scale_colour_manual(values=cbbPalette.purple.blue.red)
plot.invasion<-plot.invasion + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LnRR - 1.96*sqrt(s)), ymax = (LnRR + 1.96*sqrt(s))), width = 0.1, position = pd)
plot.invasion<-plot.invasion+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion<-plot.invasion+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion<-plot.invasion+ xlab('Species') +ylab ('Effect size')
plot.invasion<-plot.invasion+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion<-plot.invasion +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion


########## combined weeks
invasion.exp.data.lastmonth<-invasion.exp.data[with(invasion.exp.data, Week>=8), ]
head(invasion.exp.data.lastmonth)

invasion.exp.data.8to16 <-aggregate(invasion.exp.data.lastmonth, by=list(invasion.exp.data.lastmonth$Tile.ID, invasion.exp.data.lastmonth$Mesocosm), 
                                    FUN=mean, na.rm=TRUE)

invasion.exp.data.8to16$Invasives<-invasion.exp.data.12$Invasives
invasion.exp.data.8to16$CO2.Treatment<-invasion.exp.data.12$CO2.Treatment
invasion.exp.data.8to16$Treatment<-invasion.exp.data.12$Treatment
invasion.exp.data.8to16$Tile<-invasion.exp.data.12$Tile
invasion.exp.data.8to16$Tile.ID<-invasion.exp.data.12$Tile.ID
invasion.exp.data.8to16$Mesocosm<-invasion.exp.data.12$Mesocosm

head(invasion.exp.data.8to16)
invasion.exp.data.8to16<-invasion.exp.data.8to16[,-1]

data_long_invasion_comboweeks <- gather(invasion.exp.data.8to16, species, cover, c(12:47), factor_key=TRUE)
head(data_long_invasion_comboweeks)




length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

cdata.cover.8to16 <- summaryBy(cover ~ species + Treatment + Week, data=data_long_invasion_comboweeks, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length2(x), na.rm=TRUE) } )
names(cdata.cover.8to16 )[names(cdata.cover.8to16 )=="cover.length"]<-"N"
cdata.cover.8to16 $cover.se<-cdata.cover.8to16$cover.sd/sqrt(cdata.cover$N)
head(cdata.cover.8to16)

data_wide_invasion_comboweeks <- reshape(cdata.cover.8to16, 
                              timevar = "Treatment",
                              idvar = c("species", "Week"),
                              direction = "wide")
head(data_wide_invasion_comboweeks)



head(data_wide_invasion_comboweeks)
data_wide_invasion_comboweeks$cover.mean.AIRAbsent<-data_wide_invasion_comboweeks$cover.mean.AIRAbsent+0.001
data_wide_invasion_comboweeks$cover.mean.AIRPresent<-data_wide_invasion_comboweeks$cover.mean.AIRPresent+0.001
data_wide_invasion_comboweeks$cover.mean.CO2Absent<-data_wide_invasion_comboweeks$cover.mean.CO2Absent+0.001
data_wide_invasion_comboweeks$cover.mean.CO2Present<-data_wide_invasion_comboweeks$cover.mean.CO2Present+0.001


#### Calculate effect sizes

C<-data_wide_invasion_comboweeks$cover.mean.AIRAbsent
A<-data_wide_invasion_comboweeks$cover.mean.AIRPresent
B<-data_wide_invasion_comboweeks$cover.mean.CO2Absent
AB<-data_wide_invasion_comboweeks$cover.mean.CO2Present
sC<-data_wide_invasion_comboweeks$cover.sd.AIRAbsent
sA<-data_wide_invasion_comboweeks$cover.sd.AIRPresent
sB<-data_wide_invasion_comboweeks$cover.sd.CO2Absent
sAB<-data_wide_invasion_comboweeks$cover.sd.CO2Present
nC<-data_wide_invasion_comboweeks$N.AIRAbsent
nA<-data_wide_invasion_comboweeks$N.AIRPresent
nB<-data_wide_invasion_comboweeks$N.CO2Absent
nAB<-data_wide_invasion_comboweeks$N.CO2Present


data_wide_invasion_comboweeks$LnRR.invasives.overall<-log(A+AB)-log(C+B)

data_wide_invasion_comboweeks$LnRR.CO2.overall<-log(B+AB)-log(C+A)

data_wide_invasion_comboweeks$LnRR.interaction<-log(AB)-log(A)-log(B)+log(C)

data_wide_invasion_comboweeks$s.invasives.overall<-(((1/(A+AB))^2)*(((sA^2)/(nA))+((sAB^2)/nAB)))+(((1/(C+B))^2)*(((sC^2)/(nC))+((sB^2)/nB)))

data_wide_invasion_comboweeks$s.CO2.overall<-(((1/(B+AB))^2)*(((sB^2)/(nB))+((sAB^2)/nAB)))+(((1/(C+A))^2)*(((sC^2)/(nC))+((sA^2)/nA)))


data_wide_invasion_comboweeks$s.interaction<-((sA^2)/((A^2)*(nA))) + ((sB^2)/((B^2)*(nB))) + ((sAB^2)/((AB^2)*(nAB))) + ((sC^2)/((C^2)*(nC)))



head(data_wide_invasion_comboweeks)



library(metafor)

pd <- position_dodge(width = 0.4)
cbbPalette.purple.blue.red<- c( "#F8766D", "#0000FF", "#CC99CC")

theme_set(theme_classic())


length(data_wide_invasion_comboweeks)
length(data_wide_invasion_comboweeks_s)
length(data_wide_invasion_comboweeks_LnRR)
head(data_wide_invasion_comboweeks)

data_wide_invasion_comboweeks_LnRR <- gather(data_wide_invasion_comboweeks, LnRR.type, LnRR, c(23:25), factor_key=TRUE)
data_wide_invasion_comboweeks_s <- gather(data_wide_invasion_comboweeks, s.LnRR.type, s, c(26:28), factor_key=TRUE)
data_wide_invasion_comboweeks_bind <-as.data.frame(cbind(data_wide_invasion_comboweeks_LnRR, data_wide_invasion_comboweeks_s[,26:27]))

head(data_wide_invasion_comboweeks_bind)
data_wide_invasion_comboweeks_bind_CO2 <- data_wide_invasion_comboweeks_bind[data_wide_invasion_comboweeks_bind$LnRR.type=="LnRR.interaction",]
xvals_ <- data_wide_invasion_comboweeks_bind_CO2[with(data_wide_invasion_comboweeks_bind_CO2, order(LnRR)), ]$species
data_wide_invasion_comboweeks_bind$species2<-factor(data_wide_invasion_comboweeks_bind$species, levels=xvals_)


plot.invasion_comboweeks<-ggplot(data_wide_invasion_comboweeks_bind, aes(x=species,y=LnRR,col=LnRR.type))+ scale_colour_manual(values=cbbPalette.purple.blue.red)
plot.invasion_comboweeks<-plot.invasion_comboweeks + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LnRR - 1.96*sqrt(s)), ymax = (LnRR + 1.96*sqrt(s))), width = 0.1, position = pd)
plot.invasion_comboweeks<-plot.invasion_comboweeks+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_comboweeks<-plot.invasion_comboweeks+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_comboweeks<-plot.invasion_comboweeks+ xlab('Species') +ylab ('Effect size')
plot.invasion_comboweeks<-plot.invasion_comboweeks+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_comboweeks<-plot.invasion_comboweeks +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_comboweeks




###################
invasion.exp.data.lastmonth<-invasion.exp.data[with(invasion.exp.data, Week>=8), ]
head(invasion.exp.data.lastmonth)

data_long_invasion_comboweeks.sep <- gather(invasion.exp.data.lastmonth, species, cover, c(12:47), factor_key=TRUE)
head(data_long_invasion_comboweeks.sep)

### Weeks separated
cdata.cover.8to16.sep <- summaryBy(cover ~ species + Treatment + Week, data=data_long_invasion_comboweeks.sep, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length2(x), na.rm=TRUE) } )
names(cdata.cover.8to16.sep )[names(cdata.cover.8to16.sep )=="cover.length"]<-"N"
cdata.cover.8to16.sep $cover.se<-cdata.cover.8to16.sep$cover.sd/sqrt(cdata.cover$N)
head(cdata.cover.8to16.sep)


data_wide_invasion_comboweeks.sep <- reshape(cdata.cover.8to16.sep, 
                                         timevar = "Treatment",
                                         idvar = c("species", "Week"),
                                         direction = "wide")
head(data_wide_invasion_comboweeks.sep)



head(data_wide_invasion_comboweeks.sep)
data_wide_invasion_comboweeks.sep$cover.mean.AIRAbsent<-data_wide_invasion_comboweeks.sep$cover.mean.AIRAbsent+0.001
data_wide_invasion_comboweeks.sep$cover.mean.AIRPresent<-data_wide_invasion_comboweeks.sep$cover.mean.AIRPresent+0.001
data_wide_invasion_comboweeks.sep$cover.mean.CO2Absent<-data_wide_invasion_comboweeks.sep$cover.mean.CO2Absent+0.001
data_wide_invasion_comboweeks.sep$cover.mean.CO2Present<-data_wide_invasion_comboweeks.sep$cover.mean.CO2Present+0.001


#### Calculate effect sizes

C<-data_wide_invasion_comboweeks.sep$cover.mean.AIRAbsent
A<-data_wide_invasion_comboweeks.sep$cover.mean.AIRPresent
B<-data_wide_invasion_comboweeks.sep$cover.mean.CO2Absent
AB<-data_wide_invasion_comboweeks.sep$cover.mean.CO2Present
sC<-data_wide_invasion_comboweeks.sep$cover.sd.AIRAbsent
sA<-data_wide_invasion_comboweeks.sep$cover.sd.AIRPresent
sB<-data_wide_invasion_comboweeks.sep$cover.sd.CO2Absent
sAB<-data_wide_invasion_comboweeks.sep$cover.sd.CO2Present
nC<-data_wide_invasion_comboweeks.sep$N.AIRAbsent
nA<-data_wide_invasion_comboweeks.sep$N.AIRPresent
nB<-data_wide_invasion_comboweeks.sep$N.CO2Absent
nAB<-data_wide_invasion_comboweeks.sep$N.CO2Present


data_wide_invasion_comboweeks.sep$LnRR.invasives.overall<-log(A+AB)-log(C+B)

data_wide_invasion_comboweeks.sep$LnRR.CO2.overall<-log(B+AB)-log(C+A)

data_wide_invasion_comboweeks.sep$LnRR.interaction<-log(AB)-log(A)-log(B)+log(C)

data_wide_invasion_comboweeks.sep$s.invasives.overall<-(((1/(A+AB))^2)*(((sA^2)/(nA))+((sAB^2)/nAB)))+(((1/(C+B))^2)*(((sC^2)/(nC))+((sB^2)/nB)))

data_wide_invasion_comboweeks.sep$s.CO2.overall<-(((1/(B+AB))^2)*(((sB^2)/(nB))+((sAB^2)/nAB)))+(((1/(C+A))^2)*(((sC^2)/(nC))+((sA^2)/nA)))


data_wide_invasion_comboweeks.sep$s.interaction<-((sA^2)/((A^2)*(nA))) + ((sB^2)/((B^2)*(nB))) + ((sAB^2)/((AB^2)*(nAB))) + ((sC^2)/((C^2)*(nC)))



head(data_wide_invasion_comboweeks.sep)



library(metafor)

pd <- position_dodge(width = 0.4)
cbbPalette.purple.blue.red<- c( "#F8766D", "#0000FF", "#CC99CC")
cbbPalette.blue.purple.red<- c(  "#0000FF", "#CC99CC","#F8766D")

theme_set(theme_classic())


length(data_wide_invasion_comboweeks.sep)
length(data_wide_invasion_comboweeks.sep_s)
length(data_wide_invasion_comboweeks.sep_LnRR)
head(data_wide_invasion_comboweeks.sep)

data_wide_invasion_comboweeks.sep_LnRR <- gather(data_wide_invasion_comboweeks.sep, LnRR.type, LnRR, c(23:25), factor_key=TRUE)
data_wide_invasion_comboweeks.sep_s <- gather(data_wide_invasion_comboweeks.sep, s.LnRR.type, s, c(26:28), factor_key=TRUE)
data_wide_invasion_comboweeks.sep_bind <-as.data.frame(cbind(data_wide_invasion_comboweeks.sep_LnRR, data_wide_invasion_comboweeks.sep_s[,26:27]))



#model.invasions for species

head(data_wide_invasion_comboweeks.sep_bind)

library(metafor)

### botryllid
model.invasion.botryllid <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="botryllid"), struct="UN")
summary(model.invasion.botryllid)

cdata.LRR.invasion.botryllid<-as.data.frame(cbind(model.invasion.botryllid$b, model.invasion.botryllid$ci.lb, model.invasion.botryllid$ci.ub))
colnames(cdata.LRR.invasion.botryllid)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.botryllid$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.botryllid$species<-c("botryllid", "botryllid", "botryllid")
rownames(cdata.LRR.invasion.botryllid) <- 1:nrow(cdata.LRR.invasion.botryllid)
head(cdata.LRR.invasion.botryllid)

### bot.eaten
model.invasion.bot.eaten <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="bot.eaten"), struct=c("UN", "UN"))
summary(model.invasion.bot.eaten)

cdata.LRR.invasion.bot.eaten<-as.data.frame(cbind(model.invasion.bot.eaten$b, model.invasion.bot.eaten$ci.lb, model.invasion.bot.eaten$ci.ub))
colnames(cdata.LRR.invasion.bot.eaten)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.bot.eaten$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.bot.eaten$species<-c("bot.eaten", "bot.eaten", "bot.eaten")
rownames(cdata.LRR.invasion.bot.eaten) <- 1:nrow(cdata.LRR.invasion.bot.eaten)
head(cdata.LRR.invasion.bot.eaten)

### membranipora
model.invasion.membranipora <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="membranipora"), struct=c("UN", "UN"))
summary(model.invasion.membranipora)

cdata.LRR.invasion.membranipora<-as.data.frame(cbind(model.invasion.membranipora$b, model.invasion.membranipora$ci.lb, model.invasion.membranipora$ci.ub))
colnames(cdata.LRR.invasion.membranipora)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.membranipora$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.membranipora$species<-c("membranipora", "membranipora", "membranipora")
rownames(cdata.LRR.invasion.membranipora) <- 1:nrow(cdata.LRR.invasion.membranipora)
head(cdata.LRR.invasion.membranipora)

### mem.dead
model.invasion.mem.dead <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="mem.dead"), struct=c("UN", "UN"))
summary(model.invasion.mem.dead)

cdata.LRR.invasion.mem.dead<-as.data.frame(cbind(model.invasion.mem.dead$b, model.invasion.mem.dead$ci.lb, model.invasion.mem.dead$ci.ub))
colnames(cdata.LRR.invasion.mem.dead)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.mem.dead$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.mem.dead$species<-c("mem.dead", "mem.dead", "mem.dead")
rownames(cdata.LRR.invasion.mem.dead) <- 1:nrow(cdata.LRR.invasion.mem.dead)
head(cdata.LRR.invasion.mem.dead)

### mem.eaten
model.invasion.mem.eaten <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="mem.eaten"), struct=c("UN", "UN"))
summary(model.invasion.mem.eaten)

cdata.LRR.invasion.mem.eaten<-as.data.frame(cbind(model.invasion.mem.eaten$b, model.invasion.mem.eaten$ci.lb, model.invasion.mem.eaten$ci.ub))
colnames(cdata.LRR.invasion.mem.eaten)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.mem.eaten$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.mem.eaten$species<-c("mem.eaten", "mem.eaten", "mem.eaten")
rownames(cdata.LRR.invasion.mem.eaten) <- 1:nrow(cdata.LRR.invasion.mem.eaten)
head(cdata.LRR.invasion.mem.eaten)


### red.bryo
model.invasion.red.bryo <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="red.bryo"), struct=c("UN", "UN"))
summary(model.invasion.red.bryo)

cdata.LRR.invasion.red.bryo<-as.data.frame(cbind(model.invasion.red.bryo$b, model.invasion.red.bryo$ci.lb, model.invasion.red.bryo$ci.ub))
colnames(cdata.LRR.invasion.red.bryo)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.red.bryo$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.red.bryo$species<-c("red.bryo", "red.bryo", "red.bryo")
rownames(cdata.LRR.invasion.red.bryo) <- 1:nrow(cdata.LRR.invasion.red.bryo)
head(cdata.LRR.invasion.red.bryo)

###
is.factor(data_wide_invasion_comboweeks.sep_bind$species)
###
#for (i in levels(data_wide_invasion_comboweeks.sep_bind$species)){
#  print("model.invasion.[i] <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=data_wide_invasion_comboweeks.sep_bind,subset=(species=="[i]"), struct=c("UN", "UN"))
#  summary(model.invasion.[i])"
 # }

### white.bryo didn't run
### fan.bryo didn't run


### hydroid
model.invasion.hydroid <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="hydroid"), struct=c("UN", "UN"))
summary(model.invasion.hydroid)

cdata.LRR.invasion.hydroid<-as.data.frame(cbind(model.invasion.hydroid$b, model.invasion.hydroid$ci.lb, model.invasion.hydroid$ci.ub))
colnames(cdata.LRR.invasion.hydroid)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.hydroid$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.hydroid$species<-c("hydroid", "hydroid", "hydroid")
rownames(cdata.LRR.invasion.hydroid) <- 1:nrow(cdata.LRR.invasion.hydroid)
head(cdata.LRR.invasion.hydroid)


### clam did not run

### oyster didn't run

### mussel
model.invasion.mussel <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="mussel"), struct=c("UN", "UN"))
summary(model.invasion.mussel)

cdata.LRR.invasion.mussel<-as.data.frame(cbind(model.invasion.mussel$b, model.invasion.mussel$ci.lb, model.invasion.mussel$ci.ub))
colnames(cdata.LRR.invasion.mussel)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.mussel$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.mussel$species<-c("mussel", "mussel", "mussel")
rownames(cdata.LRR.invasion.mussel) <- 1:nrow(cdata.LRR.invasion.mussel)
head(cdata.LRR.invasion.mussel)

### corella
model.invasion.corella <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="corella"), struct="UN")
summary(model.invasion.corella)

cdata.LRR.invasion.corella<-as.data.frame(cbind(model.invasion.corella$b, model.invasion.corella$ci.lb, model.invasion.corella$ci.ub))
colnames(cdata.LRR.invasion.corella)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.corella$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.corella$species<-c("corella", "corella", "corella")
rownames(cdata.LRR.invasion.corella) <- 1:nrow(cdata.LRR.invasion.corella)
head(cdata.LRR.invasion.corella)


### dead.corella didn't run


### protozoa
model.invasion.protozoa <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="protozoa"), struct="UN")
summary(model.invasion.protozoa)

cdata.LRR.invasion.protozoa<-as.data.frame(cbind(model.invasion.protozoa$b, model.invasion.protozoa$ci.lb, model.invasion.protozoa$ci.ub))
colnames(cdata.LRR.invasion.protozoa)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.protozoa$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.protozoa$species<-c("protozoa", "protozoa", "protozoa")
rownames(cdata.LRR.invasion.protozoa) <- 1:nrow(cdata.LRR.invasion.protozoa)
head(cdata.LRR.invasion.protozoa)


### serpulid
model.invasion.serpulid <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="serpulid"), struct="UN")
summary(model.invasion.serpulid)

cdata.LRR.invasion.serpulid<-as.data.frame(cbind(model.invasion.serpulid$b, model.invasion.serpulid$ci.lb, model.invasion.serpulid$ci.ub))
colnames(cdata.LRR.invasion.serpulid)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.serpulid$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.serpulid$species<-c("serpulid", "serpulid", "serpulid")
rownames(cdata.LRR.invasion.serpulid) <- 1:nrow(cdata.LRR.invasion.serpulid)
head(cdata.LRR.invasion.serpulid)

### barn
model.invasion.barn <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="barn"), struct="UN")
summary(model.invasion.barn)

cdata.LRR.invasion.barn<-as.data.frame(cbind(model.invasion.barn$b, model.invasion.barn$ci.lb, model.invasion.barn$ci.ub))
colnames(cdata.LRR.invasion.barn)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.barn$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.barn$species<-c("barn", "barn", "barn")
rownames(cdata.LRR.invasion.barn) <- 1:nrow(cdata.LRR.invasion.barn)
head(cdata.LRR.invasion.barn)

### bare
model.invasion.bare <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="bare"), struct="UN")
summary(model.invasion.bare)

cdata.LRR.invasion.bare<-as.data.frame(cbind(model.invasion.bare$b, model.invasion.bare$ci.lb, model.invasion.bare$ci.ub))
colnames(cdata.LRR.invasion.bare)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.bare$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.bare$species<-c("bare", "bare", "bare")
rownames(cdata.LRR.invasion.bare) <- 1:nrow(cdata.LRR.invasion.bare)
head(cdata.LRR.invasion.bare)

### nudibranch didnt work

### nudi.eggs
model.invasion.nudi.eggs <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="nudi.eggs"), struct="UN")
summary(model.invasion.nudi.eggs)

cdata.LRR.invasion.nudi.eggs<-as.data.frame(cbind(model.invasion.nudi.eggs$b, model.invasion.nudi.eggs$ci.lb, model.invasion.nudi.eggs$ci.ub))
colnames(cdata.LRR.invasion.nudi.eggs)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.nudi.eggs$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.nudi.eggs$species<-c("nudi.eggs", "nudi.eggs", "nudi.eggs")
rownames(cdata.LRR.invasion.nudi.eggs) <- 1:nrow(cdata.LRR.invasion.nudi.eggs)
head(cdata.LRR.invasion.nudi.eggs)


### nudi.hatched didn't work

### bubble.snail didn't work

### bubble.eggs didn't work 

### chiton didn't work

### slime didn't work 

### num.species.no.bot
model.invasion.num.species.no.bot <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.species.no.bot"), struct="UN")
summary(model.invasion.num.species.no.bot)

cdata.LRR.invasion.num.species.no.bot<-as.data.frame(cbind(model.invasion.num.species.no.bot$b, model.invasion.num.species.no.bot$ci.lb, model.invasion.num.species.no.bot$ci.ub))
colnames(cdata.LRR.invasion.num.species.no.bot)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.species.no.bot$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.species.no.bot$species<-c("num.species.no.bot", "num.species.no.bot", "num.species.no.bot")
rownames(cdata.LRR.invasion.num.species.no.bot) <- 1:nrow(cdata.LRR.invasion.num.species.no.bot)
head(cdata.LRR.invasion.num.species.no.bot)

### shannon.diversity.no.bot
model.invasion.shannon.diversity.no.bot <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="shannon.diversity.no.bot"), struct="UN")
summary(model.invasion.shannon.diversity.no.bot)

cdata.LRR.invasion.shannon.diversity.no.bot<-as.data.frame(cbind(model.invasion.shannon.diversity.no.bot$b, model.invasion.shannon.diversity.no.bot$ci.lb, model.invasion.shannon.diversity.no.bot$ci.ub))
colnames(cdata.LRR.invasion.shannon.diversity.no.bot)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.shannon.diversity.no.bot$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.shannon.diversity.no.bot$species<-c("shannon.diversity.no.bot", "shannon.diversity.no.bot", "shannon.diversity.no.bot")
rownames(cdata.LRR.invasion.shannon.diversity.no.bot) <- 1:nrow(cdata.LRR.invasion.shannon.diversity.no.bot)
head(cdata.LRR.invasion.shannon.diversity.no.bot)

### num.serpulid
model.invasion.num.serpulid <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.serpulid"), struct="UN")
summary(model.invasion.num.serpulid)

cdata.LRR.invasion.num.serpulid<-as.data.frame(cbind(model.invasion.num.serpulid$b, model.invasion.num.serpulid$ci.lb, model.invasion.num.serpulid$ci.ub))
colnames(cdata.LRR.invasion.num.serpulid)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.serpulid$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.serpulid$species<-c("num.serpulid", "num.serpulid", "num.serpulid")
rownames(cdata.LRR.invasion.num.serpulid) <- 1:nrow(cdata.LRR.invasion.num.serpulid)
head(cdata.LRR.invasion.num.serpulid)

### num.barn
model.invasion.num.barn <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.barn"), struct="UN")
summary(model.invasion.num.barn)

cdata.LRR.invasion.num.barn<-as.data.frame(cbind(model.invasion.num.barn$b, model.invasion.num.barn$ci.lb, model.invasion.num.barn$ci.ub))
colnames(cdata.LRR.invasion.num.barn)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.barn$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.barn$species<-c("num.barn", "num.barn", "num.barn")
rownames(cdata.LRR.invasion.num.barn) <- 1:nrow(cdata.LRR.invasion.num.barn)
head(cdata.LRR.invasion.num.barn)

### num.red.bryo
model.invasion.num.red.bryo <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.red.bryo"), struct="UN")
summary(model.invasion.num.red.bryo)

cdata.LRR.invasion.num.red.bryo<-as.data.frame(cbind(model.invasion.num.red.bryo$b, model.invasion.num.red.bryo$ci.lb, model.invasion.num.red.bryo$ci.ub))
colnames(cdata.LRR.invasion.num.red.bryo)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.red.bryo$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.red.bryo$species<-c("num.red.bryo", "num.red.bryo", "num.red.bryo")
rownames(cdata.LRR.invasion.num.red.bryo) <- 1:nrow(cdata.LRR.invasion.num.red.bryo)
head(cdata.LRR.invasion.num.red.bryo)

### num.white.bryo
model.invasion.num.white.bryo <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.white.bryo"), struct="UN")
summary(model.invasion.num.white.bryo)

cdata.LRR.invasion.num.white.bryo<-as.data.frame(cbind(model.invasion.num.white.bryo$b, model.invasion.num.white.bryo$ci.lb, model.invasion.num.white.bryo$ci.ub))
colnames(cdata.LRR.invasion.num.white.bryo)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.white.bryo$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.white.bryo$species<-c("num.white.bryo", "num.white.bryo", "num.white.bryo")
rownames(cdata.LRR.invasion.num.white.bryo) <- 1:nrow(cdata.LRR.invasion.num.white.bryo)
head(cdata.LRR.invasion.num.white.bryo)


### num.corella
model.invasion.num.corella <- rma.mv(LnRR, s, mods = ~factor(LnRR.type) - 1,random=~factor(LnRR.type)|Week,data=na.omit(data_wide_invasion_comboweeks.sep_bind),subset=(species=="num.corella"), struct="UN")
summary(model.invasion.num.corella)

cdata.LRR.invasion.num.corella<-as.data.frame(cbind(model.invasion.num.corella$b, model.invasion.num.corella$ci.lb, model.invasion.num.corella$ci.ub))
colnames(cdata.LRR.invasion.num.corella)<- c("estimate", "lower", "upper")
cdata.LRR.invasion.num.corella$LnRR.type<-c("Invasives.overall", "CO2.overall", "Interaction")
cdata.LRR.invasion.num.corella$species<-c("num.corella", "num.corella", "num.corella")
rownames(cdata.LRR.invasion.num.corella) <- 1:nrow(cdata.LRR.invasion.num.corella)
head(cdata.LRR.invasion.num.corella)

cdata.LRR.all.species<-as.data.frame(rbind(cdata.LRR.invasion.num.corella, cdata.LRR.invasion.num.white.bryo, cdata.LRR.invasion.num.red.bryo, cdata.LRR.invasion.num.barn, cdata.LRR.invasion.num.serpulid, cdata.LRR.invasion.barn,cdata.LRR.invasion.protozoa, cdata.LRR.invasion.corella, cdata.LRR.invasion.mussel,cdata.LRR.invasion.hydroid,cdata.LRR.invasion.red.bryo, cdata.LRR.invasion.botryllid, cdata.LRR.invasion.bot.eaten,cdata.LRR.invasion.membranipora))
cdata.LRR.all.dead.species<-as.data.frame(rbind(cdata.LRR.invasion.bot.eaten,cdata.LRR.invasion.mem.dead, cdata.LRR.invasion.mem.eaten))

cdata.LRR.all.species.community<-as.data.frame(rbind(cdata.LRR.invasion.bare, cdata.LRR.invasion.num.species.no.bot, cdata.LRR.invasion.shannon.diversity.no.bot))

library(dplyr)
theme_set(theme_classic())

###
head(data_wide_invasion_comboweeks_bind)
cdata.LRR.all.species_CO2 <- cdata.LRR.all.species[cdata.LRR.all.species$LnRR.type=="CO2.overall",]
xvals_ <- cdata.LRR.all.species_CO2[with(cdata.LRR.all.species_CO2, order(estimate)), ]$species
cdata.LRR.all.species$species2<-factor(cdata.LRR.all.species$species, levels=xvals_)

plot.invasion_species_models<-ggplot(cdata.LRR.all.species, aes(x=species2,y=estimate,col=LnRR.type))+ scale_colour_manual(values=cbbPalette.blue.purple.red)
plot.invasion_species_models<-plot.invasion_species_models + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = pd)
plot.invasion_species_models<-plot.invasion_species_models+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_species_models<-plot.invasion_species_models+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_species_models<-plot.invasion_species_models+ xlab('Species') +ylab ('Effect size')
plot.invasion_species_models<-plot.invasion_species_models+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_species_models<-plot.invasion_species_models +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_species_models



cdata.LRR.all.dead.species_CO2 <- cdata.LRR.all.dead.species[cdata.LRR.all.dead.species$LnRR.type=="CO2.overall",]
xvals_ <- cdata.LRR.all.dead.species_CO2[with(cdata.LRR.all.dead.species_CO2, order(estimate)), ]$species
cdata.LRR.all.dead.species$species2<-factor(cdata.LRR.all.dead.species$species, levels=xvals_)

plot.invasion_dead.species_models<-ggplot(cdata.LRR.all.dead.species, aes(x=species2,y=estimate,col=LnRR.type))+ scale_colour_manual(values=cbbPalette.blue.purple.red)
plot.invasion_dead.species_models<-plot.invasion_dead.species_models + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = pd)
plot.invasion_dead.species_models<-plot.invasion_dead.species_models+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_dead.species_models<-plot.invasion_dead.species_models+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_dead.species_models<-plot.invasion_dead.species_models+ xlab('dead.species') +ylab ('Effect size')
plot.invasion_dead.species_models<-plot.invasion_dead.species_models+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_dead.species_models<-plot.invasion_dead.species_models +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_dead.species_models


cdata.LRR.all.species.community_CO2 <- cdata.LRR.all.species.community[cdata.LRR.all.species.community$LnRR.type=="CO2.overall",]
xvals_ <- cdata.LRR.all.species.community_CO2[with(cdata.LRR.all.species.community_CO2, order(estimate)), ]$species
cdata.LRR.all.species.community$species2<-factor(cdata.LRR.all.species.community$species, levels=xvals_)

plot.invasion_species.community_models<-ggplot(cdata.LRR.all.species.community, aes(x=species2,y=estimate,col=LnRR.type))+ scale_colour_manual(values=cbbPalette.blue.purple.red)
plot.invasion_species.community_models<-plot.invasion_species.community_models + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = pd)
plot.invasion_species.community_models<-plot.invasion_species.community_models+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_species.community_models<-plot.invasion_species.community_models+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_species.community_models<-plot.invasion_species.community_models+ xlab('species.community') +ylab ('Effect size')
plot.invasion_species.community_models<-plot.invasion_species.community_models+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_species.community_models<-plot.invasion_species.community_models +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_species.community_models











invasion.exp.data.not0<-invasion.exp.data[invasion.exp.data$Week!=0,]

##############
#Recruitment and cover together
invasion.exp.data.last.month<-rename(invasion.exp.data.last.month, c("membranipora"="% membranipora","formicula"="% formicula", "hydroid"="% hydroid", "caprellid"="% caprellid", "botryllus"="% botryllus", "didemnum"="% didemnum"))
length(invasion.exp.data.lastmonth)

data_long_invasion_rec_cover <- gather(invasion.exp.data.lastmonth, species, cover, c(8:81), factor_key=TRUE)
head(data_long_invasion_rec_cover)
View(data_long_invasion_rec_cover)
head(invasion.exp.data.lastmonth)

library(doBy)

length2 <- function (x, na.rm=FALSE) {if (na.rm) sum(!is.na(x))else length(x)}

head(data_long_invasion_rec_cover)


### This is actually useful....can't replicated within mesocosm
cdata.cover_rec_cover <- summaryBy(cover ~ species + CO2.Treatment + Invasives +Week , data=data_long_invasion_rec_cover, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length2(x)) } )
names(cdata.cover_rec_cover)[names(cdata.cover_rec_cover)=="cover.length"]<-"N"
cdata.cover_rec_cover$cover.se<-cdata.cover_rec_cover$cover.sd/sqrt(cdata.cover_rec_cover$N)
head(cdata.cover_rec_cover)

#### adding in min.10.pH here doesn't work...??! b/c need only one for each week?? Delta CO2.Treatment
### that doesn't make sense remember? 


fix(cdata.cover_rec_cover)

##### REdone April 2017

data_wide_invasion_rec_cover <- reshape(cdata.cover_rec_cover, 
                         timevar = "CO2.Treatment",
                         idvar = c("species", "Invasives", "Week"),
                         direction = "wide")

head(data_wide_invasion_rec_cover)

data_wide_invasion_rec_cover$cover.mean.Elevated<-data_wide_invasion_rec_cover$cover.mean.Elevated+0.001
data_wide_invasion_rec_cover$cover.mean.Ambient<-data_wide_invasion_rec_cover$cover.mean.Ambient+0.001


#### Calculate effect sizes
library(metafor)
data_wide_invasion_rec_cover2<-escalc(measure="ROM",m1i=cover.mean.Elevated, m2i=cover.mean.Ambient, sd1i=cover.sd.Elevated,sd2i=cover.sd.Ambient,n1i=N.Elevated,n2i=N.Ambient, data=data_wide_invasion_rec_cover,var.names=c("LRR","LRR_var"), digits=4)
head(data_wide_invasion_rec_cover2)
pd <- position_dodge(width = 0.4)
cbbPalette.all.3<- c( "#F8766D", "#00BA38", "#619CFF")
theme_set(theme_classic())
#I'm looking to reorder the x-axis by the `percent` values for category "bb"


data_wide_invasion_rec_cover2.none <- data_wide_invasion_rec_cover2[data_wide_invasion_rec_cover2$Invasives == 'None', ]
xvals <- data_wide_invasion_rec_cover2.none[with(data_wide_invasion_rec_cover2.none, order(LRR)), ]$species
data_wide_invasion_rec_cover2$species2<-factor(data_wide_invasion_rec_cover2$species, levels=xvals)
head(data_wide_invasion_rec_cover2)

View(data_wide_invasion_rec_cover2)


data_wide_invasion_rec_cover2[,12][data_wide_invasion_rec_cover2[,13]==0] <- NA
data_wide_invasion_rec_cover2[,13][data_wide_invasion_rec_cover2[,13]==0] <- NA


library(ggplot2)

## LRRvar
plot.invasion_rec_cover<-ggplot(data_wide_invasion_rec_cover2, aes(x=species2,y=LRR,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.invasion_rec_cover<-plot.invasion_rec_cover + geom_point(size=4, position = pd)+geom_errorbar(aes(ymin = (LRR - (LRR_var)), ymax = (LRR + (LRR_var))), width = 0.1, position = pd)
plot.invasion_rec_cover<-plot.invasion_rec_cover+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_rec_cover<-plot.invasion_rec_cover+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high " ~CO[2]))
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.invasion_rec_cover  + theme(legend.position="none")



##### confidence interval 95% 
plot.invasion_rec_cover<-ggplot(data_wide_invasion_rec_cover2, aes(x=species2,y=LRR,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.invasion_rec_cover<-plot.invasion_rec_cover + geom_point(size=4, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.invasion_rec_cover<-plot.invasion_rec_cover+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_rec_cover<-plot.invasion_rec_cover+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high " ~CO[2]))
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.invasion_rec_cover


#### model.invasion for this framework needs an estimate and confidence intervals... 

most.cover<-data_wide_invasion_rec_cover2[data_wide_invasion_rec_cover2$cover.mean.Ambient>2, ]
levels(most.cover$species)
View(most.cover)


####### 
############## model.invasion 
##### anova, model.invasion, stats
#### do we want to include min.10.pH in the model.invasion?? NO. Would have to be delta pH
head(data_wide_invasion_rec_cover2)

View(data_wide_invasion_rec_cover2)


###blue only
###
plot.rec_cover<-ggplot(cdata.LRR.all.species, aes(x=species2,y=estimate,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3.blue)
plot.rec_cover<-plot.rec_cover + geom_point(size=5, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = pd)
plot.rec_cover<-plot.rec_cover+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)#+coord_flip()
plot.rec_cover<-plot.rec_cover+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.rec_cover<-plot.rec_cover+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high" ~CO[2]))
plot.rec_cover<-plot.rec_cover+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.rec_cover<-plot.rec_cover +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.rec_cover<-plot.rec_cover +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.rec_cover <- plot.rec_cover + theme(legend.position="none")
plot.rec_cover 

##corrdinate flip?? 


plot.recdead<-ggplot(cdata.LRR.all.dead.species, aes(x=species,y=estimate,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.recdead<-plot.recdead + geom_point(size=5, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = pd)
plot.recdead<-plot.recdead+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)#+coord_flip()
plot.recdead<-plot.recdead+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.recdead<-plot.recdead+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high " ~CO[2]))
plot.recdead<-plot.recdead+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.recdead<-plot.recdead +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.recdead<-plot.recdead +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.recdead <- plot.recdead + theme(legend.position="none")
plot.recdead 


plot.rec_disporella<-ggplot(cdata.LRR.morphology, aes(x=species,y=estimate,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.rec_disporella<-plot.rec_disporella + geom_point(size=5, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = pd)
plot.rec_disporella<-plot.rec_disporella+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)#+coord_flip()
plot.rec_disporella<-plot.rec_disporella+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.rec_disporella<-plot.rec_disporella+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high " ~CO[2]))
plot.rec_disporella<-plot.rec_disporella+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.rec_disporella<-plot.rec_disporella +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.rec_disporella<-plot.rec_disporella +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.rec_disporella <- plot.rec_disporella + theme(legend.position="none")
plot.rec_disporella 


plot.rec_community<-ggplot(cdata.LRR.community, aes(x=species,y=estimate,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.rec_community<-plot.rec_community + geom_point(size=5, position = pd)+geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = pd)
plot.rec_community<-plot.rec_community+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)#+coord_flip()
plot.rec_community<-plot.rec_community+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.rec_community<-plot.rec_community+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high" ~CO[2]))
plot.rec_community<-plot.rec_community+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.rec_community<-plot.rec_community +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.rec_community<-plot.rec_community +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.rec_community + theme(legend.position="none")




















#In lm world: (just one week) ... but confused about the CO2.Treatment level - I need an estaimet for each food level...where is that? 
####### model.invasion that could go into this framework
M1.hydroid<-lm(hydroid~Invasives*CO2.Treatment - 1,data=invasion.exp.data.16)
summary(M1.hydroid)
anova(M1.hydroid,  L=c(1,-1, -1))
lm.beta(M1.hydroid)

coef(M1.hydroid)
M1.coef<-as.data.frame(cbind(confint(M1.hydroid), coef(M1.hydroid)))
colnames(M1.coef)<-c( "lower", "upper", "beta")


M1.hydroid$effects

##beta estimate
install.packages("lm.beta")
library(lm.beta)

lm.beta(M1.hydroid)
# standardize
lm.D9.beta <- lm.beta(M1.hydroid)
print(lm.D9.beta)
summary(lm.D9.beta)
coef(lm.D9.beta)



lmer.bottohyd.weeks<- lmer(bottohyd~ min.10.pH*factor(Invasives)*Week + (Week|Mesocosm), data = invasion.exp.data.not0)
summary(lmer.bottohyd.weeks)
######
#In model.invasion world but binomial. 

Mod.hydroid.glm.binom<- glm(cbind(hydroid, 100-hydroid)~CO2.Treatment*Invasives-1, data = invasion.exp.data.not0, family=binomial(link="logit"))
coef(Mod.hydroid.glm.binom)
confint(Mod.hydroid.glm.binom)


#### separate model.invasion for each species? Exactly confidence initervals 
### YEs... need to figure out is this the best model.invasion - or maybe should I use the bootstrapped data?? 
#probably not good to do all at once .. b/c ot independent... 
#model.invasion.all<-rma.mv(LRR,LRR_var,mods=~factor(Invasives)*species - 1,data=data_wide_invasion_rec_cover2)
#summary(model.invasion.all)

model.invasion.hydroid<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="% hydroid"))
summary(model.invasion.hydroid)
anova(model.invasion.hydroid,  L=c(1,-1, -1))

model.invasion.nudi.eggs<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="# nudibranch eggs"))
summary(model.invasion.corella)
anova(model.invasion.corella,  L=c(1, -1))

model.invasion.corella<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="# corella"))
summary(model.invasion.corella)
anova(model.invasion.corella,  L=c(1, -1))

model.invasion.corella<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="# corella"))
summary(model.invasion.corella)
anova(model.invasion.corella,  L=c(1, -1))


mod.overall.best.growth <- rma.mv(LRR, V.growth, mods = ~factor(Food.supply) - 1,random=list(~factor(Food.supply)|Paper_no, ~Unit|Paper_no), data=meta.food.growth2, struct=c("UN", "UN"))
summary(mod.overall.best.growth)
anova(mod.overall.best.growth,  L=c(1,-1))


head(data_wide_invasion_rec_cover2)
model.invasion.corella<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="# corella"))
confidence.from.model.invasion.corella<-cbind(model.invasion.corella$ci.ub, model.invasion.corella$ci.lb, "# corella")
colnames(confidence.from.model.invasion.corella)<- c("upper","lower","species")


model.invasion.corella<-rma.mv(LRR,LRR_var,mods=~factor(Invasives) - 1,data=data_wide_invasion_rec_cover2,subset=(species=="# corella"))
confidence.from.model.invasion.corella<-cbind(model.invasion.corella$ci.ub, model.invasion.corella$ci.lb, "# corella")
colnames(confidence.from.model.invasion.corella)<- c("upper","lower","species")




##### confidence interval 95% 
plot.invasion_rec_cover<-ggplot(data_wide_invasion_rec_cover2, aes(x=species2,y=LRR,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.invasion_rec_cover<-plot.invasion_rec_cover + geom_point(size=4, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.invasion_rec_cover<-plot.invasion_rec_cover+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_rec_cover<-plot.invasion_rec_cover+ xlab('Recruited species') +ylab(expression("Abundance ratio: effect of high " ~CO[2]))
plot.invasion_rec_cover<-plot.invasion_rec_cover+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_rec_cover<-plot.invasion_rec_cover +theme(legend.title = element_text(colour="black", size=16), legend.position=c(.10,.85))
plot.invasion_rec_cover









### other model.invasion options: binomial or poisson... not sure which works?? a bit odd
glm.serpulid<-glm(formula = (cover.mean.Elevated/cover.mean.Ambient) ~Invasives-1, data=data_wide_invasion_rec_cover2,subset=(species=="# serpulid"), family=poisson(link="logit"))
summary(glm.serpulid)
#normal?
qqp(data_wide_invasion_rec_cover2$LRR, "norm")
model.invasion.serpulid.lm<-lm(log(cover.mean.Elevated/cover.mean.Ambient) ~ Invasives -1, data=data_wide_invasion_rec_cover2,subset=(species=="# serpulid"))
summary(model.invasion.serpulid.lm)









####### do low quality as well

data_wide_invasion_food_rec2<-escalc(measure="ROM",m1i=cover.mean.High, m2i=cover.mean.None, sd1i=cover.sd.High,sd2i=cover.sd.None,n1i=N.High,n2i=N.None, data=data_wide_invasion_food_rec,var.names=c("LRR","LRR_var"), digits=4)

data_wide_invasion_food_rec2.none <- data_wide_invasion_food_rec2[data_wide_invasion_food_rec2$CO2.Treatment == 'Ambient', ]
xvalsfood_ <- data_wide_invasion_food_rec2.none[with(data_wide_invasion_food_rec2.none, order(LRR)), ]$species
data_wide_invasion_food_rec2$species2<-factor(data_wide_invasion_food_rec2$species, levels=xvalsfood_)



plot.not.mod.3.growth_rec<-ggplot(data_wide_invasion_food_rec2, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.red)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec+ xlab('Recruitment Species') +ylab ('Effect of Food')
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec<-plot.not.mod.3.growth_rec +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec



#######3Low quality

data_wide_invasion_food_rec3<-escalc(measure="ROM",m1i=cover.mean.Low, m2i=cover.mean.None, sd1i=cover.sd.Low,sd2i=cover.sd.None,n1i=N.Low,n2i=N.None, data=data_wide_invasion_food_rec,var.names=c("LRR","LRR_var"), digits=4)

data_wide_invasion_food_rec3.none <- data_wide_invasion_food_rec3[data_wide_invasion_food_rec3$CO2.Treatment == 'Ambient', ]
xvalsfood_ <- data_wide_invasion_food_rec3.none[with(data_wide_invasion_food_rec3.none, order(LRR)), ]$species
data_wide_invasion_food_rec3$species2<-factor(data_wide_invasion_food_rec3$species, levels=xvalsfood_)



plot.not.mod.3.growth_rec_low<-ggplot(data_wide_invasion_food_rec3, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.green)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low+ xlab('Recruitment Species') +ylab ('Effect of Food')
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec_low<-plot.not.mod.3.growth_rec_low +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec_low





##### recruitment both 

data_wide_invasion_food_rec_4<-cbind(data_wide_invasion_food_rec2, data_wide_invasion_food_rec3$LRR, data_wide_invasion_food_rec3$LRR_var)
head(data_wide_invasion_food_rec_4)
data_wide_invasion_food_rec_5 <- gather(data_wide_invasion_food_rec_4, food.level, LRR_rec_both, c(15,18), factor_key=TRUE)
head(data_wide_invasion_food_rec_5)
data_wide_invasion_food_rec_6 <-gather(data_wide_invasion_food_rec_5, food.level.2, LRR_var_rec_both, c(15,17), factor_key=TRUE)
head(data_wide_invasion_food_rec_6)

levels(data_wide_invasion_food_rec_6$food.level)[levels(data_wide_invasion_food_rec_6$food.level)=="LRR"] <- "High quality food"
levels(data_wide_invasion_food_rec_6$food.level)[levels(data_wide_invasion_food_rec_6$food.level)=="data_wide_invasion_food_rec3$LRR"] <- "Low quality food"



data_wide_invasion_food_rec_6.none <- data_wide_invasion_food_rec_6[data_wide_invasion_food_rec_6$CO2.Treatment == 'Ambient', ]
xvals_food_rec_6 <- data_wide_invasion_food_rec_6.none[with(data_wide_invasion_food_rec_6.none, order(LRR_rec_both)), ]$species
data_wide_invasion_food_rec_6$species3<-factor(data_wide_invasion_food_rec_6$species, levels=xvals_food_rec_6)
head(data_wide_invasion_food_rec_6)



plot.not.mod.3.growth_rec_both<-ggplot(data_wide_invasion_food_rec_6, aes(x=species3,y=LRR_rec_both,col=food.level, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.red.green)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR_rec_both - 1.96*sqrt(LRR_var_rec_both)), ymax = (LRR_rec_both + 1.96*sqrt(LRR_var_rec_both))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both+ xlab('Recruitment of Species') +ylab ('Effect of food')
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec_both<-plot.not.mod.3.growth_rec_both +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec_both


#### The order is not working here.... 

################# Recruitment in mesocosm

head(invasion.exp.data.inventory.12)
invasion.exp.data.inventory.12<-rename(invasion.exp.data.inventory.12, c("nudis"="nudibranch", "red.bryo"="schizoporella", "mem.bryo"="membranipora", "nudi.eggs"="nudibranch.eggs", "yellow.bryo"="cribrilina", "Disporella.total"="disporella"))

data_long_invasion_rec_meso <- gather(invasion.exp.data.inventory.12, species, cover, c(8:10, 12:17, 19,21:26, 34), factor_key=TRUE)
head(data_long_invasion_rec_meso)



cdata.cover_rec_meso <- summaryBy(cover ~ species + CO2.Treatment + Invasives, data=data_long_invasion_rec_meso, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.cover_rec_meso)[names(cdata.cover_rec_meso)=="cover.length"]<-"N"
cdata.cover_rec_meso$cover.se<-cdata.cover_rec_meso$cover.sd/sqrt(cdata.cover_rec_meso$N)
head(cdata.cover_rec_meso)
fix(cdata.cover_rec_meso)

data_wide_invasion_rec_meso <- reshape(cdata.cover_rec_meso, 
                         timevar = "CO2.Treatment",
                         idvar = c("species", "Invasives"),
                         direction = "wide")

data_wide_invasion_food_rec_meso <- reshape(cdata.cover_rec_meso, 
                              timevar = "Invasives",
                              idvar = c("species", "CO2.Treatment"),
                              direction = "wide")

head(data_wide_invasion_food)

data_wide_invasion_rec_meso$cover.mean.Elevated<-data_wide_invasion_rec_meso$cover.mean.Elevated+0.001
data_wide_invasion_rec_meso$cover.mean.Ambient<-data_wide_invasion_rec_meso$cover.mean.Ambient+0.001


#### Calculate effect sizes
library(metafor)
data_wide_invasion_rec_meso2<-escalc(measure="ROM",m1i=cover.mean.Elevated, m2i=cover.mean.Ambient, sd1i=cover.sd.Elevated,sd2i=cover.sd.Ambient,n1i=N.Elevated,n2i=N.Ambient, data=data_wide_invasion_rec_meso,var.names=c("LRR","LRR_var"), digits=4)
View(data_wide_invasion_rec_meso2)
pd <- position_dodge(width = 0.4)
cbbPalette.all.3<- c( "#F8766D", "#00BA38", "#619CFF")
theme_set(theme_classic())

data_wide_invasion_rec_meso2.none <- data_wide_invasion_rec_meso2[data_wide_invasion_rec_meso2$Invasives == 'None', ]
xvals_meso <- data_wide_invasion_rec_meso2.none[with(data_wide_invasion_rec_meso2.none, order(LRR)), ]$species
data_wide_invasion_rec_meso2$species2<-factor(data_wide_invasion_rec_meso2$species, levels=xvals_meso)



plot.invasion_rec_meso<-ggplot(data_wide_invasion_rec_meso2, aes(x=species2,y=LRR,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.invasion_rec_meso<-plot.invasion_rec_meso + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.invasion_rec_meso<-plot.invasion_rec_meso+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion_rec_meso<-plot.invasion_rec_meso+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion_rec_meso<-plot.invasion_rec_meso+ xlab('Species Recruitment') +ylab ('Effect of CO2.Treatment')
plot.invasion_rec_meso<-plot.invasion_rec_meso+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion_rec_meso<-plot.invasion_rec_meso +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion_rec_meso





data_wide_invasion_food_rec_meso2<-escalc(measure="ROM",m1i=cover.mean.High, m2i=cover.mean.None, sd1i=cover.sd.High,sd2i=cover.sd.None,n1i=N.High,n2i=N.None, data=data_wide_invasion_food_rec_meso,var.names=c("LRR","LRR_var"), digits=4)


cbbPalette.black<-c("#000000", "#000000")


data_wide_invasion_food_rec_meso2.none <- data_wide_invasion_food_rec_meso2[data_wide_invasion_food_rec_meso2$CO2.Treatment == 'Ambient', ]
xvalsfood_meso <- data_wide_invasion_food_rec_meso2.none[with(data_wide_invasion_food_rec_meso2.none, order(LRR)), ]$species
data_wide_invasion_food_rec_meso2$species2<-factor(data_wide_invasion_food_rec_meso2$species, levels=xvalsfood_meso)


plot.not.mod.3.growth_rec_meso<-ggplot(data_wide_invasion_food_rec_meso2, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.red)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso+ xlab('Percent cover of Species') +ylab ('Effect of Food')
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec_meso<-plot.not.mod.3.growth_rec_meso +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec_meso


########Low food
data_wide_invasion_food_rec_meso3<-escalc(measure="ROM",m1i=cover.mean.Low, m2i=cover.mean.None, sd1i=cover.sd.Low,sd2i=cover.sd.None,n1i=N.Low,n2i=N.None, data=data_wide_invasion_food_rec_meso,var.names=c("LRR","LRR_var"), digits=4)


cbbPalette.black<-c("#000000", "#000000")


data_wide_invasion_food_rec_meso3.none <- data_wide_invasion_food_rec_meso3[data_wide_invasion_food_rec_meso3$CO2.Treatment == 'Ambient', ]
xvalsfood_meso <- data_wide_invasion_food_rec_meso3.none[with(data_wide_invasion_food_rec_meso3.none, order(LRR)), ]$species
data_wide_invasion_food_rec_meso3$species2<-factor(data_wide_invasion_food_rec_meso3$species, levels=xvalsfood_meso)


plot.not.mod.3.growth_rec_meso_low<-ggplot(data_wide_invasion_food_rec_meso3, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.green)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low+ xlab('# of individuals') +ylab ('Effect of Food')
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec_meso_low<-plot.not.mod.3.growth_rec_meso_low +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec_meso_low


##### recruitment mesocosms both 

data_wide_invasion_food_rec_meso_4<-cbind(data_wide_invasion_food_rec_meso2, data_wide_invasion_food_rec_meso3$LRR, data_wide_invasion_food_rec_meso3$LRR_var)
head(data_wide_invasion_food_rec_meso_4)
data_wide_invasion_food_rec_meso_5 <- gather(data_wide_invasion_food_rec_meso_4, food.level, LRR_rec_meso_both, c(15,18), factor_key=TRUE)
head(data_wide_invasion_food_rec_meso_5)
data_wide_invasion_food_rec_meso_6 <-gather(data_wide_invasion_food_rec_meso_5, food.level.2, LRR_var_rec_meso_both, c(15,17), factor_key=TRUE)
head(data_wide_invasion_food_rec_meso_6)

levels(data_wide_invasion_food_rec_meso_6$food.level)[levels(data_wide_invasion_food_rec_meso_6$food.level)=="LRR"] <- "High quality food"
levels(data_wide_invasion_food_rec_meso_6$food.level)[levels(data_wide_invasion_food_rec_meso_6$food.level)=="data_wide_invasion_food_rec_meso3$LRR"] <- "Low quality food"



data_wide_invasion_food_rec_meso_6.none <- data_wide_invasion_food_rec_meso_6[data_wide_invasion_food_rec_meso_6$CO2.Treatment == 'Ambient', ]
xvals_food_rec_meso_6 <- data_wide_invasion_food_rec_meso_6.none[with(data_wide_invasion_food_rec_meso_6.none, order(LRR_rec_meso_both)), ]$species
data_wide_invasion_food_rec_meso_6$species3<-factor(data_wide_invasion_food_rec_meso_6$species, levels=xvals_food_rec_meso_6)
head(data_wide_invasion_food_rec_meso_6)



plot.not.mod.3.growth_rec_meso_both<-ggplot(data_wide_invasion_food_rec_meso_6, aes(x=species3,y=LRR_rec_meso_both,col=food.level, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.red.green)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR_rec_meso_both - 1.96*sqrt(LRR_var_rec_meso_both)), ymax = (LRR_rec_meso_both + 1.96*sqrt(LRR_var_rec_meso_both))), width = 0.1, position = pd)
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both+ xlab('Recruitment of Species') +ylab ('Effect of food')
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_rec_meso_both<-plot.not.mod.3.growth_rec_meso_both +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_rec_meso_both



###################### Community level measures

head(invasion.exp.data.16)
data_long_invasion_comm <- gather(invasion.exp.data.16, species, cover, c(44, 73, 74, 75,78 , 82,83), factor_key=TRUE)
head(data_long_invasion_comm)

cdata.cover.comm <- summaryBy(cover ~ species + CO2.Treatment + Invasives, data=data_long_invasion_comm, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.cover.comm)[names(cdata.cover.comm)=="cover.length"]<-"N"
cdata.cover.comm$cover.se<-cdata.cover.comm$cover.sd/sqrt(cdata.cover.comm$N)
head(cdata.cover.comm)
fix(cdata.cover.comm)

data_wide_invasion_comm <- reshape(cdata.cover.comm, 
                     timevar = "CO2.Treatment",
                     idvar = c("species", "Invasives"),
                     direction = "wide")

data_wide_invasion_comm_food <- reshape(cdata.cover.comm, 
                          timevar = "Invasives",
                          idvar = c("species", "CO2.Treatment"),
                          direction = "wide")

head(data_wide_invasion_comm_food)

#### Calculate effect sizes
library(metafor)
data_wide_invasion_comm2<-escalc(measure="ROM",m1i=cover.mean.Elevated, m2i=cover.mean.Ambient, sd1i=cover.sd.Elevated,sd2i=cover.sd.Ambient,n1i=N.Elevated,n2i=N.Ambient, data=data_wide_invasion_comm,var.names=c("LRR","LRR_var"), digits=4)



length(unique(data_wide_invasion_comm2$species))


pd <- position_dodge(width = 0.4)
cbbPalette.all.3<- c( "#F8766D", "#00BA38", "#619CFF")
theme_set(theme_classic())

data_wide_invasion_comm2.none <- data_wide_invasion_comm2[data_wide_invasion_comm2$Invasives == 'None', ]
xvals_ <- data_wide_invasion_comm2.none[with(data_wide_invasion_comm2.none, order(LRR)), ]$species
data_wide_invasion_comm2$species2<-factor(data_wide_invasion_comm2$species, levels=xvals_)


plot.invasion.comm<-ggplot(data_wide_invasion_comm2, aes(x=species2,y=LRR,col=Invasives))+ scale_colour_manual(values=cbbPalette.all.3)
plot.invasion.comm<-plot.invasion.comm + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.invasion.comm<-plot.invasion.comm+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.invasion.comm<-plot.invasion.comm+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.invasion.comm<-plot.invasion.comm+ xlab('Perecent cover') +ylab ('Effect of CO2.Treatment')
plot.invasion.comm<-plot.invasion.comm+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.invasion.comm<-plot.invasion.comm +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.invasion.comm





data_wide_invasion_comm_food2<-escalc(measure="ROM",m1i=cover.mean.High, m2i=cover.mean.None, sd1i=cover.sd.High,sd2i=cover.sd.None,n1i=N.High,n2i=N.None, data=data_wide_invasion_comm_food,var.names=c("LRR","LRR_var"), digits=4)

cbbPalette.all.3<- c( "#F8766D", "#00BA38", "#619CFF")
cbbPalette.red<-c("#F8766D", "#F8766D")
cbbPalette.green<-c("#00BA38", "#00BA38")


data_wide_invasion_comm_food2.none <- data_wide_invasion_comm_food2[data_wide_invasion_comm_food2$CO2.Treatment == 'Ambient', ]
xvals_food_ <- data_wide_invasion_comm_food2.none[with(data_wide_invasion_comm_food2.none, order(LRR)), ]$species
data_wide_invasion_comm_food2$species2<-factor(data_wide_invasion_comm_food2$species, levels=xvals_food_)


plot.not.mod.3.growth.comm<-ggplot(data_wide_invasion_comm_food2, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.red)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm+ xlab('Percent cover of Species') +ylab ('Effect of high quality food')
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth.comm<-plot.not.mod.3.growth.comm +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth.comm

##### Low food
head(data_wide_invasion_comm_food)
data_wide_invasion_comm_food3<-escalc(measure="ROM",m1i=cover.mean.Low, m2i=cover.mean.None, sd1i=cover.sd.Low,sd2i=cover.sd.None,n1i=N.Low,n2i=N.None, data=data_wide_invasion_comm_food,var.names=c("LRR","LRR_var"), digits=4)

data_wide_invasion_comm_food3.none <- data_wide_invasion_comm_food3[data_wide_invasion_comm_food3$CO2.Treatment == 'Ambient', ]
xvals_food_ <- data_wide_invasion_comm_food3.none[with(data_wide_invasion_comm_food3.none, order(LRR)), ]$species
data_wide_invasion_comm_food3$species2<-factor(data_wide_invasion_comm_food3$species, levels=xvals_food_)



plot.not.mod.3.growth_low_comm<-ggplot(data_wide_invasion_comm_food3, aes(x=species2,y=LRR,col=CO2.Treatment, shape=CO2.Treatment))+  scale_colour_manual(values=cbbPalette.green)+ scale_shape_manual(values=c(19,2))
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm + geom_point(size=3, position = pd)+geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1, position = pd)
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm+ xlab('Percent cover of Species') +ylab ('Effect of low quality food')
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm+theme(axis.text.x = element_text(size = 16, colour = 'black', angle=90, hjust = 1)) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.3.growth_low_comm<-plot.not.mod.3.growth_low_comm +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.3.growth_low_comm





























########## Best of all three

data_long_invasion_together <- gather(invasion.exp.data.16, species, cover, c(8:14,20,21,26, 31), factor_key=TRUE)
head(data_long_invasion)

cdata.cover <- summaryBy(cover ~ species + CO2.Treatment + Invasives, data=data_long_invasion, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.cover)[names(cdata.cover)=="cover.length"]<-"N"
cdata.cover$cover.se<-cdata.cover$cover.sd/sqrt(cdata.cover$N)
head(cdata.cover)
fix(cdata.cover)

data_wide_invasion <- reshape(cdata.cover, 
                     timevar = "CO2.Treatment",
                     idvar = c("species", "Invasives"),
                     direction = "wide")


############## CHOSEN model.invasion
# food supply in every paper so needs to accounted for ... and then unit in most papers, CO2.Treatment not in that many and already accounted for. 
#CO2.Treatment already accounted for sampling errors correlated with the V.growth part

mod.overall.best.food.growth <- rma.mv(LRR, LRR_var, mods = ~factor(CO2.Treatment_level) - 1,random=list(~factor(CO2.Treatment_level)|Paper_no,~Unit|Paper_no), data=data_wide_invasion, struct=c("UN", "CS"))
summary(mod.overall.best.food.growth)
profile(mod.overall.best.food.growth)
anova(mod.overall.best.food.growth,  L=c(1,-1))



qqnorm(residuals(mod.overall.best.food.growth,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.overall.best.food.growth,type="pearson"),col="red")
plot(fitted(mod.overall.best.food.growth) ~ residuals(mod.overall.best.food.growth))




###### Checking biases
funnel(mod.overall.best.food.growth,ylim=c(1:6,by=2),yaxis="seinv",level=c(90, 95, 99),ylab="Precision (1/SE)" ,shade=c("white", "gray", "darkgray"), refline=0)
funnel(mod.1,ylim=c(1:6,by=2),yaxis="seinv",level=c(90, 95, 99),ylab="Precision (1/SE)" ,shade=c("white", "gray", "darkgray"), refline=0)

#Publication bias with more basic model.invasion (some can't handle ram.mv)
mod.basic.growth<-rma(LRR ~ factor(Food.supply),LRR_var,data=data_wide_invasion)
funnel(mod.basic)
regtest(mod.basic.growth,model.invasion="rma",predictor="sei")
funnel(trimfill(mod.basic, side="right"))
##### The method can be used to estimate the number of studies missing from a meta-analysis due to the suppression of the most extreme results on one side of the funnel plot.
baujat(mod.basic.growth)
#### Baujat et al. (2002) proposed a diagnostic plot to detect sources of heterogeneity in meta-analytic data. 
#The plot shows the contribution of each study to the overall Q-test statistic for heterogeneity on the horizontal axis 
#versus the influence of each study (defined as the standardized squared difference between the overall estimate based 
#on a fixed-effects model.invasion with and without the ith study included in the model.invasion) on the vertical axis. An example of such a plot is shown below.

### calculate influence diagnostics
inf <- influence(mod.basic)
### plot the influence diagnostics
plot(inf, layout=c(8,1))

########## Three studies have rather large residuals and may be considered outliers. 
#Only two of those studies actually have a strong influence on the results (as reflected, for example, 
#in their Cook's distances). Removal of these studies would reduce the amount of heterogeneity quite a bit 
#and increase the precision of the estimated average outcome (i.e., r-to-z transformed correlation) from the 
#random-effects model.invasion. However, instead of just removing those studies, one should examine them in detail to 
#determine what the reason may be for their unusual results. Outliers and influential cases can actually reveal 
#patterns that may lead to new insights about study characteristics that could be acting as potential moderators (Light & Pillemer, 1984).


###Sensitivity analysis ## analgous to drop1
# RULE:If rstandard >3 AND hatvalue >2 times average of hatvalues, 
# run analysis with those cases deleted to test for sensitivity.
rs.abund.wood <- rstandard(mod.overall.best.food.growth)
hat.abund.wood <- hatvalues(mod.overall.best.food.growth)/mean(hatvalues(mod.overall.best.food.growth))
plot(hat.abund.wood, rs.abund.wood$resid, ylim = c(-4.0,4))
text(hat.abund.wood, rs.abund.wood$resid, labels = data_wide_invasion$Paper_no, cex= 1, pos = 2)
abline(h = -3)
abline(h = 3)
abline( v = 2)
#### should compare without Crook, Comeau, Hettinger outliers model.invasion fit







############# GGPLOT using summary data for each Food.supply level
library(ggplot2)

######## PLOTTING WITH A model.invasion


summary(mod.1.growth)
head(mod.1.growth)
library(plyr)
library(dplyr)


cdata.LRR.growth.byfood<-ddply(data_wide_invasion,~CO2.Treatment_level,summarise,mean=mean(LRR))
head(cdata.LRR.growth.byfood)
cdata.LRR.growth.byfood_var<-ddply(data_wide_invasion,~CO2.Treatment_level,summarise,mean=mean(LRR_var))
head(cdata.LRR.growth.byfood_var)
cdata.LRR.growth.byfood$var<-cdata.LRR.growth.byfood_var$mean
head(cdata.LRR.growth.byfood)


library(ggplot2)
##theme classic... 
theme_set(theme_classic())




plot.mod.overall.best.food.growth <-ggplot(cdata.LRR.growth.byfood,aes(x=CO2.Treatment_level,y=mod.overall.best.food.growth $b,ymax=(mod.overall.best.food.growth $ci.ub),ymin=(mod.overall.best.food.growth $ci.lb),size=2, col=CO2.Treatment_level, shape=CO2.Treatment_level))+ scale_colour_manual(values=cbbPalette.blue.only)+ scale_shape_manual(values=c(19,1))
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth +geom_pointrange(size=1)
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth  #+coord_flip()
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth +geom_hline(aes(x=0, yintercept=0), lty=2,size=1) + ylim(-1.5, 1.5)
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth +theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth + xlab('') +ylab ('LnRR Growth')
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth +theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.overall.best.food.growth <-plot.mod.overall.best.food.growth +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.mod.overall.best.food.growth


#black all around, same shapes
cbbPalette.black<-c("#000000", "#000000")

plot.mod.overall.best.food.growth2 <-ggplot(cdata.LRR.growth.byfood,aes(x=CO2.Treatment_level,y=mod.overall.best.food.growth $b,ymax=(mod.overall.best.food.growth $ci.ub),ymin=(mod.overall.best.food.growth $ci.lb),size=2, col=CO2.Treatment_level, shape=CO2.Treatment_level))+ scale_colour_manual(values=cbbPalette.black)+ scale_shape_manual(values=c(19,19))
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 +geom_pointrange(size=1)
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2  #+coord_flip()
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 +geom_hline(aes(x=0, yintercept=0), lty=2,size=1) + ylim(-1.3, 1.3)
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 +theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 + xlab('') +ylab ('LnRR effect of food')
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 +theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.overall.best.food.growth2 <-plot.mod.overall.best.food.growth2 +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.mod.overall.best.food.growth2





plot.mod.overall.best.food.growth3 <-ggplot(cdata.LRR.growth.byfood,aes(x=CO2.Treatment_level,y=mod.overall.best.food.growth $b,ymax=(mod.overall.best.food.growth $ci.ub),ymin=(mod.overall.best.food.growth $ci.lb),size=2, col=CO2.Treatment_level, shape=CO2.Treatment_level))+ scale_colour_manual(values=cbbPalette.blue.only)+ scale_shape_manual(values=c(19,19))
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 +geom_pointrange(size=1)
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3  #+coord_flip()
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 +geom_hline(aes(x=0, yintercept=0), lty=2,size=1) + ylim(-1.5, 1.5)
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 +theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 + xlab('') +ylab ('LnRR Growth')
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 +theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.overall.best.food.growth3 <-plot.mod.overall.best.food.growth3 +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.mod.overall.best.food.growth3



### plot all together.... 
require(cowplot)
theme_set(theme_classic())


### option 1 with colours
plot_grid(plot.mod.overall.best.food,plot.mod.overall.best, plot.mod.overall.best.food.growth,plot.mod.overall.best.growth, ncol=2, align='h', labels=c('(a)', '(b)', '(c)', '(d)', label_size=12))
### option 2 without colours
plot_grid(plot.mod.overall.best.food2,plot.mod.overall.best, plot.mod.overall.best.food.growth2,plot.mod.overall.best.growth2, ncol=2, align='h', labels=c('(a)', '(b)', '(c)', '(d)', label_size=12))
#### all 8 panels
plot_grid(plot.LRR.CO2.Treatment.calc.byfood,plot.LRR.CO2.Treatment.growth.byfood, plot.mod.overall.best.food2, plot.mod.overall.best.food.growth2,plot.LRR.CO2.Treatment.calc, plot.LRR.CO2.Treatment,plot.mod.overall.best,plot.mod.overall.best.growth2,ncol=2, align='h', labels=c('(a)', '(b)', '(c)', '(d)','(e)', '(f)', '(g)', '(h)', label_size=12))
### option 3 with colours
plot_grid(plot.mod.overall.best.food3,plot.mod.overall.best3, plot.mod.overall.best.food.growth3,plot.mod.overall.best.growth3, ncol=2, align='h', labels=c('(a)', '(b)', '(c)', '(d)', label_size=12))





# without the by delta CO2.Treatment # 6 panels Fig. 2 order redone
plot_grid(plot.mod.overall.best2,plot.mod.overall.best.growth2, plot.mod.overall.best.food2, plot.mod.overall.best.food.growth2,plot.LRR.CO2.Treatment.calc.byfood,plot.LRR.CO2.Treatment.growth.byfood, ncol=2, align='h', labels=c('(a)', '(b)', '(c)', '(d)','(e)', '(f)', '(g)', '(h)', label_size=12))





data_wide_invasion












#### I don't think this is an appropriate comparison b/c it's not PAIRED ... if we could pair the data then we might see a trend but for now  this counts every one as equal which is not true. 
err.growfood <- stats::predict(mod.overall.best.food.growth.CO2.Treatment, newdata=data_wide_invasion, se = TRUE)
data_wide_invasion <- as.data.frame(cbind(data_wide_invasion, err.growfood$ci.ub, err.growfood$ci.lb, err.growfood$pred))
head(err.growfood)
length(err.growfood$ci.ub)

### confidence bands not intervals
## try this: geom_ribbon(aes(ymin=lower, ymax=upper,alpha=0.05))
View(data_wide_invasion)


#### I think this is the right one... predicted fitted line, upper and lower CIs from the model.invasion. 

plot.LRR.CO2.Treatment.growth.byfood<- ggplot(data_wide_invasion, aes(x=data_wide_invasion$CO2.Treatment, y=data_wide_invasion$LRR))+ geom_point(aes(pch=CO2.Treatment_level), size=5) +scale_shape_manual(values=c(1,19), name="CO2.Treatment level")  + guides(fill=FALSE) + geom_ribbon(aes(ymin=err.growfood$ci.lb, ymax=err.growfood$ci.ub), alpha=0.15, data=data_wide_invasion)+geom_smooth(aes(y=err.growfood$pred), colour="black", size=1.5, se=FALSE, data=data_wide_invasion)
plot.LRR.CO2.Treatment.growth.byfood<- plot.LRR.CO2.Treatment.growth.byfood + theme_bw() + xlab(expression("" ~ CO[2])) + ylab(expression("LnRR effect of food"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))+geom_hline(aes(x=400, yintercept=0), lty=2,size=1)+ylim(-1,3)
plot.LRR.CO2.Treatment.growth.byfood<- plot.LRR.CO2.Treatment.growth.byfood + theme( legend.text = element_text(colour="black", size = 12))+ theme(legend.title = element_text(colour="black", size=12), legend.position=c(.80,.85))+ guides(shape=guide_legend(override.aes=list(size=5, linetype=0, fill="white")))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"))
plot.LRR.CO2.Treatment.growth.byfood<-plot.LRR.CO2.Treatment.growth.byfood + theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.LRR.CO2.Treatment.growth.byfood<-plot.LRR.CO2.Treatment.growth.byfood+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))+ theme(legend.position="none")
plot.LRR.CO2.Treatment.growth.byfood

?geom_line
?geom_smooth

mod.overall.best.food.growth.CO2.Treatment<- rma.mv(LRR, LRR_var, mods = ~CO2.Treatment,random=list(~factor(CO2.Treatment_level)|Paper_no,~Unit|Paper_no), data=data_wide_invasion, struct=c("CS", "CS"))
summary(mod.overall.best.food.growth.CO2.Treatment)
### need to have this CS... 


length(data_wide_invasion$LRR)
