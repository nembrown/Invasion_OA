repeatedm.stats<-read.csv(file.choose())
#File: All week R input .csv

levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="AIR, Absent"] <- "Control, Invasives Absent"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="CO2, Absent"] <- "Elevated CO2, Invasives Absent"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="AIR, Present"] <- "Control, Invasives Present"
levels(repeatedm.stats$Treatment)[levels(repeatedm.stats$Treatment)=="CO2, Present"] <- "Elevated CO2, Invasives Present"
levels(repeatedm.stats$CO2.Treatment)[levels(repeatedm.stats$CO2.Treatment)=="AIR"] <- "Control"
levels(repeatedm.stats$CO2.Treatment)[levels(repeatedm.stats$CO2.Treatment)=="CO2"] <- "Elevated CO2"


head(repeatedm.stats)
repeatedm.stats<-repeatedm.stats[which(repeatedm.stats$Week>0), ]
repeatedm.stats$Mesocosm<-as.factor(repeatedm.stats$Mesocosm)
repeatedm.stats$Week<-as.factor(repeatedm.stats$Week)

repeatedm.stats16<-repeatedm.stats[repeatedm.stats$Week==16, ]

repeatedm.stats14<-repeatedm.stats[repeatedm.stats$Week==14, ]

repeatedm.stats12<-repeatedm.stats[repeatedm.stats$Week==12, ]

repeatedm.stats10<-repeatedm.stats[repeatedm.stats$Week==10, ]

repeatedm.stats8<-repeatedm.stats[repeatedm.stats$Week==8, ]

repeatedm.stats6<-repeatedm.stats[repeatedm.stats$Week==6, ]

repeatedm.stats4<-repeatedm.stats[repeatedm.stats$Week==4, ]

repeatedm.stats2<-repeatedm.stats[repeatedm.stats$Week==2, ]


library(plyr)

library(ggplot2)

library(doBy)



####### LME 
library(nlme)
library(visreg)
library(lme4)

####################################################### 
###Just regular plotting with pH.uptowk

p<- ggplot(repeatedm.stats, aes(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2, colour=repeatedm.stats$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH.uptowk", y="prop.mem.eaten2 percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p



#WIth lm model
p<- ggplot(repeatedm.stats, aes(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH.uptowk", y="prop.mem.eaten2 percent cover")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p



###############################################################################################################################

### Quantile regression !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(quantreg)
library(lqmm)
#lqmm (y~x + covariates , random=~1,group=subjectID, iota=c(0.25, 0.50, 0.75),.)
# best model: prop.mem.eaten2.PQL2<- glmmPQL(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1|Mesocosm/Tile, family=poisson, data=repeatedm.stats)

num.species.lqmm<- lqmm(y~x + covariates , random=~1,group=subjectID, iota=c(0.25, 0.50, 0.75),.)



#####Yes there is quantile regression using lqmm for formulas with random effects but I'm finding difficult o visiualize residuals etc ... 
#So might not be the best ... maybe just try a glm ... try getting Week in there. 
### residuals not great -- maybe can try quantile regression?? lqmm
library(lqmm)
#lqmm(fixed, random, group, covariance = "pdDiag", tau = 0.5,
# nK = 7, type = "normal", rule = 1, data = sys.frame(sys.parent()),
#  subset, weights, na.action = na.fail, control = list(),
#   contrasts = NULL, fit = TRUE)
#what is nk=7?????number of quadrature knots. is it n+1? Read the Geraci paper 
# should I have week as a factor? Then I see it at every level?? -- if I do this does it show up at the same levels as the glmm? 

repeatedm.stats$Week<-as.factor(repeatedm.stats$Week)
# should it be random=~1|Mesocosm .... that doesn't work. #### NOT SURE ABOUT THE RANDOM PART OF THIS!!! AND GROUP??? 

zlqmm<-lqmm(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~Mesocosm, group= Tile.ID, tau = 0.5,
            nK = 7, type = "normal", rule = 1, data = repeatedm.stats, na.action = na.fail, control = list(),
            contrasts = NULL, fit = TRUE)
summary(zlqmm)

plot(residuals(zlqmm), fitted(zlqmm)) ## actually look okay ... but doesn't need this b/c not a mean. 



### can I make this work for lqmm?? ---> difficult ... 

plot(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2, cex=.25, type="n", xlab="pH", ylab="prop.mem.eaten2")
points(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2, cex=.75, col="blue")
abline(lqmm(prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1, group= Tile.ID, tau = 0.5,
            nK = 7, type = "normal", rule = 1, data = repeatedm.stats, na.action = na.fail, control = list(),
            contrasts = NULL, fit = TRUE,col="blue")
       abline(lm(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats),lty=2,col="red")
       taus <- c(.05,.1,.25,.75,.90,.95)
       for( i in 1:length(taus)){abline(lqmm((prop.mem.eaten2 ~ CO2.Treatment*Invasives*Week, random=~1, group= Tile.ID, tau = 0.5,
                                              nK = 7, type = "normal", rule = 1, data = repeatedm.stats, na.action = na.fail, control = list(),
                                              contrasts = NULL, fit = TRUE, tau=taus[i]),col="gray")}
       
       
       
       
### Quantile regression without random factors ... not useful really unless don't have Invasives as a factor
       
head(repeatedm.stats)



zbot<-lm(repeatedm.stats$prop.mem.eaten2 ~ repeatedm.stats$pH.uptowk)
       plot(zbot)
       summary(zbot)
       
       library(quantreg)
       
       zbot.rq<-rq(prop.mem.eaten2 ~ pH.uptowk, data=repeatedm.stats, tau=0.5)
       r1.pH.uptowk<-resid(zbot.rq)
       
       plot(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2, cex=.25, type="n", xlab="pH.uptowk", ylab="prop.mem.eaten2 of prop.mem.eaten2s")
       points(x=repeatedm.stats$pH.uptowk, y=repeatedm.stats$prop.mem.eaten2, cex=.75, col="blue")
       abline(rq(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats, tau=.5),col="blue")
       abline(lm(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats),lty=2,col="red")
       taus <- c(.05,.1,.25,.75,.90,.95)
       for( i in 1:length(taus)){abline(rq(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats,tau=taus[i]),col="gray")}
       
       
       summary(zbot.rq)
       summary(zbot.rq, se="nid")
       
       
       ### So quantile regression will probably only be feasible if I DON'T have invasives in there ... i.e. in things that involve botryllids
       #Becaus ethen I can discount the botryllids and then I can do just by CO2 or JUST by pHuptowk. 
       #Also will only be feasible if I do ONE WEEK 
       #Theremore might be useful for botryllid percent cover, bot. eaten at certain weeks. 
       
#### By week::: 

head(repeatedm.stats14)


library(quantreg)
       head(repeatedm.stats14)
       repeatedm.stats14.bot.only<-repeatedm.stats14[repeatedm.stats14$Invasives=="Present",]
       
       
       Rm.rq<-rq(prop.mem.eaten2 ~ pH.uptowk, data=repeatedm.stats14.bot.only, tau=0.5)
       r1.pH<-resid(Rm.rq)
       plot(Rm.rq)
       
       hist(r1.pH)
       c1.pH<-coef(Rm.rq)
       summary(Rm.rq)
       summary(Rm.rq, se="nid")
       plot(x=repeatedm.stats14.bot.only$pH.uptowk, y=repeatedm.stats14.bot.only$prop.mem.eaten2, cex=.25, type="n", xlab="pH", ylab="bot eaten")
       points(x=repeatedm.stats14.bot.only$pH.uptowk, y=repeatedm.stats14.bot.only$prop.mem.eaten2, cex=.75, col="blue")
       abline(rq(prop.mem.eaten2 ~ pH.uptowk,data=repeatedm.stats14.bot.only, tau=.5),col="blue")
       abline(lm(prop.mem.eaten2 ~ pH.uptowk,data=repeatedm.stats14.bot.only),lty=2,col="red")
       taus <- c(.05,.1,.25,.75,.90,.95)
       for( i in 1:length(taus)){abline(rq(prop.mem.eaten2~pH.uptowk,data=repeatedm.stats14.bot.only,tau=taus[i]),col="gray")}
       
       