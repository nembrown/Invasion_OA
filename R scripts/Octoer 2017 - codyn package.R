


####turnover
head(invasion.exp.data.not0)

data_long_invasion_turnover <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_turnover)
length(data_long_invasion_turnover$cover)

View(data_long_invasion_turnover)

head(invasion.exp.data.16)

install.packages("codyn")
library(codyn)
install.packages("knitr")
library(knitr)
data(knz_001d)
kable(head(knz_001d))


data_long_invasion_turnover_comm<- turnover(data_long_invasion_turnover, 
                                            time.var = "Week",
                                            species.var = "species",
                                            abundance.var = "cover", 
                                            replicate.var = "Tile.ID")



length(data_long_invasion_turnover_comm$total)

View(data_long_invasion_turnover_comm)
data_long_invasion_turnover_comm<-data_long_invasion_turnover_comm[order(data_long_invasion_turnover_comm$Week),]

View(invasion.exp.data.not0_higher2)
invasion.exp.data.not0_higher2<-invasion.exp.data.not0[invasion.exp.data.not0$Week>4, ]
invasion.exp.data.not0_higher2<-invasion.exp.data.not0_higher2[order(invasion.exp.data.not0_higher2$Week),]


head(invasion.exp.data.not0_higher2)
data_long_invasion_turnover_comm_bind<-as.data.frame(cbind(data_long_invasion_turnover_comm,invasion.exp.data.not0_higher2[, 1:10] ))
head(data_long_invasion_turnover_comm_bind)



cdata.turnover <- summaryBy(total ~ Treatment + Combined.Treatment + Invasives + CO2.Treatment + Week, data=data_long_invasion_turnover_comm_bind, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.turnover)[names(cdata.turnover)=="total.length"]<-"N"
cdata.turnover$total.se<-cdata.turnover$total.sd/sqrt(cdata.turnover$N)
head(cdata.turnover)

plot.turnover<-ggplot(cdata.turnover, aes(x=Week, y=total.mean, color = Invasives, shape=CO2.Treatment, linetype=CO2.Treatment)) + geom_errorbar(aes(ymin=total.mean-total.se, ymax=total.mean+total.se), width=0.25) + geom_line(size=1) + geom_point(size=5,aes(colour=factor(Invasives))) + scale_colour_manual(values=cbbPalette.all)+scale_shape_manual(values=c(19,2))
plot.turnover<- plot.turnover + labs(x = "Time (weeks)", y="Total species turnover")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.turnover<-plot.turnover + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.turnover <- plot.turnover + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.turnover 

plot.total<-ggplot(cdata.turnover, aes(x=Week, y=total.mean, shape=Treatment)) + geom_errorbar(aes(ymin=total.mean-total.se, ymax=total.mean+total.se), width=0.25) + geom_line(aes(linetype=CO2.Treatment), size=1) + geom_point(aes(shape=Treatment),size=5, col="#666666") +scale_shape_manual(values=c(1,19,2,17), name="Treatment: CO2, Invasives")+scale_linetype_manual(values=c("solid", "dashed"))+guides(lty=FALSE) 
plot.total<- plot.total + labs(x = "Time (weeks)", y="% cover Botryllus")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+xlim(0,16)
plot.total<-plot.total + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=12), legend.position=c(.17,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.01, "cm"))
plot.total <- plot.total+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.total 

lmer.turnover<-lmer(formula = total~ min.10.pH*Invasives*Week + (1|Mesocosm) + (0+Week|Mesocosm), data = data_long_invasion_turnover_comm_bind)
summary(lmer.turnover)
plot(lmer.turnover)
car::Anova(lmer.turnover, type="III")


##########appearance
data_long_invasion_appearance <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_appearance)

data_long_invasion_appearance_comm<- turnover(data_long_invasion_appearance, 
                                              time.var = "Week",
                                              species.var = "species",
                                              abundance.var = "cover", 
                                              replicate.var = "Tile.ID",
                                              metric="appearance")



kable(head(data_long_invasion_appearance_comm))

head(data_long_invasion_appearance_comm)
data_long_invasion_appearance_comm<-data_long_invasion_appearance_comm[order(data_long_invasion_appearance_comm$Week),]


invasion.exp.data.not0_higher2<-invasion.exp.data.not0[invasion.exp.data.not0$Week>4, ]
invasion.exp.data.not0_higher2<-invasion.exp.data.not0_higher2[order(invasion.exp.data.not0_higher2$Week),]

head(invasion.exp.data.not0_higher2)
data_long_invasion_appearance_comm_bind<-as.data.frame(cbind(data_long_invasion_appearance_comm,invasion.exp.data.not0_higher2[, 2:7] ))
head(data_long_invasion_appearance_comm_bind)


cdata.appearance <- summaryBy(appearance ~ Treatment + Combined.Treatment + Invasives + CO2.Treatment + Week, data=data_long_invasion_appearance_comm_bind, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.appearance)[names(cdata.appearance)=="appearance.length"]<-"N"
cdata.appearance$appearance.se<-cdata.appearance$appearance.sd/sqrt(cdata.appearance$N)
head(cdata.appearance)

plot.appearance<-ggplot(cdata.appearance, aes(x=Week, y=appearance.mean, shape=Treatment)) + geom_errorbar(aes(ymin=appearance.mean-appearance.se, ymax=appearance.mean+appearance.se), width=0.25) + geom_line(aes(linetype=CO2.Treatment), size=1) + geom_point(aes(shape=Treatment),size=5, col="#666666") +scale_shape_manual(values=c(1,19,2,17), name="Treatment: CO2, Invasives")+scale_linetype_manual(values=c("solid", "dashed"))+guides(lty=FALSE) 
plot.appearance<- plot.appearance + labs(x = "Time (weeks)", y="% cover Botryllus")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+xlim(0,16)
plot.appearance<-plot.appearance + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=12), legend.position=c(.17,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.01, "cm"))
plot.appearance <- plot.appearance+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.appearance

lmer.appearance<-lmer(formula = appearance~ min.10.pH*Invasives*Week + (1|Mesocosm) + (0+Week|Mesocosm), data = data_long_invasion_appearance_comm_bind)
summary(lmer.appearance)
plot(lmer.appearance)
car::Anova(lmer.appearance, type="III")


##########disappearance
data_long_invasion_disappearance <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_disappearance)

data_long_invasion_disappearance_comm<- turnover(data_long_invasion_disappearance, 
                                                 time.var = "Week",
                                                 species.var = "species",
                                                 abundance.var = "cover", 
                                                 replicate.var = "Tile.ID",
                                                 metric="disappearance")



kable(head(data_long_invasion_disappearance_comm))

head(data_long_invasion_disappearance_comm)
data_long_invasion_disappearance_comm<-data_long_invasion_disappearance_comm[order(data_long_invasion_disappearance_comm$Week),]


invasion.exp.data.not0_higher2<-invasion.exp.data.not0[invasion.exp.data.not0$Week>2, ]
invasion.exp.data.not0_higher2<-invasion.exp.data.not0_higher2[order(invasion.exp.data.not0_higher2$Week),]

head(invasion.exp.data.not0_higher2)
data_long_invasion_disappearance_comm_bind<-as.data.frame(cbind(data_long_invasion_disappearance_comm,invasion.exp.data.not0_higher2[, 2:7] ))
head(data_long_invasion_disappearance_comm_bind)


cdata.disappearance <- summaryBy(disappearance ~ Treatment + Combined.Treatment + Invasives + CO2.Treatment + Week, data=data_long_invasion_disappearance_comm_bind, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.disappearance)[names(cdata.disappearance)=="disappearance.length"]<-"N"
cdata.disappearance$disappearance.se<-cdata.disappearance$disappearance.sd/sqrt(cdata.disappearance$N)
head(cdata.disappearance)

plot.disappearance<-ggplot(cdata.disappearance, aes(x=Week, y=disappearance.mean, shape=Treatment)) + geom_errorbar(aes(ymin=disappearance.mean-disappearance.se, ymax=disappearance.mean+disappearance.se), width=0.25) + geom_line(aes(linetype=CO2.Treatment), size=1) + geom_point(aes(shape=Treatment),size=5, col="#666666") +scale_shape_manual(values=c(1,19,2,17), name="Treatment: CO2, Invasives")+scale_linetype_manual(values=c("solid", "dashed"))+guides(lty=FALSE) 
plot.disappearance<- plot.disappearance + labs(x = "Time (weeks)", y="% cover Botryllus")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+xlim(0,16)
plot.disappearance<-plot.disappearance + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=12), legend.position=c(.17,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.01, "cm"))
plot.disappearance <- plot.disappearance+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.disappearance


lmer.disappearance<-lmer(formula = disappearance~ min.10.pH*Invasives*Week + (1|Mesocosm) + (0+Week|Mesocosm), data = data_long_invasion_disappearance_comm_bind)
summary(lmer.disappearance)
plot(lmer.disappearance)
car::Anova(lmer.disappearance, type="III")

head(invasion.exp.data.not0)
##########rank_shift
data_long_invasion_rank_shift <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_rank_shift)

data_long_invasion_rank_shift_comm<- rank_shift(data_long_invasion_rank_shift, 
                                                time.var = "Week",
                                                species.var = "species",
                                                abundance.var = "cover", 
                                                replicate.var = "Tile.ID"
)



kable(head(data_long_invasion_rank_shift_comm))

head(data_long_invasion_rank_shift_comm)


invasion.exp.data.not0_higher2<-invasion.exp.data.not0[invasion.exp.data.not0$Week>2, ]
invasion.exp.data.not0_higher3<-invasion.exp.data.not0_higher2[order(invasion.exp.data.not0_higher2$Tile.ID),]

head(invasion.exp.data.not0_higher3)
data_long_invasion_rank_shift_comm_bind<-as.data.frame(cbind(data_long_invasion_rank_shift_comm,invasion.exp.data.not0_higher3[, 1:11] ))
head(data_long_invasion_rank_shift_comm_bind)


cdata.rank_shift <- summaryBy(MRS ~ Treatment + Combined.Treatment + Invasives + CO2.Treatment + Week, data=data_long_invasion_rank_shift_comm_bind, FUN=function(x) { c(mean = mean(x), sd = sd(x), length=length(x)) } )
names(cdata.rank_shift)[names(cdata.rank_shift)=="MRS.length"]<-"N"
cdata.rank_shift$MRS.se<-cdata.rank_shift$MRS.sd/sqrt(cdata.rank_shift$N)
head(cdata.rank_shift)

plot.rank_shift<-ggplot(cdata.rank_shift, aes(x=Week, y=MRS.mean, shape=Treatment)) + geom_errorbar(aes(ymin=MRS.mean-MRS.se, ymax=MRS.mean+MRS.se), width=0.25) + geom_line(aes(linetype=CO2.Treatment), size=1) + geom_point(aes(shape=Treatment),size=5, col="#666666") +scale_shape_manual(values=c(1,19,2,17), name="Treatment: CO2, Invasives")+scale_linetype_manual(values=c("solid", "dashed"))+guides(lty=FALSE) 
plot.rank_shift<- plot.rank_shift + labs(x = "Time (weeks)", y="rank shift")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+xlim(0,16)
plot.rank_shift<-plot.rank_shift + theme(legend.text = element_text(colour="black", size = 10))+theme(legend.title = element_text(colour="black", size=12), legend.position=c(.80,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.01, "cm"))
plot.rank_shift <- plot.rank_shift+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.rank_shift

#### model.invasion this one ... 
#### maybe use a glmm??admb? 

lmer.rank_shift<-lmer(formula = MRS ~ min.10.pH*Invasives*Week + (1|Mesocosm) + (0+Week|Mesocosm), data = data_long_invasion_rank_shift_comm_bind)
summary(lmer.rank_shift)
plot(lmer.rank_shift)
car::Anova(lmer.rank_shift, type="III")

anova(lmer.tile.biomass, ddf = "Kenward-Roger")



##### rate change
data_long_invasion_rate_change <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_rate_change)
View(data_long_invasion_rate_change)

(invasion.exp.data.not0[,38])


data_long_invasion_rate_change_comm<- rate_change(data_long_invasion_rate_change, 
                                                  time.var = "Week",
                                                  species.var = "species",
                                                  abundance.var = "cover", 
                                                  replicate.var = "Tile.ID"
)



kable(head(data_long_invasion_rate_change_comm))

head(data_long_invasion_rate_change_comm)

head(invasion.exp.data.16)
data_long_invasion_rate_change_comm_bind<-as.data.frame(cbind(data_long_invasion_rate_change_comm,invasion.exp.data.16[, 1:11] ))
head(data_long_invasion_rate_change_comm_bind)


plot.rate_change<- ggplot(data_long_invasion_rate_change_comm_bind, aes(x=Invasives, y=rate_change, colour=Treatment, fill=Treatment))+ scale_colour_manual(values=cbbPalette.boxplot) + geom_boxplot(position = "dodge", alpha=0.70)  + scale_fill_manual(values=cbbPalette.white, name="Treatment")
plot.rate_change<- plot.rate_change + theme_bw() + xlab(bquote('Food supply level')) + ylab(expression("Slope of rate change"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.rate_change<- plot.rate_change + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.rate_change<- plot.rate_change+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.rate_change

plot.rate_change.16.invasion<- ggplot(data_long_invasion_rate_change_comm_bind, aes(x=invasion.exp.data.16$hydrogen.concentration, y=(data_long_invasion_rate_change_comm_bind$rate_change), col="#666666")) + geom_point(size=5,aes(shape=Treatment), col="#666666") + guides(fill=FALSE) + geom_smooth(aes(lty=Invasives, fill=Invasives), colour="#666666", alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(1,19,2,17))+scale_linetype_manual(values=c("dashed", "solid"))+ scale_fill_manual(values=c( "#b2b2b2","#666666"))
plot.rate_change.16.invasion<- plot.rate_change.16.invasion + theme_bw() + xlab(bquote('Minimum 10th percentile pH')) + ylab("Rate change")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.rate_change.16.invasion<- plot.rate_change.16.invasion + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.rate_change.16.invasion<- plot.rate_change.16.invasion+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.rate_change.16.invasion <-plot.rate_change.16.invasion +  scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7))
plot.rate_change.16.invasion 


qqp(data_long_invasion_rate_change_comm_bind$rate_change, "norm")

lm.gamma.rate_change.12<- lm(formula = rate_change ~ min.10.pH*Invasives, data = data_long_invasion_rate_change_comm_bind)
plot(lm.gamma.rate_change.12)

summary(lm.gamma.rate_change.12)
car::Anova(lm.gamma.rate_change.12, type=2)

cbbPalette.boxplot<- c( "#F8766D", "#F8766D","#00BA38",  "#00BA38", "#619CFF", "#619CFF")
cbbPalette.all<- c( "#F8766D", "#00BA38",  "#619CFF", "#F8766D", "#00BA38","#619CFF")


cbbPalette.white<- c( "#F8766D", "#FFFFFF","#00BA38", "#FFFFFF","#619CFF",   "#FFFFFF")

#####stability
data_long_invasion_stability <- gather(invasion.exp.data.not0, species, cover, c(12,14,22:25,28, 32:38,43:47), factor_key=TRUE)
head(data_long_invasion_stability)


data_long_invasion_stability_comm<- community_stability(data_long_invasion_stability, 
                                                        time.var = "Week", 
                                                        abundance.var = "cover", 
                                                        replicate.var = "Tile.ID")

kable(head(data_long_invasion_stability_comm))

head(data_long_invasion_stability_comm)


data_long_invasion_stability_comm_bind<-as.data.frame(cbind(data_long_invasion_stability_comm,invasion.exp.data.16[, 1:11] ))
head(data_long_invasion_stability_comm_bind)

plot.stability.16.invasion<- ggplot(data_long_invasion_stability_comm_bind, aes(x=invasion.exp.data.16$hydrogen.concentration, y=(data_long_invasion_stability_comm_bind$stability), col="#666666")) + geom_point(size=5,aes(shape=Treatment), col="#666666") + guides(fill=FALSE) + geom_smooth(aes(lty=Invasives, fill=Invasives), colour="#666666", alpha=0.15,size=1.5,  method="glm", method.args = list(family = "Gamma")) +scale_shape_manual(values=c(1,19,2,17))+scale_linetype_manual(values=c("dashed", "solid"))+ scale_fill_manual(values=c( "#b2b2b2","#666666"))
plot.stability.16.invasion<- plot.stability.16.invasion + theme_bw() + xlab(bquote('Minimum 10th percentile pH')) + ylab("Stability")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.stability.16.invasion<- plot.stability.16.invasion + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability.16.invasion<- plot.stability.16.invasion+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability.16.invasion <-plot.stability.16.invasion +  scale_x_reverse( breaks=c(10^-7.2, 10^-7.3,10^-7.4,10^-7.5,10^-7.6, 10^-7.7), labels=c( 7.2,7.3, 7.4,7.5, 7.6, 7.7))
plot.stability.16.invasion 

plot.stability.by.bot<- ggplot(data_long_invasion_stability_comm_bind, aes(x=invasion.exp.data.10$botryllid, y=(stability), colour=Treatment)) + geom_point(size=5,aes(shape=Treatment)) + guides(fill=FALSE) + geom_smooth(aes(lty=CO2.Treatment, fill=CO2.Treatment,col=CO2.Treatment), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(1,19,2,17))+scale_linetype_manual(values=c(1,6 ))+ scale_fill_manual(values=c( "#666666","#619CFF"))+guides( lty=FALSE, col=FALSE)+scale_color_manual(values=cbbPalette.blue6)
plot.stability.by.bot<- plot.stability.by.bot + theme_bw() + xlab(bquote('% cover Botryllus')) + ylab("Proportion cover stability")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.stability.by.bot<- plot.stability.by.bot + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16), legend.position=c(.80,.80))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability.by.bot<- plot.stability.by.bot+ theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability.by.bot 


gamma.12<-fitdistr(data_long_invasion_stability_comm_bind$stability, "gamma")
qqp(data_long_invasion_stability_comm_bind$stability, "gamma", shape = gamma.12$estimate[[1]], rate = gamma.12$estimate[[2]])

qqp(data_long_invasion_stability_comm_bind$stability, "norm")

glm.gamma.stability.12<- glm(formula = stability ~ CO2.Treatment*Invasives, data = data_long_invasion_stability_comm_bind, family = "Gamma")
plot(glm.gamma.stability.12)

summary(glm.gamma.stability.12)
car::Anova(glm.gamma.stability.12, type=2)

library(car)

lm.stability<-lm(formula = stability ~ min.10.pH*Invasives, data = data_long_invasion_stability_comm_bind)
summary(lm.stability)
plot(lm.stability)
car::Anova(lm.stability, type=3)

##### Stability (using not just cover)
#####stability
head(invasion.exp.data.not0)
View(invasion.exp.data.not0)

data_long_invasion_stability_rec_cover <- gather(invasion.exp.data.not0, species, cover, c(8, 9, 11, 12, 21, 22, 25, 28,29, 31, 38, 45, 49,50, 56, 57, 60, 64, 66:71), factor_key=TRUE)
head(data_long_invasion_stability_rec_cover)
View(data_long_invasion_stability_rec_cover)
library(ggplot2)
View(invasion.exp.data.16)
View(data_long_invasion_stability_rec_cover)
install.packages("codyn")
library(codyn)
install.packages("knitr")
library(knitr)
data(knz_001d)
kable(head(knz_001d))


data_long_invasion_stability_rec_cover_comm<- community_stability(data_long_invasion_stability_rec_cover, 
                                                                  time.var = "Week", 
                                                                  abundance.var = "cover", 
                                                                  replicate.var = "Mesocosm")

kable(head(data_long_invasion_stability_rec_cover_comm))

head(data_long_invasion_stability_rec_cover_comm)
View(data_long_invasion_stability_rec_cover_comm)


data_long_invasion_stability_rec_cover_comm_bind<-as.data.frame(cbind(data_long_invasion_stability_rec_cover_comm[-46,],invasion.exp.data.16[, 2:7] ))
head(data_long_invasion_stability_rec_cover_comm_bind)
View(data_long_invasion_stability_rec_cover_comm_bind)

plot.stability_rec_cover<- ggplot(data_long_invasion_stability_rec_cover_comm_bind, aes(x=min.10.pH, y=stability, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives))) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Invasives, fill=Invasives),method="glm",  method.args = list(family = "Gamma"),alpha=0.15,size=1.5)
plot.stability_rec_cover<- plot.stability_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability_rec_cover<- plot.stability_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability_rec_cover<- plot.stability_rec_cover+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability_rec_cover

plot.stability_rec_cover<- ggplot(data_long_invasion_stability_rec_cover_comm_bind, aes(x=min.10.pH, y=stability, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives))) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Invasives, fill=Invasives),method="lm",alpha=0.15,size=1.5)
plot.stability_rec_cover<- plot.stability_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability_rec_cover<- plot.stability_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability_rec_cover<- plot.stability_rec_cover+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability_rec_cover

plot.stability_rec_cover<- ggplot(data_long_invasion_stability_rec_cover_comm_bind, aes(x=Combined.Treatment, y=stability, colour=Invasives))+ scale_colour_manual(values=cbbPalette.all) + geom_boxplot(position = "dodge") + guides(fill=CO2.Treatment) 
plot.stability_rec_cover<- plot.stability_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability_rec_cover<- plot.stability_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability_rec_cover<- plot.stability_rec_cover+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability_rec_cover



glm.gamma.stability_rec_cover.12<- glm(formula = stability ~ min.10.pH*Invasives, data = data_long_invasion_stability_rec_cover_comm_bind, family = "Gamma")
plot(glm.gamma.stability_rec_cover.12)

summary(glm.gamma.stability_rec_cover.12)
library(car)
Anova(glm.gamma.stability_rec_cover.12, type=3)


lm.stability_rec_cover<-lm(formula = stability ~ min.10.pH*Invasives, data = data_long_invasion_stability_rec_cover_comm_bind)
summary(lm.stability_rec_cover)
plot(lm.stability_rec_cover)
Anova(lm.stability_rec_cover, type=3)






##### stability.sessile with sessile
head(invasion.exp.data.not0)
data_long_invasion_stability.sessile_rec_cover <- gather(invasion.exp.data.not0, species, cover, c(8, 9, 11, 12,  28,29, 31, 38, 45, 49,50, 56, 57, 60, 64, 71), factor_key=TRUE)

View(data_long_invasion_stability.sessile_rec_cover)


View(invasion.exp.data.16)
View(data_long_invasion_stability.sessile_rec_cover)
install.packages("codyn")
library(codyn)
install.packages("knitr")
library(knitr)
data(knz_001d)
kable(head(knz_001d))


data_long_invasion_stability.sessile_rec_cover_comm<- community_stability(data_long_invasion_stability.sessile_rec_cover, 
                                                                          time.var = "Week", 
                                                                          abundance.var = "cover", 
                                                                          replicate.var = "Mesocosm")

kable(head(data_long_invasion_stability.sessile_rec_cover_comm))

head(data_long_invasion_stability.sessile_rec_cover_comm)
View(data_long_invasion_stability.sessile_rec_cover_comm)


data_long_invasion_stability.sessile_rec_cover_comm_bind<-as.data.frame(cbind(data_long_invasion_stability.sessile_rec_cover_comm[-46,],invasion.exp.data.16[, 2:7] ))
head(data_long_invasion_stability.sessile_rec_cover_comm_bind)
View(data_long_invasion_stability.sessile_rec_cover_comm_bind)

plot.stability.sessile_rec_cover<- ggplot(data_long_invasion_stability.sessile_rec_cover_comm_bind, aes(x=min.10.pH, y=stability, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives))) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Invasives, fill=Invasives),method="glm",  method.args = list(family = "Gamma"),alpha=0.15,size=1.5)
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability.sessile_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability.sessile_rec_cover

plot.stability.sessile_rec_cover<- ggplot(data_long_invasion_stability.sessile_rec_cover_comm_bind, aes(x=min.10.pH, y=stability, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives))) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Invasives, fill=Invasives),method="lm",alpha=0.15,size=1.5)
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability.sessile_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability.sessile_rec_cover

plot.stability.sessile_rec_cover<- ggplot(data_long_invasion_stability.sessile_rec_cover_comm_bind, aes(x=Combined.Treatment, y=stability, colour=Invasives))+ scale_colour_manual(values=cbbPalette.all) + geom_boxplot(position = "dodge") + guides(fill=CO2.Treatment) 
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme_bw() + xlab(bquote('min.10.pH')) + ylab(expression("Species stability.sessile_rec_cover"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.stability.sessile_rec_cover<- plot.stability.sessile_rec_cover+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.stability.sessile_rec_cover



gamma.12<-fitdistr(data_long_invasion_stability.sessile_rec_cover_comm_bind$stability + 0.01, "gamma")
qqp(data_long_invasion_stability.sessile_rec_cover_comm_bind$stability, "gamma", shape = gamma.12$estimate[[1]], rate = gamma.12$estimate[[2]])

glm.gamma.stability.sessile_rec_cover.12<- glm(formula = stability ~ min.10.pH*Invasives, data = data_long_invasion_stability.sessile_rec_cover_comm_bind, family = "Gamma")
plot(glm.gamma.stability.sessile_rec_cover.12)

summary(glm.gamma.stability.sessile_rec_cover.12)
Anova(glm.gamma.stability.sessile_rec_cover.12, type=3)


lm.stability.sessile_rec_cover<-lm(formula = stability ~ CO2.Treatment*Invasives, data = data_long_invasion_stability.sessile_rec_cover_comm_bind)
summary(lm.stability.sessile_rec_cover)
plot(lm.stability.sessile_rec_cover)
Anova(lm.stability.sessile_rec_cover, type=3)

cbbPalette.all<- c( "#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")
