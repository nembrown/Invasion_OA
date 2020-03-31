mydataweek20<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
head(mydataweek20)

fix(mydataweek20)

mydataweek20<-mydataweek20[mydataweek20$Description=="Mesocosm", ]

head(mydataweek20)
   
mydataweek20<-mydataweek20[mydataweek20$Description=="Meso",]
head(mydataweek20)

#PLOTS

p<- ggplot(mydataweek20, aes(x=mydataweek20$pH, y=mydataweek20$sum.serpulid, fill=mydataweek20$pH)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH", y="# sum.serpulid\nper mesocosm")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())



# Overlaid histograms
#ggplot(df, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5, alpHa=.5, position="identity")

# Interleaved histograms
#ggplot(df, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5, position="dodge")

#HIstorgram
p1<- ggplot(mydataweek20, aes(x=mydataweek20$sum.serpulid, fill=mydataweek20$pH)) + geom_histogram(binwidth=1, alpha=.5, position="dodge") + scale_fill_brewer(palette="Set1")

p1<-p1 + theme_bw() + labs(x="sum.serpulid", y="Count")+ theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p1 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

#Density plot
p4<- ggplot(mydataweek20, aes(x=mydataweek20$sum.serpulid, fill=mydataweek20$pH)) + geom_density(alpha=.3) + scale_fill_brewer(palette="Set1")

p4<-p4 + theme_bw() + labs(x="sum.serpulid", y="Density")+ theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p4 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())



### by pH
p<- ggplot(mydataweek20, aes(x=mydataweek20$pH, y=mydataweek20$sum.serpulid, colour=mydataweek20$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH", y="# sum.serpulid colonies")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

#WIth lm model
p<- ggplot(mydataweek20, aes(x=mydataweek20$pH, y=mydataweek20$sum.serpulid)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH", y="# sum.serpulid\nper Meso")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())





p<-ggplot(mydataweek20, aes(x=mydataweek20$pH, y=mydataweek20$sum.serpulid, colour=Treatment)) + geom_point(size=4) + scale_colour_brewer(palette="Set1")+ geom_smooth(method=lm, colour="black")
#rename my plot to p, makes it easier to work with

p<- p + labs(x="pH", y="# of serpulids\nper header\n& mesocosm system")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=14)) + theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())






#STATS
zlm<-lm(sum.serpulid ~ pH,data=mydataweek20)
plot(zlm)
summary(zlm)



#Should try quantile regression with these

library(quantreg)
head(mydataweek20)


Rm.rq<-rq(sum.serpulid ~ pH, data=mydataweek20, tau=0.5)
r1.pH<-resid(Rm.rq)
plot(Rm.rq)

hist(r1.pH)
c1.pH<-coef(Rm.rq)

summary(Rm.rq)
summary(Rm.rq, se="nid")
plot(x=mydataweek20$pH, y=mydataweek20$sum.serpulid, cex=.25, type="n", xlab="pH", ylab="sum.serpulid per Meso")
points(x=mydataweek20$pH, y=mydataweek20$sum.serpulid, cex=.75, col="blue")
abline(rq(sum.serpulid ~ pH,data=mydataweek20, tau=.5),col="blue")
abline(lm(sum.serpulid ~ pH,data=mydataweek20),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(sum.serpulid~pH,data=mydataweek20,tau=taus[i]),col="gray")}


summary((rq(sum.serpulid ~ pH, data=mydataweek20, tau=0.75)), se="nid")


###Lm with pH
z2.lm<-lm(sum.serpulid ~ pH, data=mydataweek20)
plot(z2.lm)
summary(z2.lm)

z3.glm<-glm(sum.serpulid ~ pH, family=poisson, data=mydataweek20)
plot(z3.glm)
summary(z3.glm)


############# Need to really to do poisson distribution (or any of the distributions) with pH and pH - see if it fits the distribution
# Because hard to justify - the ones iwth a lot of zeros then just apply a lot of weight to that. 
# but maybe I need to be doing summary of 75 quantile?? Doesn't seem to work very well (problems with tau + h > 1... ?)

### If I'm going to use Glm need to decide if that fits ... and if going to use qunatile need to decide if that fits
#The best. 

head(mydataweek20)

hist(mydataweek20$sum.serpulid)

z.lm<-lm(sum.serpulid ~ pH, family=poisson, data=mydataweek20)
plot(z.lm)
summary(z.lm)

z.glm<-glm(sum.serpulid ~ pH, family=poisson, data=mydataweek20)
plot(z.glm)
summary(z.glm)
#deviance:
#Deviance // model fit
#Explained deviance is 100 x null-resid/null
d=100*(316.13-299.71)/316.13
d
#5% deviance explained. 

### OVERDISPERSION
#Check for overdispersion o=D/(n-p) ... so the estimator o is the ration of Residual deviance over degrees of freedom. 

299.71/22
#overdispersed

library(AER)
dispersiontest(z.glm,trafo=1)

#If C=0 then not overdispersed, so looks like it is definitely overdispersed because c= 9.6
#try quasipoisson

z.glm.quasi<-glm(sum.serpulid ~ pH, family=quasipoisson, data=mydataweek20)
plot(z.glm.quasi)
summary(z.glm.quasi)
#doesn't look like a great model because the estimates are not correct


mydataweek20$pH<-as.factor(mydataweek20$pH)

kruskal.test(mydataweek20$sum.serpulid~mydataweek20$pH)
head(mydataweek20)
