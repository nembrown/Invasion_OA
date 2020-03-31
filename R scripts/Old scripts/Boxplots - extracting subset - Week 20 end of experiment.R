mydataweek20<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
head(mydataweek20)

fix(mydataweek20)

mydataweek20.Meso<-mydataweek20[mydataweek20$Description=="Mesocosm", ]

head(mydataweek20.Meso)
   
mydataweek20.Header<-mydataweek20[mydataweek20$Description=="Header",]
head(mydataweek20.Header)

#PLOTS

p<- ggplot(mydataweek20.Meso, aes(x=mydataweek20.Meso$Treatment, y=mydataweek20.Meso$serpulid, fill=mydataweek20.Meso$Treatment)) + geom_boxplot() + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="Treatment", y="# serpulid\nper mesocosm")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())



# Overlaid histograms
#ggplot(df, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")

# Interleaved histograms
#ggplot(df, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5, position="dodge")

#HIstorgram
p1<- ggplot(mydataweek20.Meso, aes(x=mydataweek20.Meso$serpulid, fill=mydataweek20.Meso$Treatment)) + geom_histogram(binwidth=1, alpha=.5, position="dodge") + scale_fill_brewer(palette="Set1")

p1<-p1 + theme_bw() + labs(x="serpulid", y="Count")+ theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p1 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

#Density plot
p4<- ggplot(mydataweek20.Meso, aes(x=mydataweek20.Meso$serpulid, fill=mydataweek20.Meso$Treatment)) + geom_density(alpha=.3) + scale_fill_brewer(palette="Set1")

p4<-p4 + theme_bw() + labs(x="serpulid", y="Density")+ theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p4 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())



### by pH
p<- ggplot(mydataweek20.Meso, aes(x=mydataweek20.Meso$pH, y=mydataweek20.Meso$serpulid, colour=mydataweek20.Meso$Treatment)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1")

p<- p + theme_bw() + labs(x="pH", y="# serpulid colonies")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

#WIth lm model
p<- ggplot(mydataweek20.Meso, aes(x=mydataweek20.Meso$pH, y=mydataweek20.Meso$serpulid)) + geom_point( size=5) + guides(fill=FALSE) + scale_fill_brewer(palette="Set1") + geom_smooth(method=lm)

p<- p + theme_bw() + labs(x="pH", y="# serpulid\nper Meso")  + theme(text = element_text(size=16), axis.text = element_text(size=14))+theme(axis.title.y = element_text(angle=0))

p + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())


#STATS




#Should try quantile regression with these

library(quantreg)
head(mydataweek20.Meso)


Rm.rq<-rq(serpulid ~ pH, data=mydataweek20.Meso, tau=0.5)
r1.pH<-resid(Rm.rq)
plot(Rm.rq)

hist(r1.pH)
c1.pH<-coef(Rm.rq)

summary(Rm.rq)
summary(Rm.rq, se="nid")
plot(x=mydataweek20.Meso$pH, y=mydataweek20.Meso$serpulid, cex=.25, type="n", xlab="pH", ylab="serpulid per meso")
points(x=mydataweek20.Meso$pH, y=mydataweek20.Meso$serpulid, cex=.75, col="blue")
abline(rq(serpulid ~ pH,data=mydataweek20.Meso, tau=.5),col="blue")
abline(lm(serpulid ~ pH,data=mydataweek20.Meso),lty=2,col="red")
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){abline(rq(serpulid~pH,data=mydataweek20.Meso,tau=taus[i]),col="gray")}


summary((rq(serpulid ~ pH, data=mydataweek20.Meso, tau=0.80)), se="nid")


###Lm with pH
z2.lm<-lm(serpulid ~ pH, data=mydataweek20.Meso)
plot(z2.lm)
summary(z2.lm)

z3.glm<-glm(serpulid ~ pH, family=poisson, data=mydataweek20.Meso)
plot(z3.glm)
summary(z3.glm)


############# Need to really to do poisson distribution (or any of the distributions) with pH and treatment - see if it fits the distribution
# Because hard to justify - the ones iwth a lot of zeros then just apply a lot of weight to that. 
# but maybe I need to be doing summary of 75 quantile?? Doesn't seem to work very well (problems with tau + h > 1... ?)

### If I'm going to use Glm need to decide if that fits ... and if going to use qunatile need to decide if that fits
#The best. 



hist(mydataweek20.Meso$serpulid)

z.lm<-lm(serpulid ~ Treatment, family=poisson, data=mydataweek20.Meso)
plot(z.lm)
summary(z.lm)

z.glm<-glm(serpulid ~ Treatment, family=poisson, data=mydataweek20.Meso)
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

z.glm.quasi<-glm(serpulid ~ Treatment, family=quasipoisson, data=mydataweek20.Meso)
plot(z.glm.quasi)
summary(z.glm.quasi)
#doesn't look like a great model because the estimates are not correct


mydataweek20.Meso$Treatment<-as.factor(mydataweek20.Meso$Treatment)

kruskal.test(mydataweek20.Meso$serpulid~mydataweek20.Meso$Treatment)
head(mydataweek20.Meso)
