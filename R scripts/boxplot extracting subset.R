mydata<-read.csv(file.choose(),stringsAsFactors = FALSE, na.strings = c("NA","") )
head(mydata)

fix(mydata)

mydata14<-mydata[mydata$Week==14, ]

head(mydata14)

mydata14present<-mydata14[mydata14$Invasives=="Present",]


#MUSSEL RECRUITMENT *
boxplot(mydata14present$botryllid~mydata14present$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="% botryllid")

zmuss<-lm(mydata14present$botryllid~mydata14$Treatment)
plot(zmuss)
anova(zmuss)
summary(zmuss)

#Residual vs fitted, make a funnel shape so heterogenous ... therefore need to use non-parametric. 
kruskal.test(mydata$botryllid~mydata$Treatment)
#Significant --> p=0.0432, df=1, H=4.205

#----------------------------------------------------------------------
#BRYOZOAN RECRUITMENT *
boxplot(mydata$av.bryo~mydata$CO2.Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of M. membranacea")

zbryo<-lm(mydata$av.bryo~mydata$Treatment)
plot(zbryo)
#Residuals vs fitted make a funnel shape, heterogenous ... therefore need to use non parametric
kruskal.test(mydata$av.bryo~mydata$Treatment)
#Significant --> p=0.03263, df=1, H=4.5653

#-----------------------------------------------------------------------
#BARNACLE RECRUITMENT
boxplot(mydata$av.barnacle~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of B. crenatus")

zbarnacle<-lm(mydata$av.barnacle~mydata$Treatment)
plot(zbarnacle)
#VERY not normal ... odd plot, residuals vs fitted not great
kruskal.test(mydata$av.barnacle~mydata$Treatment)
#NS --> p=0.6198, df=1, H=0.2462

#BARNACLE SIZE
boxplot(mydata$av.large.barn~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Basal area of largest B. crenatus (cm2)")
zbarnsize<-lm(mydata$av.large.barn~mydata$Treatment)
plot(zbarnsize)
#fairly normal
summary(zbarnsize)
anova(zbarnsize)
#NS --> p=0.8666, F=0.0289, df=1
#Get parameter estimates (means)
#------------------------------------------------------------------------
#BS Tunicate RECRUITMENT
boxplot(mydata$av.BS.tunicate~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of B. schlosseri")
stripchart(mydata$av.BS.tunicate~mydata$Treatment, xlab="Treatment", ylab="Percent cover of B. schlosseri", method="jitter", pch=16, vertical=TRUE)
zBschloss<-lm(mydata$av.BS.tunicate~mydata$Treatment)
plot(zBschloss)
#Fitted vs residuals make a funnel
kruskal.test(mydata$av.BS.tunicate~mydata$Treatment)
#NS ---> p=0.5147 H=0.4245 df=1

#BS SIZE
head(mydata)
boxplot(mydata$av.larg.BS.tun~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Area of largest B. schlosseri colony (cm2)")
zBSsize<-lm(mydata$av.larg.BS.tun~mydata$Treatment)
plot(zBSsize)
#NOT normal, unequal variances
kruskal.test(mydata$av.larg.BS.tun~mydata$Treatment)
#NS --> p=0.5693 H=0.3239, df=1

#------------------------------------------------------------------
#CI NATIVE RECRUITMENT
boxplot(mydata$av.CI.native.tun~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of C. inflata")
zCinflat<-lm(mydata$av.CI.native.tun~mydata$Treatment)
summary(zCinflat)
plot(zCinflat)
#No difference, all zeros

#CI SIZE
head(mydata)
boxplot(mydata$av.larg.CI.native~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Area of largest C. inflata (cm2)")
zsize.Cinflat<-lm(mydata$av.larg.CI.native~mydata$Treatment)
summary(zsize.Cinflat)
plot(zsize.Cinflat)
#residuals make a funnel.
kruskal.test(mydata$av.larg.CI.native~mydata$Treatment)
#NS --> p=0.1486, H=2.087, df=1

#---------------------------------------------------------------------
#HYDROID RECRUITMENT *
boxplot(mydata$av.hydroid~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of hydroid")
stripchart(mydata$av.hydroid~mydata$Treatment, xlab="Treatment", ylab="Percent cover of O. dichotoma", method="jitter", pch=16, vertical=TRUE)
stripchart(fitted(zhyd)~mydata$Treatment, vertical=TRUE, add=TRUE, pch="-----", method="jitter")

zhyd<-lm(mydata$av.hydroid~mydata$Treatment)
plot(zhyd)
#Residual vs fitted make a funnel --> 
kruskal.test(mydata$av.hydroid~mydata$Treatment)
#SIGNIFICANT --> p=0.0443 H=4.045, df=1

#-------------------------------------------------------------------
#BARE SPACE/ UN-USED SPACE

boxplot(mydata$av.bare.space~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of bare space")

zbare<-lm(mydata$av.bare.space~mydata$Treatment)
plot(zbare)
#residuals vs fitted make a funnel & not normal QQ
kruskal.test(mydata$av.bare.space~mydata$Treatment)
#NS --> p=0.5637, H=0.3333, df=1

#----------------------------------------------------------------
TOTAL COVER / USED SPACE

tcover<-(100-(mydata$av.bare.space))
boxplot(tcover~mydata$Treatment, cex.axis=0.8, boxwex=.5, varwidth=TRUE, xlab="Treatment", ylab="Percent cover of used space")
