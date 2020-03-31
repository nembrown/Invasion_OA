
cbbPalette.all<- c( "#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")
cbbPalette.3<- c( "#F8766D", "#00BA38", "#619CFF")

cbbPalette.2<- c( "#F8766D","#619CFF")


head(invasion.exp.data.16)
head(invasion.exp.data.16)
combined.data

rbind(invasion.exp.data.16[,c(1:24)]

head(invasion.exp.data.16)      
         
species12AIRAbsent<-invasion.exp.data.16[invasion.exp.data.16$Treatment=="AIRAbsent",c(14,22:25,28, 32:38,43:47)]
species12AIRPresent<-invasion.exp.data.16[invasion.exp.data.16$Treatment=="AIRPresent",c(14,22:25,28, 32:38,43:47)]
species12CO2Present<-invasion.exp.data.16[invasion.exp.data.16$Treatment=="CO2Present",c(14,22:25,28, 32:38,43:47)]
species12CO2Absent<-invasion.exp.data.16[invasion.exp.data.16$Treatment=="CO2Absent",c(14,22:25,28, 32:38,43:47)]

      
      ##### week 12 
      
AIRAbsent<-t(species12AIRAbsent)
AIRPresent<-t(species12AIRPresent)
CO2Present<-t(species12CO2Present)
CO2Absent<-t(species12CO2Absent)

      library(vegan)
      
      
      AIRAbsent<-decostand(AIRAbsent, method="pa")
      AIRPresent<-decostand(AIRPresent, method="pa")
      CO2Present<-decostand(CO2Present, method="pa")
      CO2Absent<-decostand(CO2Absent, method="pa")
    
      
      list.invasion.exp<-list(AIRAbsent=AIRAbsent, AIRPresent=AIRPresent, CO2Present=CO2Present, CO2Absent=CO2Absent)
      head(list.invasion.exp)
      str(list.invasion.exp)
      
      #list.invasion.exp<-as.matrix(list.invasion.exp)
      
      list.invasion.exp.CO2<-list(AIRPresent=AIRPresent, CO2Absent=CO2Absent,  NONE.CO2=NONE.CO2)
      list.invasion.exp.AIR<-list(AIRAbsent=AIRAbsent,  CO2Present=CO2Present,  NONE.AIR=NONE.AIR)
      
      list.invasion.exp.CO2.just2<-list(AIRPresent=AIRPresent,  NONE.AIR=NONE.AIR)
      
      out.raw <- iNEXT(list.invasion.exp, datatype="abundance", endpoint=20, q=0)
      
      list.invasion.exp.abun<-as.abucount(list.invasion.exp)
      
      
library(iNEXT)

out.raw <- iNEXT(list.invasion.exp, datatype="incidence_raw", endpoint=24, q=0)
out.raw.1 <- iNEXT(list.invasion.exp, datatype="incidence_raw", endpoint=24, q=1)
out.raw.2 <- iNEXT(list.invasion.exp, datatype="incidence_raw", endpoint=24, q=2)

out.raw.CO2 <- iNEXT(list.invasion.exp.CO2, datatype="incidence_raw", endpoint=24)
out.raw.AIR <- iNEXT(list.invasion.exp.AIR, datatype="incidence_raw", endpoint=24)
out.raw.CO2.just2 <- iNEXT(list.invasion.exp.CO2.just2, datatype="incidence_raw", endpoint=24, q=2)




#Species number by sample
g1.richness<-ggiNEXT(out.raw, type=1, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))
g1.diversity<-ggiNEXT(out.raw.1, type=1, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))
g1.simpson<-ggiNEXT(out.raw.2, type=1, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))


g2.richness<-ggiNEXT(out.raw, type=2, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))+  scale_linetype_manual(values=c(1, 2, 3, 2))
g2.diversity<-ggiNEXT(out.raw.1, type=2, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))
g2.simpson<-ggiNEXT(out.raw.2, type=2, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))

g3.richness<-ggiNEXT(out.raw, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))+  scale_linetype_manual(values=c(1, 2, 3, 2))
g3.diversity<-ggiNEXT(out.raw.1, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))
g3.simpson<-ggiNEXT(out.raw.2, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))


g3.richness<-ggiNEXT(out.raw, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4) + scale_shape_manual(values=c(19,1, 17,2,15,0))
g3.richness<- g3.richness + ylab(expression("Species richness across mesocosms"))+ theme_bw()
g3.richness <- g3.richness + theme(text = element_text(size=18), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
g3.richness<-  g3.richness +theme(legend.text = element_text(colour="black", size = 18))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
g3.richness<- g3.richness+ theme(legend.position="none")
g3.richness

g1.richness<-ggiNEXT(out.raw, type=1, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4) + scale_shape_manual(values=c(19,1, 17,2,15,0))
g1.richness<- g1.richness + ylab(expression("Species richness across mesocosms"))+ theme_bw()
g1.richness <- g1.richness + theme(text = element_text(size=18), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
g1.richness<-  g1.richness +theme(legend.text = element_text(colour="black", size = 18))+ theme(legend.title = element_text(colour="black", size=18))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
g1.richness<- g1.richness#+ theme(legend.position="none")
g1.richness



### is this for #12 or for all weeks??? 



plot.richness.12<- ggplot(invasion.exp.data.12, aes(x=invasion.exp.data.12$min.10.pH, y=(invasion.exp.data.12$richness), colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives), shape=factor(Treatment))) + guides(fill=FALSE) + geom_smooth(aes(col=Invasives, fill=Invasives), alpha=0.15,size=1.5, method="lm")+ scale_fill_manual(values=cbbPalette.all)
plot.richness.12<- plot.richness.12 + theme_bw() + xlab(bquote('Minimum 10th percentile pH')) + ylab(expression("Species richness per mesocosm"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ scale_y_continuous(limits = c(0, 15))
plot.richness.12<- plot.richness.12 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.richness.12<- plot.richness.12+ scale_colour_discrete(name = "Invasives")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.richness.12 <- plot.richness.12+ theme(legend.position="none")+ scale_shape_manual(values=c(19,17, 19, 17, 19, 17))
plot.richness.12




plot_grid(plot.richness.12, plot.stand.biomass.12, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))


plot_grid(g1.richness,g2.richness, g3.richness,g1.diversity, g2.diversity, g3.diversity, g1.simpson, g2.simpson, g3.simpson, ncol=3, align='h', labels=c('(a)', '(b)', label_size=12))

plot_grid(g.richness, plot.richness.12, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))

out.raw$DataInfo
out.raw$iNextEst

out.raw

ggiNEXT(out.raw, type=3, color.var="site", facet.var = "site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4) + scale_shape_manual(values=c(19,1, 17,2,15,0))

#### confused b/c it seems like abundance should be included here... 

list.invasion.exp.freq<-as.incfreq(list.invasion.exp)


?estimateD

estimateD(list.invasion.exp, datatype="incidence_raw",base="coverage", level=0.95)

estimateD(x, datatype = "abundance", base = "size", level = NULL,
          conf = 0.95)


ChaoRichness(list.invasion.exp, datatype = "incidence_raw")
ChaoShannon(list.invasion.exp, datatype = "incidence_raw",  (transform=TRUE))
ChaoShannon(list.invasion.exp, datatype = "incidence_raw")
ChaoSimpson(list.invasion.exp, datatype = "incidence_raw")

out.raw
ggiNEXT(out.raw, type=2, color.var="site")+ scale_colour_manual(values=cbbPalette.blue4) +  scale_fill_manual(values=cbbPalette.blue4)+scale_size_manual(values=c(1,4))


g5 <- g + theme_bw() + theme(legend.position = "bottom")
g6 <- g + theme_classic() + theme(legend.position = "bottom")
grid.arrange(g5, g6, ncol=2)

plot.out.raw.CO2<-ggiNEXT(out.raw.CO2, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.3) +  scale_fill_manual(values=cbbPalette.3)+ylim(0,20)
plot.out.raw.AIR<-ggiNEXT(out.raw.AIR, type=3, color.var="site")+ scale_colour_manual(values=cbbPalette.3) +  scale_fill_manual(values=cbbPalette.3)+ylim(0,20)

theme_set(theme_classic())

library(cowplot)

plot_grid(plot.out.raw.AIR, plot.out.raw.CO2, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))




#sample coverage by sampling units
ggiNEXT(out.raw, type=2)

#species diversity by sample coverage
ggiNEXT(out.raw, type=3)
values=cbbPalette.all)
cbbPalette.blue4<- c( "#F8766D","#F8766D", "#00BA38", "#00BA38", "#619CFF", "#619CFF")


plot.accum.iNEXT<-ggiNEXT(out.raw, type=3)  + scale_colour_manual(values=cbbPalette.blue4) + scale_fill_manual(values=cbbPalette.blue4) + scale_fill_manual(values=cbbPalette.blue4)
plot.accum.iNEXT<- plot.accum.iNEXT + labs(x = "Sample completeness", y="Species richness per tank")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))
plot.accum.iNEXT<-plot.accum.iNEXT + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.accum.iNEXT <- plot.accum.iNEXT +  theme(legend.position="none")+ scale_shape_manual(values=c(19,1, 17,2,15,0))+theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.accum.iNEXT

data(ant)
head(ant$h50m)
list.invasion.exp.abucount<-lapply(list.invasion.exp, as.abucount)
out.raw <- iNEXT(list.invasion.exp.abucount, datatype="abundance", endpoint=30)
ggiNEXT(out.raw)

data(ciliates)
lapply(ciliates, as.abucount)



###################

###### iNExt
install.packages("iNEXT")
library(iNEXT)

data(ciliates)
head(ciliates$CentralNamibDesert)
lapply(ciliates, as.abucount)
## example for abundance based data (list of vector)
str(spider)
head(spider$Girdled)
out1 <- iNEXT(spider, q=0, datatype="abundance")
out1$DataInfo # showing basic data information.
out1$AsyEst # showing asymptotic diversity estimates.
out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
## example for abundance based data (data.frame)
data(bird)
head(bird)
out2 <- iNEXT(bird, q=0, datatype="abundance")
ggiNEXT(out2)

head(spider)
out <- iNEXT(spider, q=c(0, 1, 2), datatype="abundance", endpoint=500)
# Sample-size-based R/E curves, separating by "site""
ggiNEXT(out, type=1, facet.var="site")
## Not run:
# Sample-size-based R/E curves, separating by "order"
ggiNEXT(out, type=1, facet.var="order")
# display black-white theme
ggiNEXT(out, type=1, facet.var="order", grey=TRUE)
## End(Not run)
The argument facet.var="site" in ggiNEXT function creates a separate plot for each site as shown below:
  # Sample-size-based R/E curves, separating by "site""
  ggiNEXT(out, type=1, facet.var="site")


##Spec pool??? 

data(dune.env)
attach(compiled.data.12)



pool <- specpool(species.12, Treatment)
pool
op <- par(mfrow=c(1,2))
boxplot(specnumber(species.12) ~ Treatment, col="hotpink", border="cyan3",
        notch=TRUE)
boxplot(specnumber(species.12)/specpool2vect(pool) ~ Treatment, col="hotpink",
        border="cyan3", notch=TRUE)
par(op)
data(BCI)
## Accumulation model
pool <- poolaccum(species.12)
summary(pool, display = "chao")
plot(pool)
## Quantitative model
estimateR(species.12[1:5,])

species.12.AIRAbsent <- specpool(just.species.12.AIRAbsent,"random")
species.12.AIRPresent <- specpool(just.species.12.AIRPresent,"random")
species.12.CO2Present <- specpool(just.species.12.CO2Present,"random")
species.12.CO2Absent <- specpool(just.species.12.CO2Absent,"random")
species.12.NONE.AIR <- specpool(just.species.12.NONE.AIR,"random")
species.12.NONE.CO2 <- specpool(just.species.12.NONE.CO2,"random")


df.12.AIRAbsent <- data.frame(species.12.AIRAbsent$richness,species.12.AIRAbsent$sites,species.12.AIRAbsent$sd)
df.12.AIRPresent <- data.frame(species.12.AIRPresent$richness,species.12.AIRPresent$sites,species.12.AIRPresent$sd)
df.12.CO2Present <- data.frame(species.12.CO2Present$richness,species.12.CO2Present$sites,species.12.CO2Present$sd)
df.12.CO2Absent <- data.frame(species.12.CO2Absent$richness,species.12.CO2Absent$sites,species.12.CO2Absent$sd)
df.12.NONE.AIR <- data.frame(species.12.NONE.AIR$richness,species.12.NONE.AIR$sites,species.12.NONE.AIR$sd)
df.12.NONE.CO2 <- data.frame(species.12.NONE.CO2$richness,species.12.NONE.CO2$sites,species.12.NONE.CO2$sd)

Make sure all colnames are the same
df.12.AIRAbsent$Invasives<-"High"
df.12.AIRPresent$Invasives<-"High"
df.12.CO2Present$Invasives<- "Low"
df.12.CO2Absent$Invasives<- "Low"
df.12.NONE.AIR$Invasives <-"None"
df.12.NONE.CO2$Invasives <- "None"

df.12.AIRAbsent$CO2<-"AIR"
df.12.AIRPresent$CO2<-"CO2"
df.12.CO2Present$CO2<- "AIR"
df.12.CO2Absent$CO2<- "CO2"
df.12.NONE.AIR$CO2 <-"AIR"
df.12.NONE.CO2$CO2<- "CO2"



colnames(df.12.AIRAbsent) <- c("richness","sites","sd","Invasives", "CO2") 
colnames(df.12.AIRPresent) <- c("richness","sites","sd","Invasives", "CO2")
colnames(df.12.CO2Present) <- c("richness","sites","sd","Invasives", "CO2")
colnames(df.12.CO2Absent) <- c("richness","sites","sd","Invasives", "CO2")
colnames(df.12.NONE.AIR) <- c("richness","sites","sd","Invasives", "CO2")
colnames(df.12.NONE.CO2) <- c("richness","sites","sd","Invasives", "CO2")
Add new column to each df as an identifier


rbind everything together

df.12 <- rbind(df.12.AIRAbsent, df.12.AIRPresent, df.12.CO2Present, df.12.CO2Absent, df.12.NONE.AIR, df.12.NONE.CO2)




plot.pool<-ggplot(df.12, aes(x=sites, y=richness, color = Invasives, shape=CO2, linetype=CO2)) + geom_line(size=1.15) + geom_point(size=0)  + geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd), width=0, size=1.15, lty=1) + scale_colour_manual(values=cbbPalette.all)
plot.pool<- plot.pool + labs(x = "Number of mesocosms samples", y="Species Richness")+ theme_bw()+ theme(text = element_text(size=16), axis.text = element_text(size=16)) + theme(axis.title.y = element_text(angle=90))+ scale_y_continuous(limits = c(0, 15))
plot.pool<-plot.pool + theme(legend.text = element_text(colour="black", size = 16))+theme(legend.title = element_text(colour="black", size=16), legend.position=c(.15,.80))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )+ theme(legend.key = element_blank())+theme(legend.key.size = unit(0.3, "cm"))
plot.pool <- plot.pool + theme(legend.position="none")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.pool


plot_grid(plot.accum, plot.richness.12, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))

