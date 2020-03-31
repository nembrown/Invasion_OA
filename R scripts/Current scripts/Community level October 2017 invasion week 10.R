
######## read in week 8 as well

library(vegan)
library(sciplot)

## For sp div and richness use nudi.combined but for community level use nudi.separated


species.invasion.10 <- read.csv(file.choose())
row.names(species.invasion.10) <- species.invasion.10$ID
species.invasion.10 <- species.invasion.10[,-1]
head(species.invasion.10)
### had to do the -1 to take off botryllid as well as row

just.species.invasion.10<- species.invasion.10
View(just.species.invasion.10)

species.invasion.10<- invasion.exp.data.10[,c(12,14,22:25,28, 32:38,43:47)]

##species.invasion.10 richness##
species.invasion.10$richness <- specnumber(just.species.invasion.10)


##species.invasion.10 accumulation##

sp2 <- specaccum(just.species.invasion.10,"random")
plot(sp2)
summary(sp2)
coef(sp2)
plot(sp2, ci.type="poly", xlab="Number of Tiles",col="blue", lwd=2, ci.lty=0, ci.col="lightblue",ylab="Number of species.invasion.10")
boxplot(sp2, col="yellow", add=TRUE, pch="+")


?specaccum
##shnnon diversity##

species.invasion.10$shannon.diversity <- diversity(just.species.invasion.10,"shannon")


write.csv(species.invasion.10, "Oct.2016.Invasives.exp.community.results.all.4.csv")

##now mds and stuff##
##read in environment dataset##
compiled.data.invasion.10<-compiled.data.invasion[compiled.data.invasion$Week==10,]


compiled.data.invasion.10 <- read.csv(file.choose())
View(compiled.data.invasion.10)
row.names(compiled.data.invasion.10) <- compiled.data.invasion.10$ID
compiled.data.invasion.10 <- compiled.data.invasion.10[,-1]

head(compiled.data.invasion.10)
compiled.data.invasion.10$Combined.Treatment<-as.factor(compiled.data.invasion.10$Combined.Treatment) 
compiled.data.invasion.10$Combined.Treatment<-relevel(compiled.data.invasion.10$Combined.Treatment, "NoneAmbient")

compiled.data.invasion.10<-compiled.data.invasion.10[compiled.data.invasion.10$Week==12,]

compiled.data.invasion.10$Week<-as.numeric(compiled.data.invasion.10$Week)

levels(compiled.data.invasion.10$Combined.Treatment)

compiled.data.invasion.10.Elevated<-compiled.data.invasion.10[compiled.data.invasion.10$CO2=="Elevated", ]
mean(compiled.data.invasion.10.Elevated$av.pH, na.omit="T")


compiled.data.invasion.10.Ambient<-compiled.data.invasion.10[compiled.data.invasion.10$CO2=="Ambient", ]
mean(compiled.data.invasion.10.Ambient$av.pH, na.omit="T")

av.diff<-mean(compiled.data.invasion.10.Elevated$av.pH, na.omit="T") - mean(compiled.data.invasion.10.Ambient$av.pH, na.omit="T")

all.data<-cbind(just.species.invasion.10,compiled.data.invasion.10)
head(all.data)

compiled.data.invasion.10$hydrogen.concentration<-Invasives.exp.data.12$hydrogen.concentration

# Make an NMDS plotting function
nmds_plot <- function(dist.matrix, colorby){
  m <- metaMDS(dist.matrix) #NMDS object
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col="grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Invasives  CO2", legend=levels(colorby), col=cols, pch = shapesies)
}


#Bray-Curtis, Jaccard, Kulczynski (most common ones)
dist <- vegdist(species.invasion.10, method = "bray") #try other methods
coldiss(dist)
nmds_plot(dist, colorby=compiled.data.invasion.10$Combined.Treatment)


###CONSTRAINED Ordination
head(species.invasion.10)
species.invasion.10 <- read.csv(file.choose())
row.names(species.invasion.10) <- compiled.data.invasion.10$ID
species.invasion.10 <- species.invasion.10[,-7]
head(species.invasion.10)

standardized.species.invasion.10<-decostand(species.invasion.10, method="total", MARGIN=2)
hellinger.species.invasion.10<-decostand(species.invasion.10, method="hellinger")
is.na(species.invasion.10)

capscale_plot.low<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2.low#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Invasives  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}

capscale_plot.control<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2.control#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Invasives  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}


?capscale
capscale_plot_invasives<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.blue4#vector of colors needed
  shapesies<-c( 2,16,2,16)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=2)
  legend("topright", title ="Invasives  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}

cbbPalette.blue4<- c( "#666666",   "#666666","#619CFF", "#619CFF")


compiled.data.invasion.10$botryllid<-invasion.exp.data.10$botryllid
head(standardized.species.invasion.10)

m_invasives<-capscale(standardized.species.invasion.10 ~ hydrogen.concentration*botryllid,compiled.data.invasion.10 , distance="bray")
capscale_plot_invasives(m_invasives, colorby=compiled.data.invasion.10$Factor)
adonis(standardized.species.invasion.10 ~ hydrogen.concentration*botryllid, method="bray", permutations = 9999, data=compiled.data.invasion.10)

m_invasives<-capscale(standardized.species.invasion.10 ~ hydrogen.concentration*Invasives,compiled.data.invasion.10 , distance="bray")
capscale_plot_invasives(m_invasives, colorby=compiled.data.invasion.10$Factor)
adonis(standardized.species.invasion.10 ~hydrogen.concentration*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)


dist.standardized.species.invasion.10 <- vegdist(standardized.species.invasion.10, method = "bray")

mod.Invasives<-betadisper(dist.standardized.species.invasion.10, compiled.data.invasion.10$Invasives, type="centroid")
bd<-anova(mod.Invasives)

mod.CO2<-betadisper(dist.standardized.species.invasion.10, compiled.data.invasion.10$Treatment, type="centroid")
anova(mod.CO2)



mod.Invasives<-betadisper(dist.standardized.species.invasion.10, compiled.data.invasion.10$Combined.Treatment, type="centroid")
anova(mod.Invasives)
plot(mod.Invasives)
bd<-betadisper(dist[[3]],groups)



# Now, we'll calculate the Jaccard index and its partitions of turnover and nestedness. We can calculate Sorensen index instead by using the argument     index.family="sorensen"    .
install.packages("betapart")
library(betapart)
head(compiled.data.invasion.10)

presabs<-ifelse(standardized.species.invasion.10>0,1,0)
dist.agg.stand.pair.jacc<-beta.pair(presabs, index.family="jaccard")
bd<-betadisper(dist.agg.stand.pair.jacc[[3]],compiled.data.invasion.10$Factor)
bd.nestedness<-betadisper(dist.agg.stand.pair.jacc[[2]],compiled.data.invasion.10$Factor)
bd.turnover<-betadisper(dist.agg.stand.pair.jacc[[1]],compiled.data.invasion.10$Factor)


plot(bd)
anova(bd)
boxplot(bd)

plot(bd.nestedness)
anova(bd.nestedness)
boxplot(bd.nestedness)


plot(bd.turnover)
anova(bd.turnover)
boxplot(bd.turnover)
############# 
dist.agg.stand.pair.bray<-bray.part(standardized.species.invasion.10)
bd.bray<-betadisper(dist.agg.stand.pair.bray[[3]],compiled.data.invasion.10$Combined.Treatment)
bd.nestedness.bray<-betadisper(dist.agg.stand.pair.bray[[2]],compiled.data.invasion.10$Combined.Treatment)
bd.turnover.bray<-betadisper(dist.agg.stand.pair.bray[[1]],compiled.data.invasion.10$Combined.Treatment)


plot(bd.bray)
anova(bd.bray)
boxplot(bd.bray)

plot(bd.nestedness.bray)
anova(bd.nestedness.bray)
boxplot(bd.nestedness.bray)

plot(bd.turnover.bray)
anova(bd.turnover.bray)
boxplot(bd.turnover.bray)



############ 
############# 
dist.agg.stand.pair.bray.co2<-bray.part(standardized.species.invasion.10)
bd.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[3]],compiled.data.invasion.10$CO2)
bd.nestedness.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[2]],compiled.data.invasion.10$CO2)
bd.turnover.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[1]],compiled.data.invasion.10$CO2)


plot(bd.bray.co2)
anova(bd.bray.co2)
boxplot(bd.bray.co2)

plot(bd.nestedness.bray.co2)
anova(bd.nestedness.bray.co2)
boxplot(bd.nestedness.bray.co2)

plot(bd.turnover.bray.co2)
anova(bd.turnover.bray.co2)
boxplot(bd.turnover.bray.co2)

####### non standardized
dist.agg.pair.bray.co2.unstand<-bray.part(species.invasion.10.agg)
bd.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[3]],compiled.data.invasion.10$CO2)
bd.nestedness.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[2]],compiled.data.invasion.10$CO2)
bd.turnover.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[1]],compiled.data.invasion.10$CO2)


plot(bd.bray.co2.unstand)
anova(bd.bray.co2.unstand)
boxplot(bd.bray.co2.unstand)

plot(bd.nestedness.bray.co2.unstand)
anova(bd.nestedness.bray.co2.unstand)
boxplot(bd.nestedness.bray.co2.unstand)

plot(bd.turnover.bray.co2.unstand)
anova(bd.turnover.bray.co2.unstand)
boxplot(bd.turnover.bray.co2.unstand)

#########
dist.agg.stand.pair.bray.Invasives<-bray.part(standardized.species.invasion.10)
bd.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[3]],compiled.data.invasion.10$Invasives)
bd.nestedness.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[2]],compiled.data.invasion.10$Invasives)
bd.turnover.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[1]],compiled.data.invasion.10$Invasives)


plot(bd.bray.Invasives)
anova(bd.bray.Invasives)
boxplot(bd.bray.Invasives)

plot(bd.nestedness.bray.Invasives)
anova(bd.nestedness.bray.Invasives)
boxplot(bd.nestedness.bray.Invasives)

plot(bd.turnover.bray.Invasives)
anova(bd.turnover.bray.Invasives)
boxplot(bd.turnover.bray.Invasives)



##### not standardized Invasives
dist.agg.pair.bray.Invasives.unstand<-bray.part(species.invasion.10.agg)
bd.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[3]],compiled.data.invasion.10$Invasives)
bd.nestedness.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[2]],compiled.data.invasion.10$Invasives)
bd.turnover.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[1]],compiled.data.invasion.10$Invasives)


plot(bd.bray.Invasives.unstand)
anova(bd.bray.Invasives.unstand)
boxplot(bd.bray.Invasives.unstand)

plot(bd.nestedness.bray.Invasives.unstand)
anova(bd.nestedness.bray.Invasives.unstand)
boxplot(bd.nestedness.bray.Invasives.unstand)

plot(bd.turnover.bray.Invasives.unstand)
anova(bd.turnover.bray.Invasives.unstand)
boxplot(bd.turnover.bray.Invasives.unstand)

# To get the pairwise Jaccard index turnover partition between communities, type: dist[[1]]. To get nestedness partition, type: dist[[2]]. To get all beta diversity: dist[[3]].

# If we want to compare the beta diversities of communities aggregated by the treatments of "undisturbed" and "disturbed", we can use "betadisper" analysis.

bd<-betadisper(dist[[3]],groups)

plot(bd)













require(cowplot)
theme_set(theme_classic())

plot_grid(m.8to12.pH.stand.plot, m.8to12.pH.stand.plot, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))

par(mfrow=c(1,1))
###standardized species.invasion.10 agg might be the way to go....? 


##### 

capscale_plot.none(m, colorby=compiled.data.invasion.10$Combined.Treatment)
capscale_plot.low(m, colorby=compiled.data.invasion.10$Combined.Treatment)
capscale_plot.high(m, colorby=compiled.data.invasion.10$Combined.Treatment)

adonis(standardized.species.invasion.10 ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)


adonis(standardized.species.invasion.10 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)
#can do by CO2 or by min.6.pH

adonis(standardized.species.invasion.10 ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)


adonis(standardized.species.invasion.10 ~ Invasives*min.6.pH, method="bray", permutations = 9999, data=compiled.data.invasion.10)



### Betadisperion

dist.standardized.6 <- vegdist(standardized.species.invasion.10, method = "bray")




mod.CO2<-betadisper(dist.standardized.6, compiled.data.invasion.10$CO2, type="centroid")
anova(mod.CO2)


mod.pH<-betadisper(dist.standardized.6, compiled.data.invasion.10$min.10.pH, type="centroid")
anova(mod.pH)

mod.Invasives<-betadisper(dist.standardized.6, compiled.data.invasion.10$Invasives, type="centroid")
anova(mod.Invasives)

mod.Combined<-betadisper(dist.standardized.6, compiled.data.invasion.10$Combined.Treatment, type="centroid")
anova(mod.Combined)

pmod <- permutest(mod.Combined, permutations = 99, pairwise = TRUE)


(mod.HSD <- TukeyHSD(mod.Combined))
plot(mod.HSD)








######## Sessile only
head(standardized.species.invasion.10)
standardized.species.invasion.10.sessile<-standardized.species.invasion.10[,-(10:12)]
standardized.species.invasion.10.sessile<-standardized.species.invasion.10.sessile[,-18]
head(standardized.species.invasion.10.sessile)

m.sessile<-capscale(standardized.species.invasion.10.sessile ~ CO2*Invasives,compiled.data.invasion.10 , distance="bray")
capscale_plot(m.sessile, colorby=compiled.data.invasion.10$Combined.Treatment)
adonis(standardized.species.invasion.10.sessile ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)



###pa presence absence
m.pa<-capscale(standardized.species.invasion.10 ~ CO2*Invasives,compiled.data.invasion.10 , distance="jaccard")

capscale_plot(m.pa, colorby=compiled.data.invasion.10$Combined.Treatment)



########Correlation analysis #################################################

species.invasion.10.sum<-apply(species.invasion.10, 2, sum)
species.invasion.10.sorted<-species.invasion.10[,order(species.invasion.10.sum, decreasing=TRUE)]
head(species.invasion.10.sorted)
species.invasion.10.small<-species.invasion.10.sorted[,1:15]
species.invasion.10.very.small<-species.invasion.10.sorted[,1:7]
species.invasion.10.top.ten<-species.invasion.10.sorted[,1:10]
species.invasion.10.top.five<-species.invasion.10.sorted[,1:5]


head(species.invasion.10.very.small)

species.invasion.10.small.std<-decostand(species.invasion.10.small, method="total", MARGIN=2)
species.invasion.10.small.t<-t(species.invasion.10.small.std)
head(species.invasion.10.small.t)

species.invasion.10.t.kmeans.casc<-cascadeKM(species.invasion.10.small.t, inf.gr=2, sup.gr=14, iter=100, criterion="calinski")
plot(species.invasion.10.t.kmeans.casc, sortg=TRUE)
species.invasion.10.t.kmeans.casc$results
species.invasion.10.t.kmeans.casc$partition

clusters<-species.invasion.10.t.kmeans.casc$partition[,1]
clusters

species.invasion.10.kendall.global<-kendall.global(species.invasion.10.small.hel, clusters)
species.invasion.10.kendall.global

species.invasion.10.kendall.post<-kendall.post(species.invasion.10.small.hel, clusters, nperm=9999)
species.invasion.10.kendall.post


################### Standardized and small 

species.invasion.10.small.t<-t(species.invasion.10.top.five)
head(species.invasion.10.small.t)

species.invasion.10.standardized.t<-decostand(species.invasion.10.small.t, method="total", MARGIN=1)

species.invasion.10.standardized.t.dist<-dist(species.invasion.10.standardized.t)

coldiss(species.invasion.10.standardized.t.dist, diag=TRUE)


############# exact same result if you standardize before or after so that's good. 

species.invasion.10.standardized.small<-decostand(species.invasion.10.very.small, method="total", MARGIN=2)

species.invasion.10.standardized.small.t<-t(species.invasion.10.standardized.small)

species.invasion.10.standardized.t.dist<-dist(species.invasion.10.standardized.small.t)

coldiss(species.invasion.10.standardized.t.dist, diag=TRUE)



dist.small <- vegdist(species.invasion.10.small.t, method = "bray")

coldiss(dist.small, diag=TRUE)


##################
head(species.invasion.10)

scaled.species.invasion.10<-scale(species.invasion.10)
fix(species.invasion.10)
head(scaled.species.invasion.10)

standardized.species.invasion.10<-decostand(species.invasion.10, method="total", MARGIN=2)



head(standardized.species.invasion.10.2)
wisconsin.species.invasion.10<- wisconsin(species.invasion.10)
head(wisconsin.species.invasion.10)
head(species.invasion.10)




#species.invasion.10.2<-species.invasion.10.2[, -11]

head(compiled.data.invasion.10)

adonis(species.invasion.10.2 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)

adonis(standardized.species.invasion.10.2 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)


adonis(standardized.species.invasion.10.3 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)
adonis(standardized.species.invasion.10.3 ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)
########


adonis(standardized.species.invasion.10.2 ~ CO2*Invasives, method="jaccard", permutations = 9999, data=compiled.data.invasion.10)




adonis(species.invasion.10 ~ min.6.pH, method="bray", permutations = 9999, data=compiled.data.invasion.10)

adonis(species.invasion.10 ~ av.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.10)



dist.standardized.2 <- vegdist(standardized.species.invasion.10.2, method = "bray")


mod.CO2<-betadisper(dist.standardized.2, compiled.data.invasion.10$CO2, type="centroid")
anova(mod.CO2)

mod.Invasives<-betadisper(dist.standardized.2, compiled.data.invasion.10$Invasives, type="centroid")
anova(mod.Invasives)

mod.Combined<-betadisper(dist.standardized.2, compiled.data.invasion.10$Combined.Treatment, type="centroid")
anova(mod.Combined)


######### SIMPER - type analysis?
#### simper is contribution to dissimilarity between sites
#### indv val is indicative of fidelity of that species.invasion.10 to a particular site
#### 1 - mv abunda - 2- indval? 

##### do I do CO2 or Invasives as the measure?? Too complicated to do both...?? 

sim <- with(compiled.data.invasion.10, simper(species.invasion.10,Invasives))
summary(sim)

simper(species.invasion.10,compiled.data.invasion.10$CO2, permutations = 0, trace = FALSE)


library(mvabund)
?decostand
tasmvabund <- mvabund(Tasmania$copepods)
plot(tasmvabund ~ treatment, col = as.numeric(block))



species.invasion.10.mvabund<-mvabund(standardized.species.invasion.10)
plot(species.invasion.10.mvabund ~ compiled.data.invasion.10$Combined.Treatment)
#### not sure how useful this is. 

##### Indval
###"The indval approach looks for species.invasion.10 that are
#both necessary and sufficient, i.e. if you find that species.invasion.10 you should be in that type,
#and if you are in that type you should find that species.invasion.10
install.packages("labdsv")
library(labdsv)
iva<-indval(standardized.species.invasion.10,compiled.data.invasion.10$Combined.Treatment)


iva$relfrq
iva$relabu
iva$indval

gr<- iva$maxcls[iva$pval<=0.05]
iv<- iva$indcls[iva$pval<=0.05]
pv<- iva$pval[iva$pval<=0.05]
fr<-apply(standardized.species.invasion.10>0, 2, sum)[iva$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg<-fidg[order(fidg$group,-fidg$indval), ]
fidg


isamic(standardized.species.invasion.10,compiled.data.invasion.10$CO2, sort=TRUE)



####


species.invasion.10.sum<-apply(species.invasion.10, 2, sum)
species.invasion.10.sorted<-species.invasion.10[,order(species.invasion.10.sum, decreasing=TRUE)]
species.invasion.10.3<-species.invasion.10.sorted[,1:22]
head(species.invasion.10.3)

standardized.species.invasion.10.3<-decostand(species.invasion.10.3, method="total", MARGIN=2)

iva<-indval(standardized.species.invasion.10.3,compiled.data.invasion.10$CO2)



#####
library(vegan)
set.seed(2)


treat=c(rep("High",19),rep("high",20), rep("Low", 20))
ordiplot(example_NMDS,type="n")
ordiellipse(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(example_NMDS,display="species.invasion.10",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)


groups = c(rep("group1", 36), rep("group2", 64))

orditorp(example_NMDS, display = "species.invasion.10")
