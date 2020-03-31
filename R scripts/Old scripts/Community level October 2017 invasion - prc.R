
######## read in week 8 as well

library(vegan)
library(sciplot)

## For sp div and richness use nudi.combined but for community level use nudi.separated

### not this
species.invasion <- read.csv(file.choose())
row.names(species.invasion) <- species.invasion$ID
species.invasion <- species.invasion[,-1]
head(species.invasion)

species.invasion<- invasion.exp.data[,c(12,14,22:25,28, 32:38,43:47)]
row.names(species.invasion) <- compiled.data.invasion$ID
just.species.invasion<- species.invasion


##species.invasion richness##
species.invasion$richness <- specnumber(just.species.invasion)



##shnnon diversity##

species.invasion$shannon.diversity <- diversity(just.species.invasion,"shannon")


write.csv(species.invasion, "Oct.2016.invasion.exp.community.results.all.4.csv")

##now mds and stuff##
##read in environment dataset##

compiled.data.invasion <- read.csv(file.choose())
View(compiled.data.invasion)
row.names(compiled.data.invasion) <- compiled.data.invasion$ID
compiled.data.invasion <- compiled.data.invasion[,-1]

head(compiled.data.invasion)
compiled.data.invasion$Factor<-as.factor(compiled.data.invasion$Factor) 
compiled.data.invasion$Factor<-relevel(compiled.data.invasion$Factor, "AIR, Absent")

compiled.data.invasion.16<-compiled.data.invasion[compiled.data.invasion$Week==16,]



compiled.data.invasion$hydrogen.concentration<-invasion.exp.data$hydrogen.concentration

all.data.invasion<-cbind(just.species.invasion,compiled.data.invasion)
head(all.data.invasion)


head(species.invasion)

all.data.invasion.16<-all.data.invasion[all.data.invasion$Week==16,]
all.data.invasion.14<-all.data.invasion[all.data.invasion$Week==14,]
all.data.invasion.12<-all.data.invasion[all.data.invasion$Week==12,]
all.data.invasion.10<-all.data.invasion[all.data.invasion$Week==10,]
all.data.invasion.8<-all.data.invasion[all.data.invasion$Week==8,]
all.data.invasion.6<-all.data.invasion[all.data.invasion$Week==6,]
all.data.invasion.4<-all.data.invasion[all.data.invasion$Week==4,]
all.data.invasion.2<-all.data.invasion[all.data.invasion$Week==2,]
all.data.invasion.0<-all.data.invasion[all.data.invasion$Week==0,]

species.invasion.nozero<-rbind(all.data.invasion.12[,c(2:19)],all.data.invasion.10[,c(2:19)], all.data.invasion.8[,c(2:19)], all.data.invasion.6[,c(2:19)], all.data.invasion.4[,c(2:19)], all.data.invasion.2[,c(2:19)])

head(species.invasion.8to12)


standardized.all.data.invasion.12<-decostand(all.data.invasion.12[,c(2:19)], method="total", MARGIN=2)
head(standardized.all.data.invasion.12)
standardized.all.data.invasion.10<-decostand(all.data.invasion.10[,c(2:19)], method="total", MARGIN=2)
standardized.all.data.invasion.8<-decostand(all.data.invasion.8[,c(2:19)], method="total", MARGIN=2)
standardized.all.data.invasion.6<-decostand(all.data.invasion.6[,c(2:19)], method="total", MARGIN=2)
standardized.all.data.invasion.4<-decostand(all.data.invasion.4[,c(2:19)], method="total", MARGIN=2)
standardized.all.data.invasion.2<-decostand(all.data.invasion.2[,c(2:19)], method="total", MARGIN=2)
standardized.all.data.invasion.0<-decostand(all.data.invasion.0[,c(2:19)], method="total", MARGIN=2)


###hellinger
head(standardized.all.data.invasion.hellinger12)
standardized.all.data.invasion.hellinger16<-decostand(all.data.invasion.16[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger14<-decostand(all.data.invasion.14[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger12<-decostand(all.data.invasion.12[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger10<-decostand(all.data.invasion.10[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger8<-decostand(all.data.invasion.8[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger6<-decostand(all.data.invasion.6[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger4<-decostand(all.data.invasion.4[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger2<-decostand(all.data.invasion.2[,c(2:19)], method="hellinger")
standardized.all.data.invasion.hellinger0<-decostand(all.data.invasion.0[,c(2:19)], method="hellinger")


species.invasion.standardized.all.data.invasion.hellinger.nozero<-rbind(standardized.all.data.invasion.hellinger16,standardized.all.data.invasion.hellinger14,standardized.all.data.invasion.hellinger12, standardized.all.data.invasion.hellinger10, standardized.all.data.invasion.hellinger8, standardized.all.data.invasion.hellinger6, standardized.all.data.invasion.hellinger4, standardized.all.data.invasion.hellinger2)

###################
compiled.data.invasion.no.zero<-compiled.data.invasion[compiled.data.invasion$Week!=0,]
head(compiled.data.invasion.no.zero)
compiled.data.invasion.no.zero$Week<-as.factor(compiled.data.invasion.no.zero$Week)
compiled.data.invasion.no.zero$Factor<-as.factor(compiled.data.invasion.no.zero$Factor)

compiled.data.invasion.no.zero$Factor<-relevel(compiled.data.invasion.no.zero$Factor, "AIR, Absent")


compiled.data.invasion.no.zero$Factor<-as.factor(compiled.data.invasion.no.zero$Factor)


##### no zero hellinger #this is for publication

species.invasion.standardized.all.data.invasion.hellinger.nozero<-plyr::rename(species.invasion.standardized.all.data.invasion.hellinger.nozero, c("alive.bot"="botryllus", "schizo"="schizoporella", "alive.mem"="membranipora", "alive.barn"="balanus"))


pyr_prc_no.zero.hellinger<- prc(response = species.invasion.standardized.all.data.invasion.hellinger.nozero, treatment = compiled.data.invasion.no.zero$Factor, time = compiled.data.invasion.no.zero$Week, scaling=2)
plot(pyr_prc_no.zero.hellinger,xlab="Week", lty=c(1,2,1), col = c(  "#666666",  "#619CFF", "#619CFF"), cex=0.8, lwd=3, type="b", pch=c(19,2, 19),  scaling=2) 
segments(y0 =1.71,  x1=8.8, x0=0.2, y1=1.71, cex=0.8, lwd=3,lty=2, col="#666666")
points(col="#666666", xy.coords(8.8, 1.71),pch=19)
points(col="#666666", xy.coords(0, 1.71),pch=19)

ID<-as.factor(compiled.data.invasion.no.zero$ID)
legpos=NA



require(vegan)
 
###stat = week - means what to permute within 	An integer vector or factor specifying the strata for permutation. If supplied, observations are permuted only within the specified strata.
#Assess only the significance of the first constrained eigenvalue 
anova(pyr_prc_no.zero.hellinger, strata = compiled.data.invasion.no.zero$Week, first=TRUE)



overall.CO2<-anova(pyr_prc_no.zero.hellinger.CO2, strata = compiled.data.invasion.no.zero$Week, first=TRUE)
#overall CO2 effect
overall.invasion<-anova(pyr_prc_no.zero.hellinger.Invasives, strata = compiled.data.invasion.no.zero$Week, first=TRUE)
# overall no invasion effect

noinvasion.anova<-anova(pyr_prc_no.zero.hellinger.noinvasion.CO2, strata = compiled.data.invasion.no.zero.noinvasion$Week, first=TRUE)

lowinvasion.anova<-anova(pyr_prc_no.zero.hellinger.lowinvasion.CO2, strata = compiled.data.invasion.no.zero.lowinvasion$Week, first=TRUE)
# CO2 effect in low invasion
highinvasion.anova<-anova(pyr_prc_no.zero.hellinger.highinvasion.CO2, strata = compiled.data.invasion.no.zero.highinvasion$Week, first=TRUE)
# overall CO2 effect in high invasion

p<-rbind(highinvasion.anova$`Pr(>F)`, lowinvasion.anova$`Pr(>F)`, noinvasion.anova$`Pr(>F)`, overall.invasion$`Pr(>F)`, overall.CO2$`Pr(>F)`)
colnames(p)<-c("pvalue", "type")


### can't quite figure this out
p[1:5, 2]<-c("highinvasion", "lowinvasion", "noinvasion", "overall.invasion", "overall.co2")

p$type<- factor(p$type,labels = c("highinvasion", "lowinvasion", "noinvasion", "overall.invasion", "overall.co2"))



p.adjust(p[,1], method = "bonferroni", n = length(p[,1]))

#### this works,.... 




p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")




install.packages("permute", repos="http://R-Forge.R-project.org")
library(permute)

control = how(plots = Plots(strata = ID, type = "free"), within = Within(type = "none"), 
              nperm = 199)
set.seed(1234)
permutations <- shuffleSet(nrow(species.invasion.standardized.all.data.invasion.hellinger.nozero), control = control)

mod_perm <- permutest(pyr_prc_no.zero.hellinger, permutations = permutations, first = TRUE)
mod_perm
#### this is just overall?? I can get that from ANOVA??? 
require(multcomp)
df <- data.frame(dose = compiled.data.invasion.no.zero$Factor, week = compiled.data.invasion.no.zero$Week)
prc(response = pyrifos, treatment = dose, time = week)


out_willi <- NULL

for (i in levels(week)) {
  take_spec <-species.invasion.standardized.all.data.invasion.hellinger.nozero[week == i, ]
  pca <- rda(take_spec)  # Compute PCA
  pca_scores <- scores(pca, display = "sites", choices = 1)  # scores of first principle component
  
  out_willi[[i]] <- summary(glht(aov(pca_scores ~ dose, data = df[week == i, ]), alternative = "t", linfct = mcp(dose = "Williams")))
}
# extract p-values
result <- lapply(out_willi, function(x) data.frame(comp = levels(df$dose)[-1], 
                                                   pval = x$test$pvalues, sig = x$test$pvalues < 0.05))
# shows the results of Dunnett-Test on PCA-scores for week 1:
result[["8"]]

##### not sure this is right??? might have something to do with why the scale on weeksis fucked.... 

############ 



### can I make the hline blue solid? + geom_hline(aes(x=2, yintercept=0), lty=2,size=1)
?prc

pyr_prc <- prc(response = species.invasion.standardized.all, treatment = compiled.data.invasion$Factor, time = compiled.data.invasion$Week)
plot(pyr_prc, lty=c(1,2,1,2,2), col = c( "#F8766D", "#F8766D", "#00BA38","#00BA38", "#619CFF"), cex=0.8, lwd=2.5, type="b", pch=19)

summary(pyr_prc)
anova(pyr_prc, strata = compiled.data.invasion$Week, first=TRUE)

pyr_prc_CO2 <- prc(response = species.invasion.standardized.all, treatment = compiled.data.invasion$Invasives, time = compiled.data.invasion$Week)
plot(pyr_prc_CO2, lty=c(1,2,1,2,2), col = c( "#F8766D", "#00BA38", "#619CFF"), cex=0.8, lwd=2.5, type="b", pch=19)


#### Standardized ... unclear what this is showing differently
standardized.species.invasion<-decostand(species.invasion, method="total", MARGIN=2)

pyr_prc_standardized <- prc(response = species.invasion.standardized.all, treatment = compiled.data.invasion$Factor, time = compiled.data.invasion$Week)
plot(pyr_prc_standardized, legpos=NA, lty=c(1,2,1,2,2), col = c( "#F8766D", "#F8766D", "#00BA38","#00BA38", "#619CFF"), cex=0.8, lwd=2.5, type="b", pch=19)


pyr_prc_standardized.nozero <- prc(response = species.invasion.standardized.all.nozero, treatment = compiled.data.invasion.no.zero$Factor, time = compiled.data.invasion.no.zero$Week)
plot(pyr_prc_standardized.nozero, legpos=NA, lty=c(1,2,1,2,2), col = c( "#F8766D", "#F8766D", "#00BA38","#00BA38", "#619CFF"), cex=0.8, lwd=2.5, type="b", pch=19, xlab="Week", scaling=2)

##### Need to really carefully think this through ... standardizing each species.invasion ... so each species.invasion not on same scale ... but also gives
# lot of weight to the rare species.invasion ... 

#### Also needs to be standardized by tile too!! Maybe first?? **** b/c only species.invasion ... not bare space... 

###Standardized by tile first then by species.invasion
pyr_prc_standardized.both <- prc(response = species.invasion.standardized.all.both, treatment = compiled.data.invasion$Factor, time = compiled.data.invasion$Week)
plot(pyr_prc_standardized.both, lty=c(1,2,1,2,2), col = c( "#F8766D", "#F8766D", "#00BA38","#00BA38", "#619CFF"), cex=0.8, lwd=2.5, type="b", pch=19)



summary(pyr_prc)
anova(pyr_prc, strata = compiled.data.invasion$Week, first=TRUE)



###
cbbPalette.all.2<- c( "#F8766D", "#F8766D", "#00BA38" , "#00BA38", "#619CFF", "#619CFF")
cbbPalette.all.2.none<- c( "#FFFFFF", "#FFFFFF", "#FFFFFF" , "#FFFFFF", "#619CFF", "#619CFF")
cbbPalette.all.2.low<- c( "#FFFFFF", "#FFFFFF", "#00BA38" , "#00BA38","#FFFFFF" , "#FFFFFF")
cbbPalette.all.2.high<- c( "#F8766D", "#F8766D", "#FFFFFF", "#FFFFFF", "#FFFFFF" , "#FFFFFF")


cbbPalette.all.2.control<- c( "#F8766D", "#FFFFFF", "#00BA38" , "#FFFFFF", "#619CFF", "#FFFFFF")

##text(m, labels = row.names(as.matrix(dist.matrix)), cex = 0.6) 

# Make an NMDS plotting function
nmds_plot <- function(dist.matrix, colorby){
  m <- metaMDS(dist.matrix) #NMDS object
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col="grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="invasion  CO2", legend=levels(colorby), col=cols, pch = shapesies)
}


#Bray-Curtis, Jaccard, Kulczynski (most common ones)
dist <- vegdist(species.invasion, method = "bray") #try other methods
coldiss(dist)
nmds_plot(dist, colorby=compiled.data.invasion$Factor)


###CONSTRAINED Ordination
head(species.invasion)
species.invasion.12 <- read.csv(file.choose())
row.names(species.invasion.12) <- species.invasion.12$ID
species.invasion.12 <- species.invasion.12[,-1]
head(species.invasion.12)

standardized.species.invasion.12<-decostand(species.invasion.12, method="total", MARGIN=2)
hellinger.species.invasion.12<-decostand(species.invasion.12, method="hellinger")
is.na(species.invasion)

capscale_plot.low<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2.low#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="invasion  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}

capscale_plot.control<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2.control#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="invasion  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}


?capscale
capscale_plot<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordiellipse(m,groups=new.compiled.12$CO2,draw="polygon",col= "grey90",label=F)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=2)
  legend("topright", title ="invasion  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.9)
}

m<-capscale(standardized.species.invasion ~ CO2*Invasives,compiled.data.invasion , distance="bray")
capscale_plot(m, colorby=compiled.data.invasion$Factor)

head(species.invasion.8to12)

species.invasion.8to12<-mean()


?capscale
?prc
library(vegan)

m1<-capscale(standardized.species.invasion ~ min.10.pH*Invasives,compiled.data.invasion.12 , distance="bray")
capscale_plot(m1, colorby=compiled.data.invasion.12$Factor)
adonis(standardized.species.invasion ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)



m2<-capscale(standardized.species.invasion.12 ~ hydrogen.concentration*Invasives,compiled.data.invasion.12 , distance="bray")
capscale_plot(m2, colorby=compiled.data.invasion.12$Factor)
adonis(standardized.species.invasion.12 ~ hydrogen.concentration*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)


m3<-capscale(hellinger.species.invasion.12 ~ hydrogen.concentration*Invasives,compiled.data.invasion.12 , distance="euclidean")
capscale_plot(m3, colorby=compiled.data.invasion.12$Factor)
adonis(standardized.species.invasion.12 ~ hydrogen.concentration*Invasives, method="euclidean", permutations = 9999, data=compiled.data.invasion.12)

?capscale
species.invasion8to12

species.invasion8to12HIGH.AIR<-aggdata[aggdata$Factor=="HighAmbient",]
species.invasion8to12HIGH.CO2<-aggdata[aggdata$Factor=="HighElevated",]
species.invasion8to12LOW.AIR<-aggdata[aggdata$Factor=="LowAmbient",]
species.invasion8to12LOW.CO2<-aggdata[aggdata$Factor=="LowElevated",]
species.invasion8to12NONE.AIR<-aggdata[aggdata$Factor=="AIR, Absent",]
species.invasion8to12NONE.CO2<-aggdata[aggdata$Factor=="NoneElevated",]


m.8to12<-capscale(species.invasion.agg ~ CO2*Invasives,compiled.data.invasion.12.agg, distance="bray")
capscale_plot(m.8to12, colorby=compiled.data.invasion$Factor)

library(vegan)


length(species.invasion.agg$alive.bot)
species.invasion.agg<-species.invasion.agg[,-11]

species.invasion.agg.hell<-decostand(species.invasion.agg, method="hellinger")

m.8to12.pH<-capscale(species.invasion.agg ~min.10.pH*Invasives,compiled.data.invasion.12, distance="bray")
capscale_plot(m.8to12.pH, colorby=compiled.data.invasion$Factor)
adonis(species.invasion.agg ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)


m.8to12.pH.hell<-capscale(species.invasion.agg.hell ~min.10.pH*Invasives,compiled.data.invasion.12, distance="euclidean")
capscale_plot(m.8to12.pH.hell, colorby=compiled.data.invasion$Factor)
adonis(species.invasion.agg.hell ~ min.10.pH*Invasives, method="euclidean", permutations = 9999, data=compiled.data.invasion.12)


m.8to12.CO2.hell<-capscale(species.invasion.agg.hell ~CO2*Invasives,compiled.data.invasion.12, distance="euclidean")
capscale_plot(m.8to12.CO2.hell, colorby=compiled.data.invasion$Factor)
adonis(species.invasion.agg.hell ~ CO2*Invasives, method="euclidean", permutations = 9999, data=compiled.data.invasion.12)





adonis(species.invasion.agg ~ agg.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)
adonis(species.invasion.agg ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)

standardized.species.invasion.agg<-decostand(species.invasion.agg, method="total", MARGIN=2)
standardized.species.invasion.agg.sessile<-decostand(species.invasion.agg.sessile, method="total", MARGIN=2)
2=aaq
bd<-betadisper(dist[[3]],groups)

m.8to12.CO2.stand<-capscale(standardized.species.invasion.agg ~CO2*Invasives,compiled.data.invasion.12, distance="bray")
capscale_plot(m.8to12.CO2.stand, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.agg  ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)


m.8to12.pH.stand<-capscale(standardized.species.invasion.agg ~min.10.pH*Invasives,compiled.data.invasion.12, distance="bray")
capscale_plot(m.8to12.pH.stand, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.agg  ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)



m.8to12.pH.stand.jac<-capscale(standardized.species.invasion.agg ~min.10.pH*Invasives,compiled.data.invasion.12, distance="jaccard")
capscale_plot(m.8to12.pH.stand.jac, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.agg  ~ min.10.pH*Invasives, method="jaccard", permutations = 9999, data=compiled.data.invasion.12)

m.8to12.pH.stand.sessile<-capscale(standardized.species.invasion.agg.sessile ~min.10.pH*Invasives,compiled.data.invasion.12, distance="jaccard")
capscale_plot(m.8to12.pH.stand.sessile, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.agg.sessile  ~ min.10.pH*Invasives, method="jaccard", permutations = 9999, data=compiled.data.invasion.12)


m.8to12.pH.stand.sessile<-capscale(standardized.species.invasion.agg.sessile ~min.10.pH*Invasives,compiled.data.invasion.12, distance="bray")
capscale_plot(m.8to12.pH.stand.sessile, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.agg.sessile  ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)






capscale_plot.control(m.8to12.pH.stand, colorby=compiled.data.invasion$Factor)


m.8to12.pH.stand.plot<-capscale_plot(m.8to12.pH.stand, colorby=compiled.data.invasion$Factor)

adonis(standardized.species.invasion.agg  ~ agg.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion.12)



dist.standardized.species.invasion.agg <- vegdist(standardized.species.invasion.agg, method = "bray")

mod.CO2<-betadisper(dist.standardized.species.invasion.agg, compiled.data.invasion.12$CO2, type="centroid")
anova(mod.CO2)

###
head(compiled.data.invasion.12)

### pH doesn't work b/c all diff groups... 
mod.pH<-betadisper(dist.standardized.species.invasion.agg, compiled.data.invasion.12$min.10.pH, type="centroid")
anova(mod.pH)

mod.invasion<-betadisper(dist.standardized.species.invasion.agg, compiled.data.invasion.12$Invasives, type="centroid")
bd<-anova(mod.invasion)

mod.invasion<-betadisper(dist.standardized.species.invasion.agg, compiled.data.invasion.12$Factor, type="centroid")
anova(mod.invasion)



mod.invasion<-betadisper(dist.standardized.species.invasion.agg, compiled.data.invasion.12$Factor, type="centroid")
anova(mod.invasion)
plot(mod.invasion)
bd<-betadisper(dist[[3]],groups)


dist.standardized.species.invasion.agg <- vegdist(standardized.species.invasion.agg, method = "bray")

# Now, we'll calculate the Jaccard index and its partitions of turnover and nestedness. We can calculate Sorensen index instead by using the argument     index.family="sorensen"    .
install.packages("betapart")
library(betapart)

presabs<-ifelse(standardized.species.invasion.agg>0,1,0)
dist.agg.stand.pair.jacc<-beta.pair(presabs, index.family="jaccard")
bd<-betadisper(dist.agg.stand.pair.jacc[[3]],compiled.data.invasion.12$Factor)
bd.nestedness<-betadisper(dist.agg.stand.pair.jacc[[2]],compiled.data.invasion.12$Factor)
bd.turnover<-betadisper(dist.agg.stand.pair.jacc[[1]],compiled.data.invasion.12$Factor)


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
dist.agg.stand.pair.bray<-bray.part(standardized.species.invasion.agg)
bd.bray<-betadisper(dist.agg.stand.pair.bray[[3]],compiled.data.invasion.12$Factor)
bd.nestedness.bray<-betadisper(dist.agg.stand.pair.bray[[2]],compiled.data.invasion.12$Factor)
bd.turnover.bray<-betadisper(dist.agg.stand.pair.bray[[1]],compiled.data.invasion.12$Factor)


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
dist.agg.stand.pair.bray.co2<-bray.part(standardized.species.invasion.agg)
bd.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[3]],compiled.data.invasion.12$CO2)
bd.nestedness.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[2]],compiled.data.invasion.12$CO2)
bd.turnover.bray.co2<-betadisper(dist.agg.stand.pair.bray.co2[[1]],compiled.data.invasion.12$CO2)


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
dist.agg.pair.bray.co2.unstand<-bray.part(species.invasion.agg)
bd.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[3]],compiled.data.invasion.12$CO2)
bd.nestedness.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[2]],compiled.data.invasion.12$CO2)
bd.turnover.bray.co2.unstand<-betadisper(dist.agg.pair.bray.co2.unstand[[1]],compiled.data.invasion.12$CO2)


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
dist.agg.stand.pair.bray.Invasives<-bray.part(standardized.species.invasion.agg)
bd.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[3]],compiled.data.invasion.12$Invasives)
bd.nestedness.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[2]],compiled.data.invasion.12$Invasives)
bd.turnover.bray.Invasives<-betadisper(dist.agg.stand.pair.bray.Invasives[[1]],compiled.data.invasion.12$Invasives)


plot(bd.bray.Invasives)
anova(bd.bray.Invasives)
boxplot(bd.bray.Invasives)

plot(bd.nestedness.bray.Invasives)
anova(bd.nestedness.bray.Invasives)
boxplot(bd.nestedness.bray.Invasives)

plot(bd.turnover.bray.Invasives)
anova(bd.turnover.bray.Invasives)
boxplot(bd.turnover.bray.Invasives)



##### not standardized invasion
dist.agg.pair.bray.Invasives.unstand<-bray.part(species.invasion.agg)
bd.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[3]],compiled.data.invasion.12$Invasives)
bd.nestedness.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[2]],compiled.data.invasion.12$Invasives)
bd.turnover.bray.Invasives.unstand<-betadisper(dist.agg.pair.bray.Invasives.unstand[[1]],compiled.data.invasion.12$Invasives)


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
###standardized species.invasion agg might be the way to go....? 


##### 

capscale_plot.none(m, colorby=compiled.data.invasion$Factor)
capscale_plot.low(m, colorby=compiled.data.invasion$Factor)
capscale_plot.high(m, colorby=compiled.data.invasion$Factor)

adonis(standardized.species.invasion ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)


adonis(standardized.species.invasion ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)
#can do by CO2 or by min.6.pH

adonis(standardized.species.invasion ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)


adonis(standardized.species.invasion ~ Invasives*min.6.pH, method="bray", permutations = 9999, data=compiled.data.invasion)



### Betadisperion

dist.standardized.6 <- vegdist(standardized.species.invasion, method = "bray")




mod.CO2<-betadisper(dist.standardized.6, compiled.data.invasion$CO2, type="centroid")
anova(mod.CO2)


mod.pH<-betadisper(dist.standardized.6, compiled.data.invasion$min.10.pH, type="centroid")
anova(mod.pH)

mod.invasion<-betadisper(dist.standardized.6, compiled.data.invasion$Invasives, type="centroid")
anova(mod.invasion)

mod.Combined<-betadisper(dist.standardized.6, compiled.data.invasion$Factor, type="centroid")
anova(mod.Combined)

pmod <- permutest(mod.Combined, permutations = 99, pairwise = TRUE)


(mod.HSD <- TukeyHSD(mod.Combined))
plot(mod.HSD)








######## Sessile only
head(standardized.species.invasion)
standardized.species.invasion.sessile<-standardized.species.invasion[,-(10:12)]
standardized.species.invasion.sessile<-standardized.species.invasion.sessile[,-18]
head(standardized.species.invasion.sessile)

m.sessile<-capscale(standardized.species.invasion.sessile ~ CO2*Invasives,compiled.data.invasion , distance="bray")
capscale_plot(m.sessile, colorby=compiled.data.invasion$Factor)
adonis(standardized.species.invasion.sessile ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)



###pa presence absence
m.pa<-capscale(standardized.species.invasion ~ CO2*Invasives,compiled.data.invasion , distance="jaccard")

capscale_plot(m.pa, colorby=compiled.data.invasion$Factor)



########Correlation analysis #################################################

species.invasion.sum<-apply(species.invasion, 2, sum)
species.invasion.sorted<-species.invasion[,order(species.invasion.sum, decreasing=TRUE)]
head(species.invasion.sorted)
species.invasion.small<-species.invasion.sorted[,1:15]
species.invasion.very.small<-species.invasion.sorted[,1:7]
species.invasion.top.ten<-species.invasion.sorted[,1:10]
species.invasion.top.five<-species.invasion.sorted[,1:5]


head(species.invasion.very.small)

species.invasion.small.std<-decostand(species.invasion.small, method="total", MARGIN=2)
species.invasion.small.t<-t(species.invasion.small.std)
head(species.invasion.small.t)

species.invasion.t.kmeans.casc<-cascadeKM(species.invasion.small.t, inf.gr=2, sup.gr=14, iter=100, criterion="calinski")
plot(species.invasion.t.kmeans.casc, sortg=TRUE)
species.invasion.t.kmeans.casc$results
species.invasion.t.kmeans.casc$partition

clusters<-species.invasion.t.kmeans.casc$partition[,1]
clusters

species.invasion.kendall.global<-kendall.global(species.invasion.small.hel, clusters)
species.invasion.kendall.global

species.invasion.kendall.post<-kendall.post(species.invasion.small.hel, clusters, nperm=9999)
species.invasion.kendall.post


################### Standardized and small 

species.invasion.small.t<-t(species.invasion.top.five)
head(species.invasion.small.t)

species.invasion.standardized.t<-decostand(species.invasion.small.t, method="total", MARGIN=1)

species.invasion.standardized.t.dist<-dist(species.invasion.standardized.t)

coldiss(species.invasion.standardized.t.dist, diag=TRUE)


############# exact same result if you standardize before or after so that's good. 

species.invasion.standardized.small<-decostand(species.invasion.very.small, method="total", MARGIN=2)

species.invasion.standardized.small.t<-t(species.invasion.standardized.small)

species.invasion.standardized.t.dist<-dist(species.invasion.standardized.small.t)

coldiss(species.invasion.standardized.t.dist, diag=TRUE)



dist.small <- vegdist(species.invasion.small.t, method = "bray")

coldiss(dist.small, diag=TRUE)


##################
head(species.invasion)

scaled.species.invasion<-scale(species.invasion)
fix(species.invasion)
head(scaled.species.invasion)

standardized.species.invasion<-decostand(species.invasion, method="total", MARGIN=2)



head(standardized.species.invasion.2)
wisconsin.species.invasion<- wisconsin(species.invasion)
head(wisconsin.species.invasion)
head(species.invasion)




#species.invasion.2<-species.invasion.2[, -11]

head(compiled.data.invasion)

adonis(species.invasion.2 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)

adonis(standardized.species.invasion.2 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)


adonis(standardized.species.invasion.3 ~ CO2*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)
adonis(standardized.species.invasion.3 ~ min.6.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)
########


adonis(standardized.species.invasion.2 ~ CO2*Invasives, method="jaccard", permutations = 9999, data=compiled.data.invasion)




adonis(species.invasion ~ min.6.pH, method="bray", permutations = 9999, data=compiled.data.invasion)

adonis(species.invasion ~ av.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.invasion)



dist.standardized.2 <- vegdist(standardized.species.invasion.2, method = "bray")


mod.CO2<-betadisper(dist.standardized.2, compiled.data.invasion$CO2, type="centroid")
anova(mod.CO2)

mod.invasion<-betadisper(dist.standardized.2, compiled.data.invasion$Invasives, type="centroid")
anova(mod.invasion)

mod.Combined<-betadisper(dist.standardized.2, compiled.data.invasion$Factor, type="centroid")
anova(mod.Combined)


######### SIMPER - type analysis?
#### simper is contribution to dissimilarity between sites
#### indv val is indicative of fidelity of that species.invasion to a particular site
#### 1 - mv abunda - 2- indval? 

##### do I do CO2 or invasion as the measure?? Too complicated to do both...?? 

sim <- with(compiled.data.invasion, simper(species.invasion,Invasives))
summary(sim)

simper(species.invasion,compiled.data.invasion$CO2, permutations = 0, trace = FALSE)


library(mvabund)
?decostand
tasmvabund <- mvabund(Tasmania$copepods)
plot(tasmvabund ~ treatment, col = as.numeric(block))



species.invasion.mvabund<-mvabund(standardized.species.invasion)
plot(species.invasion.mvabund ~ compiled.data.invasion$Factor)
#### not sure how useful this is. 

##### Indval
###"The indval approach looks for species.invasion that are
#both necessary and sufficient, i.e. if you find that species.invasion you should be in that type,
#and if you are in that type you should find that species.invasion
install.packages("labdsv")
library(labdsv)
iva<-indval(standardized.species.invasion,compiled.data.invasion$Factor)


iva$relfrq
iva$relabu
iva$indval

gr<- iva$maxcls[iva$pval<=0.05]
iv<- iva$indcls[iva$pval<=0.05]
pv<- iva$pval[iva$pval<=0.05]
fr<-apply(standardized.species.invasion>0, 2, sum)[iva$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg<-fidg[order(fidg$group,-fidg$indval), ]
fidg


isamic(standardized.species.invasion,compiled.data.invasion$CO2, sort=TRUE)



####


species.invasion.sum<-apply(species.invasion, 2, sum)
species.invasion.sorted<-species.invasion[,order(species.invasion.sum, decreasing=TRUE)]
species.invasion.3<-species.invasion.sorted[,1:22]
head(species.invasion.3)

standardized.species.invasion.3<-decostand(species.invasion.3, method="total", MARGIN=2)

iva<-indval(standardized.species.invasion.3,compiled.data.invasion$CO2)



#####
library(vegan)
set.seed(2)


treat=c(rep("High",19),rep("high",20), rep("Low", 20))
ordiplot(example_NMDS,type="n")
ordiellipse(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(example_NMDS,display="species.invasion",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)


groups = c(rep("group1", 36), rep("group2", 64))

orditorp(example_NMDS, display = "species.invasion")
