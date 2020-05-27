#read in packages
library(vegan)
library(ggplot2)
library(betapart)
library(bipartite)
library(car)
library(fitdistrplus)

## read in data from both mesocosms and tiles
  
  
invasion.exp.data<-read.csv("C:Biological data/allweeks_cover_counts_without_pres.csv",stringsAsFactors = FALSE, na.strings = c("NA","") )
head(invasion.exp.data)
invasion.exp.data$num.nudi<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.eggs+invasion.exp.data$nudi.hatched
invasion.exp.data$nudi.cover<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.hatched

invasion.exp.data$disporella<-invasion.exp.data$white.bryo+invasion.exp.data$fan.bryo+invasion.exp.data$erect.bryo

invasion.exp.data$hydroid.001<-(0.01*(invasion.exp.data$hydroid))+0.01
invasion.exp.data$botryllid.001<-(0.01*(invasion.exp.data$botryllid))+0.01
invasion.exp.data$membranipora.001<-(0.01*(invasion.exp.data$membranipora))+0.01
invasion.exp.data$mussel.001<-(0.01*(invasion.exp.data$mussel))+0.01
invasion.exp.data$didemnum<-invasion.exp.data$white.bryo
invasion.exp.data$num.red.bryo<-invasion.exp.data$red.bryo
invasion.exp.data$folliculina<-invasion.exp.data$protozoa
invasion.exp.data$folliculina.001<-(0.01*(invasion.exp.data$folliculina))+0.01
invasion.exp.data$didemnum.001<-(0.01*(invasion.exp.data$didemnum))+0.01
invasion.exp.data$occupied.space<-(100 - invasion.exp.data$bare)
invasion.exp.data$occupied.space.001<-(0.01*(invasion.exp.data$occupied.space))+0.01
invasion.exp.data$native.occupied.space<-(100 - invasion.exp.data$botryllid  -invasion.exp.data$bot.eaten -  invasion.exp.data$bare)
invasion.exp.data$native.occupied.space.001<-(0.01*(invasion.exp.data$native.occupied.space))+0.01
invasion.exp.data.16<-invasion.exp.data %>% filter(Week==16)
invasion.exp.data.8<-invasion.exp.data %>% filter(Week==8)
  
 #other bryo = cribrillina 
    
#Select data to be used: 
# some cover and some counts: 

#counts
native_sp_names<-c("Tile.ID", 
                          "num.nudi" ,
                          "mussel" ,"hydroid","membranipora", "mem.eaten", "folliculina",
                          "num.corella" ,
                          "num.red.bryo" ,
                          "num.white.bryo" ,
                          "num.serpulid" ,
                          "num.barn" ,
                          "other.bryo" ,
                          "clam")
                          

invasion.exp.data.8.tile.selected<-invasion.exp.data.8[,colnames(invasion.exp.data.8) %in% native_sp_names]
head(invasion.exp.data.8.tile.selected)


#community level data combined counts and percent cover
species.rec_cover.8 <-invasion.exp.data.8.tile.selected

just.species.rec_cover.8<- species.rec_cover.8[,-1]
head(just.species.rec_cover.8)



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity

names_cover_food_exp_tile<-c("Tile.ID",
                             "disporella",
                             "mussel" ,"hydroid","membranipora","mem.dead", "mem.eaten", "folliculina",
                             "corella" ,"dead.corella",
                             "red.bryo" ,
                             "serpulid" ,
                             "barn" ,
                             "other.bryo" ,
                             "clam", "oyster", "disporella", "nudi.cover", "nudi eggs", "bubble snail", "bubble eggs", "slime", "chiton")



species.cover <- invasion.exp.data.8[,colnames(invasion.exp.data.8) %in% names_cover_food_exp_tile]
head(species.cover)
just.species.cover<-species.cover[,-1]

species.cover$richness<-specnumber(just.species.cover)
species.cover$shannon.diversity<-diversity(just.species.cover, index="shannon")
species.cover$evenness<-species.cover$shannon.diversity/(log(invasion.exp.data.8$num.species.no.bot))

#evenness has to be created from just cover data
invasion.exp.data.8$evenness<-species.cover$evenness

#richness can be from count data - it's just pres/abs of given species
invasion.exp.data.8$richness <- specnumber(just.species.rec_cover.8)

head(  invasion.exp.data.8)
# MDS ---------------------------------------------------------------------

compiled.data.8 <- invasion.exp.data.8 %>% dplyr::select(Mesocosm, Tile.ID, Invasives, CO2.Treatment, Treatment, pH, av.pH.all, pH.uptowk, min.10.pH, Week)
head(compiled.data.8)

row.names(compiled.data.8)<-compiled.data.8$Tile.ID

head(compiled.data.8)
compiled.data.8$Treatment<-as.factor(compiled.data.8$Treatment) 

#Combinging species and environment
all.data.rec_cover<-merge(species.rec_cover.8,compiled.data.8)
head(all.data.rec_cover)

colorset_invasives = c("Present"="#A20226" ,"Absent"="#818392")
theme_set(theme_classic(base_size = 6))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

colorset_treatment<-c("AIRPresent"="#A20226" ,"AIRAbsent"="#818392", "CO2.TreatmentPresent"="#A20226" ,"CO2.TreatmentAbsent"="#818392")

###CONSTRAINED Ordination

capscale_plot<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- colorset_treatment#vector of colors needed
  shapesies<-c( 16,16,2,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordisurf(m ~ min.10.pH, data=compiled.data.8_zscores, method = "REML", select = TRUE)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Invasives CO2.Treatment", legend=levels(colorby), col=cols, pch = shapesies, cex=.5)
}


# need to have zscores for pH ... otherwise evaluating at ph=0 which is not meaningful 
compiled.data.8_zscores<-compiled.data.8
compiled.data.8_zscores$min.10.pH<-scale(compiled.data.8$min.10.pH, center=TRUE, scale=TRUE)
head(compiled.data.8_zscores)



# Bray-Curtis Capscale / constrained ordination -------------------------------------------------------------


#Standardizing by total of the species either percent or count
# This makes each species on their own scale, so the Tile.ID got x % of the total mussels for eg.
standardized.species.rec_cover.8<-decostand(just.species.rec_cover.8, method="total", MARGIN=2)
head(standardized.species.rec_cover.8)

model.meso.bray<-capscale(standardized.species.rec_cover.8 ~ min.10.pH*Invasives,compiled.data.8_zscores , distance="bray")
capscale_plot(model.meso.bray, colorby=compiled.data.8$Treatment)

model.meso.bray.sf<-ordisurf(model.meso.bray ~ min.10.pH, data=compiled.data.8_zscores, method = "REML", select = TRUE)
summary(model.meso.bray.sf)

adonis(standardized.species.rec_cover.8 ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.8_zscores)
summary(model.meso.bray)

model.meso.bray.scores<- as.data.frame(scores(model.meso.bray)$sites)
head(model.meso.bray.scores)
model.meso.bray.scores$Tile.ID<-species.rec_cover.8$Tile.ID

model.meso.bray.scores.CAP<-merge(model.meso.bray.scores, compiled.data.8, by="Tile.ID")
head(model.meso.bray.scores.CAP)

write.csv(model.meso.bray.scores,"C:Biological data//model.meso.bray.scores.csv", row.names=FALSE)

invasion.exp.data.8.community<-model.meso.bray.scores[,1:3]

head(invasion.exp.data.8.community)

# betadispersion partitioned ----------------------------------------------


dist.part.bray<-bray.part(standardized.species.rec_cover.8)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
bd.bray<-betadisper(dist.part.bray[[3]],compiled.data.8_zscores$Treatment, type="centroid" )
bd.nestedness.bray<-betadisper(dist.part.bray[[2]],compiled.data.8_zscores$Treatment, type="centroid")
bd.turnover.bray<-betadisper(dist.part.bray[[1]],compiled.data.8_zscores$Treatment, type="centroid")

head(bd.bray)
plot(bd.bray, hull=FALSE, ellipse = TRUE)
anova(bd.bray)
boxplot(bd.bray)
permutest(bd.bray)

plot(bd.nestedness.bray)
anova(bd.nestedness.bray)
boxplot(bd.nestedness.bray)

plot(bd.turnover.bray)
anova(bd.turnover.bray)
boxplot(bd.turnover.bray)


bd.overall.bray.distances<- as.data.frame(bd.bray$distances)
head(bd.overall.bray.distances)
bd.overall.bray.distances$distcentroid<-bd.overall.bray.distances$`bd.bray$distances`
bd.overall.bray.distances$Tile.ID<-species.rec_cover.8$Tile.ID



bd.overall.bray.distances.2<-merge(bd.overall.bray.distances, compiled.data.8, by="Tile.ID")
head(bd.overall.bray.distances.2)

plot.overall.distcentroid.12.hydrogen<- ggplot(bd.overall.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives), shape=CO2.Treatment)) + guides(fill=FALSE) + scale_fill_manual(values=colorset_invasives)+ geom_smooth(aes(fill=Invasives), method="lm") +scale_shape_manual(values=c(19,17))
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Distance to centroid")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))
plot.overall.distcentroid.12.hydrogen
write.csv(bd.overall.bray.distances.2,"C:Data//Tile.ID inventory data/bd.overall.bray.distances.2.csv",row.names=FALSE )

invasion.exp.data.8.community<-merge(bd.overall.bray.distances,invasion.exp.data.8.community, by="Tile.ID")
head(invasion.exp.data.8.community)

write.csv(invasion.exp.data.8.community,"C:Biological data/invasion.exp.data.8.community.csv", row.names=FALSE)




############## Week 16
invasion.exp.data.16.tile.selected<-invasion.exp.data.16[,colnames(invasion.exp.data.16) %in% native_sp_names]
head(invasion.exp.data.16.tile.selected)


#community level data combined counts and percent cover
species.rec_cover.16 <-invasion.exp.data.16.tile.selected

just.species.rec_cover.16<- species.rec_cover.16[,-1]
head(just.species.rec_cover.16)



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity

names_cover_food_exp_tile<-c("Tile.ID",
                             "num.nudi" ,"disporella",
                             "mussel" ,"hydroid","membranipora", "mem.eaten", "folliculina",
                             "corella" ,"dead.corella",
                             "red.bryo" ,
                             "nserpulid" ,
                             "barn" ,
                             "other.bryo" ,
                             "clam" )



species.cover <- invasion.exp.data.16[,colnames(invasion.exp.data.16) %in% names_cover_food_exp_tile]
head(species.cover)
just.species.cover<-species.cover[,-1]

species.cover$richness<-specnumber(just.species.cover)
species.cover$shannon.diversity<-diversity(just.species.cover, index="shannon")
species.cover$evenness<-species.cover$shannon.diversity/(log(invasion.exp.data.16$num.species.no.bot))

#evenness has to be created from just cover data
invasion.exp.data.16$evenness<-species.cover$evenness

#richness can be from count data - it's just pres/abs of given species
invasion.exp.data.16$richness <- specnumber(just.species.rec_cover.16)

head(  invasion.exp.data.16)
# MDS ---------------------------------------------------------------------

compiled.data.16 <- invasion.exp.data.16 %>% dplyr::select(Mesocosm, Tile.ID, Invasives, CO2.Treatment, Treatment, pH, av.pH.all, pH.uptowk, min.10.pH, Week)
head(compiled.data.16)

row.names(compiled.data.16)<-compiled.data.16$Tile.ID

head(compiled.data.16)
compiled.data.16$Treatment<-as.factor(compiled.data.16$Treatment) 

#Combinging species and environment
all.data.rec_cover<-merge(species.rec_cover.16,compiled.data.16)
head(all.data.rec_cover)

colorset_invasives = c("Present"="#A20226" ,"Absent"="#16116392")
theme_set(theme_classic(base_size = 6))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))

colorset_treatment<-c("AIRPresent"="#A20226" ,"AIRAbsent"="#16116392", "CO2.TreatmentPresent"="#A20226" ,"CO2.TreatmentAbsent"="#16116392")

###CONSTRAINED Ordination

capscale_plot<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- colorset_treatment#vector of colors needed
  shapesies<-c( 16,16,2,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordisurf(m ~ min.10.pH, data=compiled.data.16_zscores, method = "REML", select = TRUE)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Invasives CO2.Treatment", legend=levels(colorby), col=cols, pch = shapesies, cex=.5)
}


# need to have zscores for pH ... otherwise evaluating at ph=0 which is not meaningful 
compiled.data.16_zscores<-compiled.data.16
compiled.data.16_zscores$min.10.pH<-scale(compiled.data.16$min.10.pH, center=TRUE, scale=TRUE)
head(compiled.data.16_zscores)



# Bray-Curtis Capscale / constrained ordination -------------------------------------------------------------


#Standardizing by total of the species either percent or count
# This makes each species on their own scale, so the Tile.ID got x % of the total mussels for eg.
standardized.species.rec_cover.16<-decostand(just.species.rec_cover.16, method="total", MARGIN=2)
head(standardized.species.rec_cover.16)

model.meso.bray<-capscale(standardized.species.rec_cover.16 ~ min.10.pH*Invasives,compiled.data.16_zscores , distance="bray")
capscale_plot(model.meso.bray, colorby=compiled.data.16$Treatment)

model.meso.bray.sf<-ordisurf(model.meso.bray ~ min.10.pH, data=compiled.data.16_zscores, method = "REML", select = TRUE)
summary(model.meso.bray.sf)

adonis(standardized.species.rec_cover.16 ~ min.10.pH*Invasives, method="bray", permutations = 9999, data=compiled.data.16_zscores)
summary(model.meso.bray)

model.meso.bray.scores<- as.data.frame(scores(model.meso.bray)$sites)
head(model.meso.bray.scores)
model.meso.bray.scores$Tile.ID<-species.rec_cover.16$Tile.ID

model.meso.bray.scores.CAP<-merge(model.meso.bray.scores, compiled.data.16, by="Tile.ID")
head(model.meso.bray.scores.CAP)

write.csv(model.meso.bray.scores,"C:Biological data//model.meso.bray.scores.csv", row.names=FALSE)

invasion.exp.data.16.community<-model.meso.bray.scores[,1:3]

head(invasion.exp.data.16.community)

# betadispersion partitioned ----------------------------------------------


dist.part.bray<-bray.part(standardized.species.rec_cover.16)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
bd.bray<-betadisper(dist.part.bray[[3]],compiled.data.16_zscores$Treatment, type="centroid" )
bd.nestedness.bray<-betadisper(dist.part.bray[[2]],compiled.data.16_zscores$Treatment, type="centroid")
bd.turnover.bray<-betadisper(dist.part.bray[[1]],compiled.data.16_zscores$Treatment, type="centroid")

head(bd.bray)
plot(bd.bray, hull=FALSE, ellipse = TRUE)
anova(bd.bray)
boxplot(bd.bray)
permutest(bd.bray)

plot(bd.nestedness.bray)
anova(bd.nestedness.bray)
boxplot(bd.nestedness.bray)

plot(bd.turnover.bray)
anova(bd.turnover.bray)
boxplot(bd.turnover.bray)


bd.overall.bray.distances<- as.data.frame(bd.bray$distances)
head(bd.overall.bray.distances)
bd.overall.bray.distances$distcentroid<-bd.overall.bray.distances$`bd.bray$distances`
bd.overall.bray.distances$Tile.ID<-species.rec_cover.16$Tile.ID



bd.overall.bray.distances.2<-merge(bd.overall.bray.distances, compiled.data.16, by="Tile.ID")
head(bd.overall.bray.distances.2)

plot.overall.distcentroid.12.hydrogen<- ggplot(bd.overall.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Invasives)) + geom_point(size=5,aes(colour=factor(Invasives), shape=CO2.Treatment)) + guides(fill=FALSE) + scale_fill_manual(values=colorset_invasives)+ geom_smooth(aes(fill=Invasives), method="lm") +scale_shape_manual(values=c(19,17))
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Distance to centroid")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))
plot.overall.distcentroid.12.hydrogen
write.csv(bd.overall.bray.distances.2,"C:Data//Tile.ID inventory data/bd.overall.bray.distances.2.csv",row.names=FALSE )

invasion.exp.data.16.community<-merge(bd.overall.bray.distances,invasion.exp.data.16.community, by="Tile.ID")
head(invasion.exp.data.16.community)

write.csv(invasion.exp.data.16.community,"C:Biological data/invasion.exp.data.16.community.csv", row.names=FALSE)

