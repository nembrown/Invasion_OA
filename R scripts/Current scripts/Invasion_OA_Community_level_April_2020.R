#read in packages
library(vegan)
library(ggplot2)
library(betapart)
library(bipartite)
library(car)
library(fitdistrplus)

## read in data from both mesocosms and tiles
  

invasion.exp.data<-read.csv("C:Biological data/allweeks_cover_counts_with_pres.csv",stringsAsFactors = FALSE, na.strings = c("NA","") )
head(invasion.exp.data)
#ordered and unordered factors
invasion.exp.data$oTreatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=TRUE)
invasion.exp.data$Treatment<-factor(invasion.exp.data$Treatment, levels=c("AIRAbsent",  "CO2Absent", "CO2Present", "AIRPresent"), ordered=FALSE)

invasion.exp.data$oInvasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=TRUE)
invasion.exp.data$Invasives<-factor(invasion.exp.data$Invasives, levels=c("Absent", "Present"), ordered=FALSE)

invasion.exp.data$oCO2.Treatment<-factor(invasion.exp.data$CO2.Treatment, levels=c("AIR", "CO2"), ordered=TRUE)
invasion.exp.data$CO2.Treatment<-factor(invasion.exp.data$CO2.Treatment, levels=c("AIR", "CO2"), ordered=FALSE)

#New variables
# invasion.exp.data$num.nudi<-invasion.exp.data$nudibranch+invasion.exp.data$nudi.eggs+invasion.exp.data$nudi.hatched

invasion.exp.data$bot.total<-invasion.exp.data$botryllid + invasion.exp.data$bot.eaten
invasion.exp.data$mem.total<-invasion.exp.data$membranipora + invasion.exp.data$mem.eaten + invasion.exp.data$mem.dead
invasion.exp.data$corella.total<-invasion.exp.data$corella + invasion.exp.data$dead.corella

invasion.exp.data$prop.mem.dead<-invasion.exp.data$mem.dead/invasion.exp.data$membranipora
invasion.exp.data$prop.mem.eaten<-invasion.exp.data$mem.eaten/invasion.exp.data$membranipora


invasion.exp.data$hydroid.001<-(0.01*(invasion.exp.data$hydroid))+0.01
invasion.exp.data$botryllid.001<-(0.01*(invasion.exp.data$botryllid))+0.01
invasion.exp.data$bot.eaten.001<-(0.01*(invasion.exp.data$bot.eaten))+0.01
invasion.exp.data$mem.eaten.001<-(0.01*(invasion.exp.data$mem.eaten))+0.01
invasion.exp.data$mem.total.001<-(0.01*(invasion.exp.data$mem.total))+0.01

invasion.exp.data$membranipora.001<-(0.01*(invasion.exp.data$membranipora))+0.01
invasion.exp.data$mussel.001<-(0.01*(invasion.exp.data$mussel))+0.01
invasion.exp.data$didemnum<-invasion.exp.data$white.bryo+invasion.exp.data$fan.bryo
invasion.exp.data$folliculina<-invasion.exp.data$protozoa
invasion.exp.data$folliculina.001<-(0.01*(invasion.exp.data$folliculina))+0.01
invasion.exp.data$didemnum.001<-(0.01*(invasion.exp.data$didemnum))+0.01
invasion.exp.data$occupied.space<-(100 - invasion.exp.data$bare)
invasion.exp.data$occupied.space.001<-(0.01*(invasion.exp.data$occupied.space))+0.01
invasion.exp.data$native.occupied.space<-(100 - invasion.exp.data$botryllid  -invasion.exp.data$bot.eaten -  invasion.exp.data$bare)
invasion.exp.data$native.occupied.space.001<-(0.01*(invasion.exp.data$native.occupied.space))+0.01



head(invasion.exp.data)
#splitting into mid and end
invasion.exp.data.16<-invasion.exp.data %>% filter(Week==16)
invasion.exp.data.8<-invasion.exp.data %>% filter(Week==8)


 #other bryo = cribrillina 
    
#Select data to be used: 
# some cover and some counts: 

#counts & cover
native_sp_names<-c("Tile.ID", 
                          "num.nudi.all" ,
                          "mussel" ,"hydroid","mem.total", "folliculina",
                          "num.corella" ,
                          "num.red.bryo" ,
                          "num.white.bryo" ,
                          "num.serpulid" ,
                          "num.barn" ,
                          "other.bryo" ,
                          "clam", "oyster", "bubble snail", "bubble eggs", "slime", "chiton", "yellow.sponge")
                          

invasion.exp.data.8.tile.selected<-invasion.exp.data.8[,colnames(invasion.exp.data.8) %in% native_sp_names]
head(invasion.exp.data.8.tile.selected)


#community level data combined counts and percent cover
species.rec_cover.8 <-invasion.exp.data.8.tile.selected

just.species.rec_cover.8<- species.rec_cover.8[,-1]
head(just.species.rec_cover.8)
#richness can be from count data - it's just pres/abs of given species
invasion.exp.data.8$num.species.no.bot <- specnumber(just.species.rec_cover.8)



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity
names_cover_food_exp_tile<-c("Tile.ID",
                             "disporella",
                             "mussel" ,"hydroid","membranipora","mem.dead", "mem.eaten", "folliculina",
                             "corella" ,"dead.corella",
                             "red.bryo" ,
                             "serpulid" ,
                             "barn" ,
                             "other.bryo" ,
                             "clam", "oyster", "disporella", "nudi.cover", "nudi eggs", "bubble snail", "bubble eggs", "slime", "chiton", "yellow.sponge")



species.cover.8 <- invasion.exp.data.8[,colnames(invasion.exp.data.8) %in% names_cover_food_exp_tile]
head(species.cover.8)
just.species.cover.8<-species.cover[,-1]

species.cover.8$shannon.diversity<-diversity(just.species.cover.8, index="shannon")
species.cover.8$evenness<-species.cover.8$shannon.diversity/(log(invasion.exp.data.8$num.species.no.bot))

#evenness has to be created from just cover data
invasion.exp.data.8$evenness<-species.cover.8$evenness



head(  invasion.exp.data.8)
# MDS ---------------------------------------------------------------------

compiled.data.8 <- invasion.exp.data.8 %>% dplyr::select(Mesocosm, Tile.ID, Invasives, CO2.Treatment, Treatment, pH, av.pH.all, pH.uptowk, min.10.pH, Week, bot.total)
head(compiled.data.8)

row.names(compiled.data.8)<-compiled.data.8$Tile.ID

head(compiled.data.8)
compiled.data.8$Treatment<-as.factor(compiled.data.8$Treatment) 

#Combinging species and environment
all.data.rec_cover.8<-merge(species.rec_cover.8,compiled.data.8)
head(all.data.rec_cover.8)

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


summary(model.meso.bray.inv.8)
model.meso.bray.inv.8<-capscale(standardized.species.rec_cover.8 ~ min.10.pH*Invasives,compiled.data.8_zscores , distance="bray")
capscale_plot(model.meso.bray.inv.8, colorby=compiled.data.8$Treatment)
model.meso.bray.scores.8.inv<- as.data.frame(scores(model.meso.bray.inv.8)$sites)
model.meso.bray.scores.8.inv$Tile.ID<-species.rec_cover.8$Tile.ID
model.meso.bray.scores.CAP.inv.8<-merge(model.meso.bray.scores.8.inv, compiled.data.8, by="Tile.ID")
invasion.exp.data.8.community.inv<-model.meso.bray.scores.8.inv[,1:3]

model.meso.bray.8<-capscale(standardized.species.rec_cover.8 ~ bot.total*CO2.Treatment,compiled.data.8_zscores , distance="bray")
capscale_plot(model.meso.bray.8, colorby=compiled.data.8$Treatment)
model.meso.bray.scores.8<- as.data.frame(scores(model.meso.bray.8)$sites)
model.meso.bray.scores.8$Tile.ID<-species.rec_cover.8$Tile.ID
model.meso.bray.scores.CAP.inv.8<-merge(model.meso.bray.scores.8, compiled.data.8, by="Tile.ID")
invasion.exp.data.8.community<-model.meso.bray.scores.8[,1:3]


# betadispersion partitioned ----------------------------------------------


dist.part.bray.8<-bray.part(standardized.species.rec_cover.8)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
bd.bray.8<-betadisper(dist.part.bray.8[[3]],compiled.data.8_zscores$Treatment, type="centroid" )
bd.nestedness.bray.8<-betadisper(dist.part.bray.8[[2]],compiled.data.8_zscores$Treatment, type="centroid")
bd.turnover.bray.8<-betadisper(dist.part.bray.8[[1]],compiled.data.8_zscores$Treatment, type="centroid")

bd.overall.bray.distances.8<- as.data.frame(bd.bray.8$distances)
bd.overall.bray.distances.8$distcentroid<-bd.overall.bray.distances.8$`bd.bray.8$distances`
bd.overall.bray.distances.8$Tile.ID<-species.rec_cover.8$Tile.ID



bd.overall.bray.distances.2.8<-merge(bd.overall.bray.distances.8, compiled.data.8, by="Tile.ID")
head(bd.overall.bray.distances.2.8)

invasion.exp.data.8.community<-merge(bd.overall.bray.distances.8,invasion.exp.data.8.community, by="Tile.ID")
head(invasion.exp.data.8.community)

invasion.exp.data.8.community.inv$CAP1.inv<-invasion.exp.data.8.community.inv$CAP1

invasion.exp.data.8.community<-merge(invasion.exp.data.8.community, invasion.exp.data.8.community.inv[,3:4])
write.csv(invasion.exp.data.8.community,"C:Biological data/invasion.exp.data.8.community.csv", row.names=FALSE)


# Week 16 -----------------------------------------------------------------



#other bryo = cribrillina 

#Select data to be used: 
# some cover and some counts: 

#counts & cover
native_sp_names<-c("Tile.ID", 
                   "num.nudi.all" ,
                   "mussel" ,"hydroid","mem.total", "folliculina",
                   "num.corella" ,
                   "num.red.bryo" ,
                   "num.white.bryo" ,
                   "num.serpulid" ,
                   "num.barn" ,
                   "other.bryo" ,
                   "clam", "oyster", "bubble snail", "bubble eggs", "slime", "chiton", "yellow.sponge")


invasion.exp.data.16.tile.selected<-invasion.exp.data.16[,colnames(invasion.exp.data.16) %in% native_sp_names]
head(invasion.exp.data.16.tile.selected)


#community level data combined counts and percent cover
species.rec_cover.16 <-invasion.exp.data.16.tile.selected

just.species.rec_cover.16<- species.rec_cover.16[,-1]
head(just.species.rec_cover.16)
#richness can be from count data - it's just pres/abs of given species
invasion.exp.data.16$num.species.no.bot <- specnumber(just.species.rec_cover.16)



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity
names_cover_food_exp_tile<-c("Tile.ID",
                             "disporella",
                             "mussel" ,"hydroid","membranipora","mem.dead", "mem.eaten", "folliculina",
                             "corella" ,"dead.corella",
                             "red.bryo" ,
                             "serpulid" ,
                             "barn" ,
                             "other.bryo" ,
                             "clam", "oyster", "disporella", "nudi.cover", "nudi eggs", "bubble snail", "bubble eggs", "slime", "chiton", "yellow.sponge")



species.cover.16 <- invasion.exp.data.16[,colnames(invasion.exp.data.16) %in% names_cover_food_exp_tile]
head(species.cover.16)
just.species.cover.16<-species.cover[,-1]

species.cover.16$shannon.diversity<-diversity(just.species.cover.16, index="shannon")
species.cover.16$evenness<-species.cover.16$shannon.diversity/(log(invasion.exp.data.16$num.species.no.bot))

#evenness has to be created from just cover data
invasion.exp.data.16$evenness<-species.cover.16$evenness



head(  invasion.exp.data.16)
# MDS ---------------------------------------------------------------------

compiled.data.16 <- invasion.exp.data.16 %>% dplyr::select(Mesocosm, Tile.ID, Invasives, CO2.Treatment, Treatment, pH, av.pH.all, pH.uptowk, min.10.pH, Week, bot.total)
head(compiled.data.16)

row.names(compiled.data.16)<-compiled.data.16$Tile.ID

head(compiled.data.16)
compiled.data.16$Treatment<-as.factor(compiled.data.16$Treatment) 

#Combinging species and environment
all.data.rec_cover.16<-merge(species.rec_cover.16,compiled.data.16)
head(all.data.rec_cover.16)

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

model.meso.bray.inv.16<-capscale(standardized.species.rec_cover.16 ~ min.10.pH*Invasives,compiled.data.16_zscores , distance="bray")
capscale_plot(model.meso.bray.inv.16, colorby=compiled.data.16$Treatment)
adonis(standardized.species.rec_cover.16 ~ min.10.pH*Invasives,compiled.data.16_zscores , distance="bray")
model.meso.bray.scores.16.inv<- as.data.frame(scores(model.meso.bray.inv.16)$sites)
model.meso.bray.scores.16.inv$Tile.ID<-species.rec_cover.16$Tile.ID
invasion.exp.data.16.community.inv<-model.meso.bray.scores.16.inv[,1:3]

model.meso.bray.16<-capscale(standardized.species.rec_cover.16 ~ bot.total*CO2.Treatment,compiled.data.16_zscores , distance="bray")
capscale_plot(model.meso.bray.16, colorby=compiled.data.16$Treatment)
model.meso.bray.scores.16<- as.data.frame(scores(model.meso.bray.16)$sites)
model.meso.bray.scores.16$Tile.ID<-species.rec_cover.16$Tile.ID
invasion.exp.data.16.community<-model.meso.bray.scores.16[,1:3]

summary(model.meso.bray.inv.16)

# betadispersion partitioned ----------------------------------------------


dist.part.bray.16<-bray.part(standardized.species.rec_cover.16)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
bd.bray.16<-betadisper(dist.part.bray.16[[3]],compiled.data.16_zscores$Treatment, type="centroid" )
bd.nestedness.bray.16<-betadisper(dist.part.bray.16[[2]],compiled.data.16_zscores$Treatment, type="centroid")
bd.turnover.bray.16<-betadisper(dist.part.bray.16[[1]],compiled.data.16_zscores$Treatment, type="centroid")

bd.overall.bray.distances.16<- as.data.frame(bd.bray.16$distances)
bd.overall.bray.distances.16$distcentroid<-bd.overall.bray.distances.16$`bd.bray.16$distances`
bd.overall.bray.distances.16$Tile.ID<-species.rec_cover.16$Tile.ID



bd.overall.bray.distances.2.16<-merge(bd.overall.bray.distances.16, compiled.data.16, by="Tile.ID")
head(bd.overall.bray.distances.2.16)

invasion.exp.data.16.community<-merge(bd.overall.bray.distances.16,invasion.exp.data.16.community, by="Tile.ID")
head(invasion.exp.data.16.community)


invasion.exp.data.16.community.inv$CAP1.inv<-invasion.exp.data.16.community.inv$CAP1

invasion.exp.data.16.community<-merge(invasion.exp.data.16.community, invasion.exp.data.16.community.inv[,3:4])
write.csv(invasion.exp.data.16.community,"C:Biological data/invasion.exp.data.16.community.csv", row.names=FALSE)

