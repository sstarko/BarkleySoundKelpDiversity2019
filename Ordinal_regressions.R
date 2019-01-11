##Ordinal Logistic Regressions
library(MASS)
library(ordinal)
library(VGAM)

setwd("~/Desktop/Kelp Diversity Final/ordinal logistic regressions")
dat<-read.csv("Abundance_allyears_trimmed_44sites.csv")

#Make data ordinal
dat$Alaria_marginata<-ordered(dat$Alaria_marginata)
str(dat$Alaria_marginata)
dat$Costaria_costata<-ordered(dat$Costaria_costata)
dat$Egregia_menzieii<-ordered(dat$Egregia_menzieii)
dat$Macrocystis_pyrifera<-ordered(dat$Macrocystis_pyrifera)
dat$Nereocystis_luetkeana<-ordered(dat$Nereocystis_luetkeana)
dat$Saccharina_sessilis<-ordered(dat$Saccharina_sessilis)
dat$Ecklonia_arborea<-ordered(dat$Ecklonia_arborea)
dat$Laminaria_setchellii<-ordered(dat$Laminaria_setchellii)
dat$Lessoniopsis_littoralis<-ordered(dat$Lessoniopsis_littoralis)

dat$Fucus_spp<-ordered(dat$Fucus_spp)
dat$Sargassum_muticum<-ordered(dat$Sargassum_muticum)
dat$Phyllospadix_spp<-ordered(dat$Phyllospadix_spp)


##Run proportional odds models
anova(vglm(Costaria_costata~WA*Year2, family=propodds, data=dat))
anova(vglm(Alaria_marginata~WA*Year2, family=propodds, data=dat))
anova(vglm(Egregia_menzieii~WA*Year2, family=propodds, data=dat))
anova(vglm(Egregia_menzieii~WA*Year2, family=propodds, data=dat))

anova(vglm(Macrocystis_pyrifera~WA*Year2, family=propodds, data=dat))
anova(vglm(Nereocystis_luetkeana~WA*Year2, family=propodds, data=dat))
anova(vglm(Saccharina_sessilis~WA*Year2, family=propodds, data=dat))
anova(vglm(Ecklonia_arborea~WA*Year2, family=propodds, data=dat))
anova(vglm(Laminaria_setchellii~WA*Year, family=propodds, data=dat))
anova(vglm(Lessoniopsis_littoralis~WA*Year, family=propodds, data=dat))

anova(vglm(Egregia_menzieii~WA*Year2, family=propodds, data=dat_no1993))

dat_no1993<-dat[dat$Year>1993,]

dat1<-dat[dat$WA==1,]
dat3<-dat[dat$WA==3,]

#Non-kelps
anova(vglm(Fucus_spp~WA*Year2, family=propodds, data=dat))
summary(vglm(Fucus_spp~WA*Year2, family=propodds, data=dat3))
anova(vglm(Sargassum_muticum~WA*Year2, family=propodds, data=dat))
summary(vglm(Sargassum_muticum~Year2, family=propodds, data=dat1))
anova(vglm(Phyllospadix_spp~WA*Year2, family=propodds, data=dat))

####This reordering is for plotting only
dat$Alaria_marginata<-ordered(dat$Alaria_marginata, levels=c("1","2","3","0"))
str(dat$Alaria_marginata)
dat$Costaria_costata<-ordered(dat$Costaria_costata, levels=c("1","2","3","0"))
dat$Egregia_menzieii<-ordered(dat$Egregia_menzieii, levels=c("1","2","3","0"))
dat$Macrocystis_pyrifera<-ordered(dat$Macrocystis_pyrifera, levels=c("1","2","3","0"))
dat$Nereocystis_luetkeana<-ordered(dat$Nereocystis_luetkeana, levels=c("1","2","3","0"))
dat$Saccharina_sessilis<-ordered(dat$Saccharina_sessilis, levels=c("1","2","3","0"))
dat$Ecklonia_arborea<-ordered(dat$Ecklonia_arborea, levels=c("1","2","3","0"))
dat$Laminaria_setchellii<-ordered(dat$Laminaria_setchellii, levels=c("1","2","3","0"))
dat$Fucus_spp<-ordered(dat$Fucus_spp, levels=c("1","2","3","0"))
dat$Sargassum_muticum<-ordered(dat$Sargassum_muticum, levels=c("1","2","3","0"))
dat$Phyllospadix_spp<-ordered(dat$Phyllospadix_spp, levels=c("1","2","3","0"))
dat$Lessoniopsis_littoralis<-ordered(dat$Lessoniopsis_littoralis, levels=c("1","2","3","0"))
dat$Pterygophora_californica<-ordered(dat$Pterygophora_californica, levels=c("1","2","3","0"))
dat$Pleurophycus_gardneri<-ordered(dat$Pleurophycus_gardneri, levels=c("1","2","3","0"))


WA1<-dat[dat$WA==1,]
WA2<-dat[dat$WA==2,]
WA3<-dat[dat$WA==3,]

par(mfrow=c(1,3), mar=c(4,4,1,1))
plot(Sargassum_muticum~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA3, ylab="", las=1, yaxt="n")
plot(Sargassum_muticum~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="", las=1)
plot(Sargassum_muticum~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)

plot(Fucus_spp~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",las=1)
plot(Fucus_spp~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA2, ylab="",las=1)
plot(Fucus_spp~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)

plot(Phyllospadix_spp~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",las=1)
plot(Phyllospadix_spp~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA2, ylab="",las=1)
plot(Phyllospadix_spp~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)


par(mfrow=c(1,3), mar=c(4,4,1,1))
plot(Alaria_marginata~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="")
plot(Alaria_marginata~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Alaria_marginata~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA1, ylab="")

par(mfrow=c(1,3), mar=c(4,4,1,1))
plot(Costaria_costata~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="")
plot(Costaria_costata~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Costaria_costata~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA1, ylab="")

plot(Ecklonia_arborea~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA1, ylab="")
plot(Ecklonia_arborea~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA2, ylab="")
plot(Ecklonia_arborea~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA3, ylab="")

plot(Egregia_menzieii~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA1, ylab="")
plot(Egregia_menzieii~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA2, ylab="")
plot(Egregia_menzieii~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA3, ylab="")


plot(Laminaria_setchellii~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA1, ylab="")
plot(Laminaria_setchellii~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Laminaria_setchellii~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA3, ylab="")


plot(Macrocystis_pyrifera~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",yaxt="n")
plot(Macrocystis_pyrifera~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Macrocystis_pyrifera~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA1, ylab="")

plot(Nereocystis_luetkeana~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",yaxt="n")
plot(Nereocystis_luetkeana~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Nereocystis_luetkeana~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA1, ylab="")


plot(Saccharina_sessilis~as.factor(Year2),col=c( "lightblue", "blue","navy", "white"), data=WA1, ylab="")
plot(Saccharina_sessilis~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA2, ylab="")
plot(Saccharina_sessilis~as.factor(Year2),col=c( "lightblue", "blue","navy","white"), data=WA3, ylab="")



plot(Lessoniopsis_littoralis~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",las=1)
plot(Lessoniopsis_littoralis~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA2, ylab="",las=1)
plot(Lessoniopsis_littoralis~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)

plot(Pleurophycus_gardneri~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",las=1)
plot(Pleurophycus_gardneri~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA2, ylab="",las=1)
plot(Pleurophycus_gardneri~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)

plot(Pterygophora_californica~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA3, ylab="",las=1)
plot(Pterygophora_californica~as.factor(Year2),col=c("lightblue", "blue","navy", "white"), data=WA2, ylab="",las=1)
plot(Pterygophora_californica~as.factor(Year2),col=c("lightblue", "blue","navy","white"), data=WA1, ylab="",las=1)

