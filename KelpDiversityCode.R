setwd("~/Path_to_data")
kelp.richness<-read.csv("KelpRichness2.csv")
ab<-read.csv("KelpAbundance.csv")
library(bayou)
library(onewaytests)
library(dunn.test)

add.error.bars <- function(X,Y,SE,w,col=1){
  X0 = X; Y0 = (Y-SE); X1 =X; Y1 = (Y+SE);
  arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col);
}

#Figure 3; 2017 and 2018 data are averaged, 1993-1994 data are averaged
par(mfrow=c(1,2))
boxplot(kelp.richness$HistoricRichness2[kelp.richness$WA2==3], kelp.richness$Richnessavg[kelp.richness$WA2==3],kelp.richness$HistoricRichness2[kelp.richness$WA2==2], kelp.richness$Richnessavg[kelp.richness$WA2==2], kelp.richness$HistoricRichness2[kelp.richness$WA2==1], kelp.richness$Richnessavg[kelp.richness$WA2==1], col=c(makeTransparent("blue",200), makeTransparent("orange",200),makeTransparent("blue",200), makeTransparent("orange",200),makeTransparent("blue",200), makeTransparent("orange",200)), ylab="Richness", ylim=c(0,12), las=1 )
boxplot(ab$average.abundanceHist[ab$WA2==3], ab$average.abundanceMod[ab$WA2==3],ab$average.abundanceHist[ab$WA2==2], ab$average.abundanceMod[ab$WA2==2], ab$average.abundanceHist[ab$WA2==1], ab$average.abundance2017[ab$WA2==1], col=c(makeTransparent("blue",200), makeTransparent("orange",200),makeTransparent("blue",200), makeTransparent("orange",200),makeTransparent("blue",200), makeTransparent("orange",200)), ylab="Average Abundance of Species", ylim=c(0,3), las=1 )
kelp.richness$RichRatio<-kelp.richness$Richnessavg/kelp.richness$HistoricRichness2
z<-kruskal.test(RichRatio~as.factor(WA2), data=kelp.richness)
z
dunn.test(kelp.richness$RichRatio,kelp.richness$WA2, method="Bonferroni")

t.test(kelp.richness$RichRatio[kelp.richness$WA==4]-1)

ab$AvRatio.rank<-ab$average.abundanceMod/ab$average.abundanceHist
z<-aov(rank(ab$AvRatio.rank)~as.factor(ab$WA2))
z<-kruskal.test((ab$AvRatio.rank)~as.factor(ab$WA2))
summary(z)
z
dunn.test(ab$AvRatio.rank,as.factor(ab$WA2), method="Bonferroni")

par(mar=c(4,4,1,1), mfrow=c(1,2))
boxplot(kelp.richness$RichRatio~kelp.richness$WA2,col=makeTransparent("navy", 65), las=1, ylim=c(0,2.5))
points(kelp.richness$RichRatio~kelp.richness$WA2, pch=19)
abline(1,0, lty=3, lwd=3, col="red")
boxplot(ab$AvRatio~ab$WA2, las=1, col=makeTransparent("navy", 65), ylim=c(0,2))
points(jitter(ab$AvRatio,10)~ab$WA2,pch=19)
abline(1,0, lty=3, lwd=3, col="red")

#Figure 4 (Species accumulation curves)
##All sites
library(vegan)
sac.data2<-read.csv("SAccRichness3.csv")
sac.data1993<-sac.data2[sac.data2$Year==1993,]
sac.data1995<-sac.data2[sac.data2$Year==1995,]
sac.data2017<-sac.data2[sac.data2$Year==2017,]
sac.data2018<-sac.data2[sac.data2$Year==2018,]

head(sac.data2)
z<-specaccum(sac.data1993[,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(z, ylim=c(0,20))
j<-specaccum(sac.data1995[,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(j, ylim=c(0,15))
j$richness
k<-specaccum(sac.data2017[,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(k, ylim=c(0,15))
k$richness
l<-specaccum(sac.data2018[,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(l, ylim=c(0,15))
l$richness





z1<-specaccum(sac.data1993[sac.data1993$WA==2,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(z1, ylim=c(0,15))
j1<-specaccum(sac.data1995[sac.data1995$WA==2,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(j1, ylim=c(0,15))
j1$richness
j1
k1<-specaccum(sac.data2017[sac.data2017$WA==2,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(j1, ylim=c(0,15))
k1$richness
k1
l1<-specaccum(sac.data2018[sac.data2018$WA==2,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(l1, ylim=c(0,15))
l1$richness
l1


z2<-specaccum(sac.data1993[sac.data1993$WA==3,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(z2, ylim=c(0,15))
j2<-specaccum(sac.data1995[sac.data1995$WA==3,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(j2, ylim=c(0,15))
j2$richness
j2
k2<-specaccum(sac.data2017[sac.data2017$WA==3,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(k2, ylim=c(0,15))
k2$richness
k2
l2<-specaccum(sac.data2018[sac.data2018$WA==3,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(l2, ylim=c(0,15))
l2$richness
l2


z3<-specaccum(sac.data1993[sac.data1993$WA==4,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(z3, ylim=c(0,15))
j3<-specaccum(sac.data1995[sac.data1995$WA==4,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(j3, ylim=c(0,15))
j3$richness
j3
k3<-specaccum(sac.data2017[sac.data2017$WA==4,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(k3, ylim=c(0,15))
k3$richness
k3
l3<-specaccum(sac.data2018[sac.data2018$WA==4,4:20], method="random", permutations=10000, ylim=c(1,15))
plot(l3, ylim=c(0,15))
l3$richness
l3

###Plot Figure 4
quartz()
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(z,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,15))
title("D",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(j,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(k,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(l,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))


plot(z3,ci=1,ci.type="polygon", col="black", lwd=2, xlab="",las=1, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,15))
title("D",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(k3,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(j3,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(l3,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))


plot(z2,ci=1,ci.type="polygon", col="black", lwd=2, xlab="",las=1, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,12))
title("D",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(j2,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(k2,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(l2,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))


plot(z1,ci=1,ci.type="polygon", col="black", lwd=2, xlab="",las=1, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,12))
title("D",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(j1,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(k1,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(l1,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))



####Estimating regional species pool
##ALL SITES
com93<-specpool(sac.data1993[,4:20])
com95<-specpool(sac.data1995[,4:20])
com17<-specpool(sac.data2017[,4:20])
com18<-specpool(sac.data2018[,4:20])


par(mar=c(4,4,1,1), mfrow=c(1,4))
X <- 1;
Y <- com93$boot[1];
SE <- (com93$boot.se[1]*1.96);
plot(X,Y, las=1,pch=19,cex=3,col=makeTransparent("white",150), xlim=c(0,5), ylim=c(0,20));
add.error.bars(X,Y,SE,0.1,col="black");
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",150))
points(X,Y, las=1,, cex=3, col="black")

X <- 2;
Y <- com95$boot[1];
SE <- (com95$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",200));
points(X,Y, las=1,, cex=3, col="black")

X <- 3;
Y <- com17$boot[1];
SE <- (com17$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",150));
points(X,Y, las=1,, cex=3, col="black")

X <- 4;
Y <- com18$boot[1];
SE <- (com18$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",200));
points(X,Y, las=1,, cex=3, col="black")


###Exposed
com93<-specpool(sac.data1993[sac.data1993$WA==4,4:20])
com95<-specpool(sac.data1995[sac.data1995$WA==4,4:20])
com17<-specpool(sac.data2017[sac.data2017$WA==4,4:20])
com18<-specpool(sac.data2018[sac.data2018$WA==4,4:20])


X <- 1;
Y <- com93$boot[1];
SE <- (com93$boot.se[1]*1.96);
plot(X,Y, las=1,pch=19,cex=3,col=makeTransparent("white",150), xlim=c(0,5), ylim=c(0,20));
add.error.bars(X,Y,SE,0.1,col="black");
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",150))
points(X,Y, las=1,, cex=3, col="black")

X <- 2;
Y <- com95$boot[1];
SE <- (com95$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",200));
points(X,Y, las=1,, cex=3, col="black")

X <- 3;
Y <- com17$boot[1];
SE <- (com17$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",150));
points(X,Y, las=1,, cex=3, col="black")

X <- 4;
Y <- com18$boot[1];
SE <- (com18$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",200));
points(X,Y, las=1,, cex=3, col="black")

###Moderate
com93<-specpool(sac.data1993[sac.data1993$WA==3,4:20])
com95<-specpool(sac.data1995[sac.data1995$WA==3,4:20])
com17<-specpool(sac.data2017[sac.data2017$WA==3,4:20])
com18<-specpool(sac.data2018[sac.data2018$WA==3,4:20])


X <- 1;
Y <- com93$boot[1];
SE <- (com93$boot.se[1]*1.96);
plot(X,Y, las=1,pch=19,cex=3,col=makeTransparent("white",150), xlim=c(0,5), ylim=c(0,20));
add.error.bars(X,Y,SE,0.1,col="black");
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",150))
points(X,Y, las=1,, cex=3, col="black")

X <- 2;
Y <- com95$boot[1];
SE <- (com95$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",200));
points(X,Y, las=1, cex=3, col="black")

X <- 3;
Y <- com17$boot[1];
SE <- (com17$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",150));
points(X,Y, las=1,, cex=3, col="black")

X <- 4;
Y <- com18$boot[1];
SE <- (com18$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",200));
points(X,Y, las=1,, cex=3, col="black")

###Sheltered
com93<-specpool(sac.data1993[sac.data1993$WA==2,4:20])
com95<-specpool(sac.data1995[sac.data1995$WA==2,4:20])
com17<-specpool(sac.data2017[sac.data2017$WA==2,4:20])
com18<-specpool(sac.data2018[sac.data2018$WA==2,4:20])


X <- 1;
Y <- com93$boot[1];
SE <- (com93$boot.se[1]*1.96);
plot(X,Y, las=1,pch=19,cex=3,col=makeTransparent("white",150), xlim=c(0,5), ylim=c(0,20));
add.error.bars(X,Y,SE,0.1,col="black");
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",150))
points(X,Y, las=1,, cex=3, col="black")

X <- 2;
Y <- com95$boot[1];
SE <- (com95$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("blue",200));
points(X,Y, las=1,, cex=3, col="black")

X <- 3;
Y <- com17$boot[1];
SE <- (com17$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",150));
points(X,Y, las=1,, cex=3, col="black")

X <- 4;
Y <- com18$boot[1];
SE <- (com18$boot.se[1]*1.96);
add.error.bars(X,Y,SE,0.1,col="black")
points(X,Y, las=1,pch=19, cex=3, col=makeTransparent("orange",200));
points(X,Y, las=1,, cex=3, col="black")



###SW facing sites - Exposure Validation
SW<-read.csv("SWSites.csv")
kruskal.test(logExp~as.factor(WA2),data=SW)
dunn.test(SW$logExp,SW$WA2)



##Moran's I

BS.GPSdata <- read.table("MoransIcalc.csv", sep=",", header=T)
BS.GPSdata[1:49,]->BS.GPSdata
head(BS.GPSdata, n=10)
dists <- as.matrix(dist(cbind(BS.GPSdata$Lon, BS.GPSdata$Lat)))

dists.inv <- 1/dists
diag(dists.inv) <- 0

dists.inv[1:5, 1:5]
Moran.I(BS.GPSdata$Ratio, dists)
