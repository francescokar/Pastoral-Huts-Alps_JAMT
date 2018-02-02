## Load required packages ##
library(spatstat)
library(gstat)
library(sp)
library(pgirmess)


# Importing table

header<-c("ID","OBJECT","MATERIAL","DIMENSION","PRESERVATION","X","Y")
RR_FB018A<-read.table(file.choose(),header=F,col.names=header)


# Creating Region of Interest and dataset

FB018A<-owin(poly=list(x=c(-156.62,-18.54,-9.07,47.09,126.14,62.12,158.49,257.79,188.31,11.73),
                       y=c(476.33,180.96,79.96,11.78,29.12,119.95,163.63,274.80,364.47,554.45)))
RR_FB018A_COORD<-ppp(RR_FB018A[,6],RR_FB018A[,7],window=FB018A)


# Global Homogeneous L-Function (9999 permutations)

RR_FB018A_LEST<-envelope(RR_FB018A_COORD,fun=Lest,nsim=9999,correction="Ripley")
plot(RR_FB018A_LEST,main="FB018A",xlab="distance(cm)")


# Moran's I correlogram

RR_FB018A_XY<-cbind(RR_FB018A$X,RR_FB018A$Y)
RR_FB018A_MORAN_CORRELOG<-correlog(coords=RR_FB018A_XY,z=log(RR_FB018A$DIMENSION),method="Moran")
RR_FB018A_MORAN_CORRELOG_TABLE<-as.data.frame(RR_FB018A_MORAN_CORRELOG,row.names=NULL)
plot(RR_FB018A_MORAN_CORRELOG_TABLE$dist.class,RR_FB018A_MORAN_CORRELOG_TABLE$coef,type="b",pch=0,cex=2,
     xlim=c(0,250),ylim=c(-0.30,0.50),xlab="Distance (cm)",ylab="Moran's I",main="FB018A: Moran's I Correlogram")
par(new=T)
plot(RR_FB018A_MORAN_CORRELOG_TABLE$dist.class[1],RR_FB018A_MORAN_CORRELOG_TABLE$coef[1],
     pch=15,cex=2,xlim=c(0,250),ylim=c(-0.30,0.50),ann=F)


# Variogram

RR_FB018A_SPDF<-SpatialPointsDataFrame(coords=RR_FB018A[6:7],data=log(RR_FB018A[4]))
RR_FB018A_VARIOGRAM<-variogram(DIMENSION~1,RR_FB018A_SPDF,cutoff=135,width=13)
plot(gamma~dist,RR_FB018A_VARIOGRAM,pch=16,xlab="distance(cm)",ylab="semivariance",col="black")
lines(RR_FB018A_VARIOGRAM$dist,RR_FB018A_VARIOGRAM$gamma,lty=2,lwd=1.5)

# Directional variogram

RR_FB018A_DIRECTVARIOGRAM<-variogram(DIMENSION~1,RR_FB018A_SPDF,alpha=c(0,45,90,135),cutoff=135,width=13)
plot(RR_FB018A_DIRECTVARIOGRAM,as.table=T,pch=16,lty=2,lwd=1.5,xlab="distance(cm)",col="black")

# Variogram map

RR_FB018A_VARIOGRAMMAP<-variogram(DIMENSION~1,RR_FB018A_SPDF,cutoff=135,width=13,map=T)
plot(RR_FB018A_VARIOGRAMMAP,ylab="dy(cm)",xlab="dx(cm)",main="FB018A",col.regions=heat.colors(100))