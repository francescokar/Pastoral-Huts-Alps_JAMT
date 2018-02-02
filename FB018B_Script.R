## Load Required Packages ##
library(spatstat)
library(MuMIn)
library(pgirmess)
library(spdep)
library(gstat)
library(sp)


# Importing table

header<-c("ID","OBJECT","MATERIAL","DIMENSION","PRESERVATION","X","Y")
RR_FB018B<-read.table(file.choose(),header=F,col.names=header)


# Creating Region of Interest and dataset

FB018B<-owin(poly=list(x=c(53.14,-82.03,-108.79,-122.66,-157.03,-184.77,-105.27,-65.23,-104.45,1.96,28.33,84.20,103.34,164.87,209.99,219.37,133.03),
                       y=c(578.26,349.53,309.87,264.56,222.56,183.10,135.83,123.53,17.07,1.25,94.62,155.37,245.61,312.22,483.52,497.59,558.14)))
RR_FB018B_COORD<-ppp(RR_FB018B[,6],RR_FB018B[,7],window=FB018B)


# Global Homogeneous L-Function

RR_FB018B_LEST<-envelope(RR_FB018B_COORD,fun=Lest,nsim=9999,correction="Ripley")
plot(RR_FB018B_LEST,main="FB018B",xlab="distance(cm)")


# Local L-Function

RR_FB018B_LOCALL_5<-localL(RR_FB018B_COORD,correction="Ripley",rvalue=5)
RR_FB018B_LOCALL_15<-localL(RR_FB018B_COORD,correction="Ripley",rvalue=15)
RR_FB018B_LOCALL_30<-localL(RR_FB018B_COORD,correction="Ripley",rvalue=30)
RR_FB018B_LOCALL_40<-localL(RR_FB018B_COORD,correction="Ripley",rvalue=40)

# Local L maps

RR_FB018B_LOCALL<-data.frame(ID=RR_FB018B$ID,X=RR_FB018B$X,Y=RR_FB018B$Y,LOCALL5=RR_FB018B_LOCALL_5,
                             LOCALL15=RR_FB018B_LOCALL_15,LOCALL30=RR_FB018B_LOCALL_30,LOCALL40=RR_FB018B_LOCALL_40)

RR_FB018B_LOCALL_FILTER5<-subset(RR_FB018B_LOCALL,LOCALL5>5)
RR_FB018B_COORD_LOCALL05<-ppp(RR_FB018B_LOCALL_FILTER5[,2],RR_FB018B_LOCALL_FILTER5[,3],
                              window=FB018B,marks=RR_FB018B_LOCALL_FILTER5[,4])
plot(RR_FB018B_COORD_LOCALL05,markscale=0.3,cols=c("red"),main=" ") 

RR_FB018B_LOCALL_FILTER15<-subset(RR_FB018B_LOCALL,LOCALL15>15)
RR_FB018B_COORD_LOCALL15<-ppp(RR_FB018B_LOCALL_FILTER15[,2],RR_FB018B_LOCALL_FILTER15[,3],
                              window=FB018B,marks=RR_FB018B_LOCALL_FILTER15[,5])
plot(RR_FB018B_COORD_LOCALL15,markscale=0.3,cols=c("red"),main=" ")

RR_FB018B_LOCALL_FILTER30<-subset(RR_FB018B_LOCALL,LOCALL30>30)
RR_FB018B_COORD_LOCALL30<-ppp(RR_FB018B_LOCALL_FILTER30[,2],RR_FB018B_LOCALL_FILTER30[,3],
                              window=FB018B,marks=RR_FB018B_LOCALL_FILTER30[,6])
plot(RR_FB018B_COORD_LOCALL30,markscale=0.3,cols=c("red"),main=" ")

RR_FB018B_LOCALL_FILTER40<-subset(RR_FB018B_LOCALL,LOCALL40>40)
RR_FB018B_COORD_LOCALL40<-ppp(RR_FB018B_LOCALL_FILTER40[,2],RR_FB018B_LOCALL_FILTER40[,3],
                              window=FB018B,marks=RR_FB018B_LOCALL_FILTER40[,7])
plot(RR_FB018B_COORD_LOCALL40,markscale=0.2,cols=c("red"),main=" ")


# Inhomogeneous Possion model

X<-c(-180.42)
Y<-c(184.14)
FIREPLACE_FB018B<-ppp(X,Y,window=FB018B)
DIST_FIREPLACE_FB018B<-distmap(FIREPLACE_FB018B)
FB018B_POISSON_MODEL<-ppm(RR_FB018B_COORD,~D,covariates=list(D=DIST_FIREPLACE_FB018B))


# Testing Inhomogeneity

FB018B_NULL_MODEL<-ppm(RR_FB018B_COORD,~1)
AICc(FB018B_POISSON_MODEL,FB018B_NULL_MODEL)
Weights(AICc(FB018B_POISSON_MODEL,FB018B_NULL_MODEL))


# L-Function conditioned on the Inhomogeneous Poisson model

RR_FB018B_LINHOM<-envelope(FB018B_POISSON_MODEL,Lest,nsim=9999,correction="Ripley")
plot(RR_FB018B_LINHOM,main="FB018B: Inhomogeneous Poisson",xlab="distance(cm)")


# Moran's I correlogram

RR_FB018B_XY<-cbind(RR_FB018B$X,RR_FB018B$Y)
RR_FB018B_MORAN_CORRELOG<-correlog(coords=RR_FB018B_XY,z=log(RR_FB018B$DIMENSION),method="Moran")
RR_FB018B_MORAN_CORRELOG_TABLE<-as.data.frame(RR_FB018B_MORAN_CORRELOG,row.names=NULL)
plot(RR_FB018B_MORAN_CORRELOG_TABLE$dist.class,RR_FB018B_MORAN_CORRELOG_TABLE$coef,type="b",pch=0,cex=2,
     xlim=c(0,290),ylim=c(-0.12,0.27),xlab="Distance (cm)",ylab="Moran's I",main="FB018B: Moran's I Correlogram")
par(new=T)
plot(RR_FB018B_MORAN_CORRELOG_TABLE$dist.class[2],RR_FB018B_MORAN_CORRELOG_TABLE$coef[2],pch=15,cex=2,
     xlim=c(0,290),ylim=c(-0.12,0.27),ann=F)


# Local Moran's I test of RR

RR_FB018B_KNEARNEIGH<-knearneigh(RR_FB018B_XY,k=4,RANN=F)
RR_FB018B_KNN2NB<-knn2nb(RR_FB018B_KNEARNEIGH)
RR_FB018B_LISTW<-nb2listw(RR_FB018B_KNN2NB)
RR_FB018B_LOCALMORAN<-localmoran(log(RR_FB018B$DIMENSION),RR_FB018B_LISTW)
xy<-RR_FB018B[6:7]
moran<-as.data.frame(RR_FB018B_LOCALMORAN[1:74])
RR_FB018B_SPDF_MORAN<-SpatialPointsDataFrame(coords=xy,data=moran)
bubble(RR_FB018B_SPDF_MORAN,"RR_FB018B_LOCALMORAN.1.74.",main="FB018B Local Moran's I")


# Variogram

RR_FB018B_SPDF<-SpatialPointsDataFrame(coords=RR_FB018B[6:7],data=log(RR_FB018B[4]))
RR_FB018B_VARIOGRAM<-variogram(DIMENSION~1,RR_FB018B_SPDF,cutoff=120,width=13)
plot(gamma~dist,RR_FB018B_VARIOGRAM,pch=16,xlab="distance(cm)",ylab="semivariance",col="black")
lines(RR_FB018B_VARIOGRAM$dist,RR_FB018B_VARIOGRAM$gamma,lty=2,lwd=1.5)

# Directional variogram

RR_FB018B_DIRECTVARIOGRAM<-variogram(DIMENSION~1,RR_FB018B_SPDF,alpha=c(0,45,90,135),cutoff=120,width=13)
plot(RR_FB018B_DIRECTVARIOGRAM,as.table=T,pch=16,lty=2,lwd=1.5,xlab="distance(cm)",col="black")

# Variogram map

RR_FB018B_VARIOGRAMMAP<-variogram(DIMENSION~1,RR_FB018B_SPDF,cutoff=120,width=13,map=T)
plot(RR_FB018B_VARIOGRAMMAP,ylab="dy(cm)",xlab="dx(cm)",main="RR FB018B",col.regions=heat.colors(100))