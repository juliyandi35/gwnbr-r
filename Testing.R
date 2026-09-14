# Input DBD
data <- read.csv("Pendidikan1.csv", sep = ";")
data
View(data)
summary(data[2:9])

# Mencari Nilai Variansi
attach(data) #Mengambil DBD per variabel dari objek DBD
var(data[2:9])
var(data[7])

#Mencari nila standar deviasi
sd(data$Putus)
sd(data$Miskin)
sd(data$IPM)
sd(data$RLS)
sd(data$HLS)
sd(data$APS)
sd(data$TPT)
sd(data$PDRB)

# Pengujian Multikolinearitas
library(car)
multiko=vif(lm(data$Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB))
multiko
cor(data[3:9])

#REGRESI POISSON#
modelpoisson=glm(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB,family=poisson,data=data)
summary(modelpoisson) 

#Menaksir Theta 
library(MASS)
nb=glm.nb(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB, data=data) 
summary(nb)

#REGRESI BINOMIAL NEGATIF# 
modelnegbin=glm(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB,family=
                  negative.binomial(2.042), data=data) 
summary(modelnegbin) 

#BPTES# 
library(lmtest)
depen=lm(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB, data=data) 
bptest(depen) 

#MORANS I# 
library(ape)
tugas.dists=as.matrix(dist(cbind(data$U, data$V))) 
tugas.dists.inv=1/tugas.dists
tugas.dists.inv[is.infinite(tugas.dists.inv)] <- 0
Moran.I(x=data$Putus, tugas.dists.inv) 

#Bandwitdh
library(spgwr)
bdwtBisquare=ggwr.sel(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB,data=data,
                      coords=cbind(data$U,data$V),adapt=TRUE,gweight=gwr.bisquare) 
GRTGB=ggwr(Putus~Miskin+IPM+RLS+HLS+APS+TPT+PDRB,data=data,
           coords=cbind(data$U,data$V),adapt=bdwtBisquare,gweight=gwr.bisquare) 
GRTGB$bandwidth

# Jarak Euclidean
u =data[,10]
u <- as.matrix(u)
i <- nrow(u)
v <- data[,11]
v <- as.matrix(v)
j <- nrow(v)
library(fields)
jarak <- matrix(nrow = 34, ncol = 34)
for(i in 1:34){
  for(j in 1:34){jarak[i,j]=sqrt((u[i,]-u[j,])**2+(v[i,]-v[j,])**2)}
}
write.table(jarak, file = "C:/Users/JULI YANDI RAHMAN/Downloads/Kerjaan/Selesai/GWNBR DL RStudio Sabtu/jarak.csv", sep = ";")
jarak

#PEMBOBOT# 
bdwtBisquare<- GRTGB$bandwidth
bdwtBisquare<- as.matrix(bdwtBisquare) 
bdwtBisquare 
i<-nrow(bdwtBisquare) 
pembobotB<-matrix(nrow=34,ncol=34) 
for(i in 1:34){
  for(j in 1:34) {
  pembobotB[i,j]=(1-(jarak[i,j]/bdwtBisquare[i,])^2)^2 
  pembobotB[i,j]<-ifelse(jarak[i,j]<bdwtBisquare[i,],pembobotB[i,j],0)
  } 
} 
write.table(pembobotB,file="C:/Users/JULI YANDI RAHMAN/Downloads/Kerjaan/Selesai/GWNBR DL RStudio Sabtu/pembobotB1.csv",sep=";")


library(MASS)
gemes2=function(X,y,W1,phi1,b1){ 
  beta=matrix(c(0),9,9,byrow=T) 
  beta[1,1]=phi1
  for(i in 1:9){ 
    satu<-rep(1,34) 
    satu<-as.matrix(satu) 
    b01<-rbind(c(phi1,beta[i,2:9])) 
    Xb1<-as.matrix(X)%*%as.matrix(beta[i,2:9]) 
    mu1<-exp(Xb1) 
    delta11<-((log(1+phi1*mu1)-digamma(y+(1/phi1))+digamma(1/phi1))/phi1^2)+((y-mu1)/((1+phi1*mu1)*phi1)) 
    delta11<-as.matrix(delta11) 
    p11<-t(satu)%*%W1%*%delta11 
    delta21<-(y-mu1)/(1+phi1*mu1) 
    delta21<-as.matrix(delta21) 
    p21<-t(X)%*%as.matrix(W1)%*%delta21 
    p21<-as.matrix(p21) 
    gt1<-rbind(p11,p21) 
    delta31<-((trigamma(y+(1/phi1))-
                 trigamma(1/phi1))/phi1^4)+((2*digamma(y+(1/phi1))-2*digamma(1/phi1)-2*log(1+phi1*mu1))/phi1^3)+
      ((2*mu1)/(phi1^2*(1+phi1*mu1)))+(((y+(1/phi1))*mu1^2)/(1+phi1*mu1)^2)-
      (y/phi1^2) 
    delta31<-as.matrix(delta31) 
    p31<-t(satu)%*%W1%*%delta31 
    p31<-as.matrix(p31) 
    delta41<-mu1*(mu1-y)/(1+phi1*mu1)^2
    delta41<-as.matrix(delta41) 
    p41<-t(X)%*%W1%*%delta41 
    p41<-as.matrix(p41) 
    h11<-rbind(p31,p41) 
    delta51<-mu1*(phi1*y+1)/(1+phi1*mu1)^2 
    delta51<-t(delta51) 
    delta51<-c(delta51) 
    delta51<-as.matrix(diag(delta51)) 
    p51<-t(X)%*%as.matrix(W1)%*%delta51%*%as.matrix(X) 
    p51<--1*p51 
    p51<-as.matrix(p51) 
    h21<-rbind(t(p41),p51) 
    H1<-cbind(h11,h21) 
    HI1<-ginv(H1) 
    beta[i,]<-(t(b01)-HI1%*%gt1) 
  } 
  return(list(beta=beta,hessian=H1)) 
} 
gwnbr1 <- function(x,y,W,teta){ 
  beta <- ginv(t(x)%*%x)%*%t(x)%*%y 
  param <- matrix(c(0),nrow(x),ncol(x)+1, byrow=T) 
  zhit <- matrix(c(0),nrow(x),ncol(x), byrow=T) 
  for(i in 1:34){ 
    ww <- as.matrix(diag(W[i,])) 
    hit <- gemes2(x,y,ww,teta,beta) 
    param[i,] <- hit$beta[9,] 
    write.csv(hit$hessian,file=paste("hessian",i,".csv")) 
    invh <- -ginv(as.matrix(hit$hessian)) 
    for(j in 1:ncol(x)){ 
      zhit[i,j] <- param[i,j] / invh[j+1,j+1] 
    } 
  }
  return(list(koefisien=param,Z_hitung=zhit)) 
} 


#Memanggil GWNBR
xx=data[,3:9]
y=data[,2]
x=as.matrix(cbind(1,xx))
mod=gwnbr1(x,y,pembobotB, 2.042) 
mod$Z_hitung
mod$koefisien 
write.table(mod$Z_hitung,file="C:/Users/JULI YANDI RAHMAN/Downloads/Kerjaan/Selesai/GWNBR DL RStudio Sabtu/z_hitungF1.csv",sep=";") 
write.table(mod$koefisien,file="C:/Users/JULI YANDI RAHMAN/Downloads/Kerjaan/Selesai/GWNBR DL RStudio Sabtu/koefisienF1.csv",sep=";")

datay=y
datax=x
DNB=26
mod=mod
tetagw<-as.matrix(mod$koefisien[,1])
betagw<-as.matrix(mod$koefisien[,2:9])
muogw<-as.matrix(exp(mod$koefisien[,2]))
muwgw<-as.matrix(exp(apply(datax*betagw,1,sum)))
slr<-matrix(0,nrow(data),1)
for (i in 1:nrow(data)) {
  Ri=datay[i]-1
  for (r in 0:Ri) {slr[i]<- log(1+r*tetagw)
  }
}
Lwgw <- sum(slr-log(factorial(log(datay)))+datay*log(muwgw)-(datay+1/tetagw)*log(tetagw*muwgw+1))
Logw <- sum(slr-log(factorial(log(datay)))+datay*log(muogw)-(datay+1/tetagw)*log(tetagw*muogw+1))
DGW <- 2*(Logw-Lwgw)
DGW

Fhit = 1/(DGW/26)
Fhit
