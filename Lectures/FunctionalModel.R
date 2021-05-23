### Statistical Modeling 2
### Functional Data Models
### 2021

############################
### Example 5.1


library(refund)
data(DTI)
Y<-t(DTI$cca[-c(125,126,130,131,319,321),]) # removal of missing values and transpose of the data
t<-seq(0,1,length=dim(Y)[1])
matplot(t,Y, type="n")
matlines(t,Y, col="red")

library(fda)

# a)


## polynomial basis

basis <- create.monomial.basis(c(0,1),nbasis=10)
plot(basis)
Y.fd<-smooth.basis(t,Y,basis)
plot(Y.fd)

## B-spline basis

basis<-create.bspline.basis(c(0,1),nbasis=10, norder=4)
plot(basis)
Y.fd<-smooth.basis(t,Y, basis)
plot(Y.fd)

C<-Y.fd$fd$coef
C

fitted<-eval.fd(t,Y.fd$fd)

plot(fitted[,1],Y[,1])


Phi<-eval.basis(t,basis)
H.Phi<-Phi%*%solve(t(Phi)%*%Phi)%*%t(Phi)
fitted.2<-H.Phi%*%Y


# b)

basis<-create.bspline.basis(c(0,1),nbasis=100, norder=4)
basisPar<-fdPar(basis, Lfd=2, lambda=0.01)
Y.fd<-smooth.basis(t,Y, basisPar)
plot(Y.fd, col="red")

eval.fd(0.5, Y.fd$fd)[1]

c1<-Y.fd$fd$coef[,1]
phi.vector<-eval.basis(0.5,basis)
phi.vector%*%c1

sum(Y.fd$gcv)


# c)

loglambda<-seq(-8, 0, 0.25)
Gcvsave<-rep(NA, length(loglambda))
names(Gcvsave)<-loglambda
Dfsave<-Gcvsave
for(i in 1:length(loglambda)){
  basisPar<-fdPar(basis, Lfdobj=2, lambda=10^loglambda[i])
  Y.fd<-smooth.basis(t,Y, basisPar)
  Gcvsave[i]<-sum(Y.fd$gcv)
  Dfsave[i]<-Y.fd$df
}

Gcvsave

basis<-create.bspline.basis(c(0,1),nbasis=100, norder=4)
basisPar<-fdPar(basis, Lfd=2, lambda=10^(-7.25))
Y.fd<-smooth.basis(t,Y, basisPar)
plot(Y.fd, col="red")


############################
### Example 5.2

library(fda)
data(gait)
Y<-gait[,,"Knee Angle"]
t<-seq(0.025,0.975, 0.05)

matplot(t,Y, type="n")
matlines(t,Y)

basis<-create.fourier.basis(c(0.025,0.975), nbasis=19)
plot(basis)
Y.fd<-smooth.basis(t,Y,basis)
plot(Y.fd)


basis<-create.fourier.basis(c(0.025,0.975), nbasis=19)
Lcoef<-c(0,(2*pi/diff(c(0.025,0.975)))^2,0)
harmaccelLfd<-vec2Lfd(Lcoef, c(0.025,0.975))
fdParobj<-fdPar(basis, harmaccelLfd, lambda=10^(-8))
Y.fd<-smooth.basis(t,Y,fdParobj)
plot(Y.fd)


loglambda<-seq(-8, 0, 0.25)
Gcvsave<-rep(NA, length(loglambda))
names(Gcvsave)<-loglambda
Dfsave<-Gcvsave
for(i in 1:length(loglambda)){
 	fdParobj<-fdPar(basis, harmaccelLfd, lambda=10^loglambda[i])
	Y.fd<-smooth.basis(t,Y,fdParobj)
  	Gcvsave[i]<-sum(Y.fd$gcv)
}

Gcvsave

basis<-create.fourier.basis(c(0.025,0.975), nbasis=19)
Lcoef<-c(0,(2*pi/diff(c(0.025,0.975)))^2,0)
harmaccelLfd<-vec2Lfd(Lcoef, c(0.025,0.975))
fdParobj<-fdPar(basis, harmaccelLfd, lambda=10^(-7.75))
Y.fd<-smooth.basis(t,Y,fdParobj)
plot(Y.fd)

fitted<-eval.fd(t,Y.fd$fd)
plot(fitted[,1],Y[,1])
points(fitted[,3],Y[,3], pch=19, col="red")

time<-seq(0.025,0.975, 0.001)
D.fit<-eval.fd(time, Y.fd$fd,1) # First Derivates


plot(c(0,1), c(-600,600), type="n", main="Derivates")
for(j in 1:dim(D.fit)[2]){
lines(time,D.fit[,j], col="orange", lwd=2)
}

par(mfrow=c(1,2))
plot(Y.fd)
plot(c(0,1), c(-600,600), type="n", main="Derivates")
for(j in 1:dim(D.fit)[2]){
lines(time,D.fit[,j], col="orange", lwd=2)
}



############################
### Example 5.3

library(fda)

X<-as.matrix(read.table("stagediameter.txt", sep="\t", header=TRUE))
Y<-as.matrix(read.table("stageheight.txt", sep="\t", header=TRUE))
t<-seq(10,110,10)
head(X)
head(Y)
matplot(t,X,type="n")
matlines(t,X,col="blue")
matplot(t,Y,type="n")
matlines(t,Y,col="red")

## a)

basis<-create.bspline.basis(c(10,110),nbasis=4, norder=4)
Y.fd<-smooth.basis(t,Y,basis)
X.fd<-smooth.basis(t,X,basis)

plot(X.fd, col="blue", add=TRUE)
plot(Y.fd, col="red", add=TRUE)

par(mfrow=c(1,2))
matplot(t,Y,type="n")
matlines(t,Y,col="red")
plot(Y.fd, col="red", add=TRUE)

## b)

mean.fd<-mean(Y.fd$fd)
plot(Y.fd)
lines(mean.fd, col="red", lwd=4)

eval.fd(100, mean.fd)

## c)

time.f<-seq(110,120,0.1)
eval.fd(time.f, Y.fd$fd)

basis.f<-create.bspline.basis(c(10,120),nbasis=4, norder=4)
Yf.fd<-smooth.basis(t,Y,basis.f)

Phi.f<-eval.basis(time.f,basis.f)

C.hat<-Yf.fd$y2cMap%*%Y

Pred<-Phi.f%*%C.hat

plot(c(10,120),c(0,70), type="n")
lines(Yf.fd)

for(i in 1: dim(Y)[2]){

lines(time.f,Pred[,i], col="brown", lwd=3)

}


## d)


y1f<-Y[1:4,1]

mean.fd<-mean(Y.fd$fd)
mu<-eval.fd(t, mean.fd)
mu1<-mu[1:4]
mu2<-mu[5:11]
ones<-rbind(rep(1,dim(Y)[2]))

V<-Y-(mu%*%ones)
Sigma<-V%*%t(V)/(dim(Y)[1]-1)

Sigma11<-Sigma[1:4,1:4]
Sigma12<-Sigma[1:4,5:11]
Sigma21<-Sigma[5:11,1:4]
Sigma22<-Sigma[5:11,5:11]

y2f.hat<-mu2+Sigma21%*%solve(Sigma11)%*%(y1f-mu1)

yf.fd<-smooth.basis(t,c(y1f,y2f.hat), basis)

plot(yf.fd)
lines(mean.fd, col="red")


############################
### Example 5.4

library(fda)
library(refund)
data(DTI)
Y<-t(DTI$cca[-c(125,126,130,131,319,321),]) # removal of missing values and transpose of the data
t<-seq(0,1,length=dim(Y)[1])
matplot(t,Y, type="n")
matlines(t,Y, col="red")

## a)

basis<-create.bspline.basis(c(0,1),nbasis=100, norder=4)
basisPar<-fdPar(basis, Lfd=2, lambda=10^(-7.25))
Y.fd<-smooth.basis(t,Y, basisPar)
plot(Y.fd, col="red")

X1<-DTI$case[-c(125,126,130,131,319,321)]
X2<-as.numeric(DTI$sex[-c(125,126,130,131,319,321)]=="female")  # denoting females

## b)

X.list<-vector("list",3)
X.list[[1]]<-rep(1,dim(Y)[2])
X.list[[2]]<-X1
X.list[[3]]<-X2

beta0.basis<-create.bspline.basis(c(0,1),nbasis=50, norder=6)
beta1.basis<-create.constant.basis(rangeval=c(0,1), names="const")
beta2.basis<-create.constant.basis(rangeval=c(0,1), names="const")

beta.list<-vector("list",3)
beta.list[[1]]<-beta0.basis
beta.list[[2]]<-beta1.basis
beta.list[[3]]<-beta2.basis

model<-fRegress(Y.fd$fd,X.list,beta.list)

Yhat.fd<-model$yhatfd

plot(Y.fd, type="n")
lines(Y.fd, col="red")
lines(Yhat.fd, col="black",lwd=4)

beta0.fd<-model$betaestlist[[1]]
beta1.fd<-model$betaestlist[[2]]
beta2.fd<-model$betaestlist[[3]]

plot(beta0.fd)
plot(beta1.fd)
plot(beta2.fd)

time<-seq(0,1,0.001)
fit.11<-eval.fd(time,beta0.fd$fd)+eval.fd(time,beta1.fd$fd)+eval.fd(time,beta2.fd$fd)
fit.00<-eval.fd(time,beta0.fd$fd)
plot(Y.fd, type="n")
lines(Y.fd, col="red")
lines(time, fit.11, col="black", lwd=5)
lines(time, fit.00, col="blue", lwd=5)

fitFD.11<-beta0.fd$fd+beta1.fd$fd+beta2.fd$fd

############################
### Example 5.5


X<-as.matrix(read.table("stagediameter.txt", sep="\t", header=TRUE))
Y<-as.matrix(read.table("stageheight.txt", sep="\t", header=TRUE))
t<-seq(10,110,10)
head(X)
head(Y)
matplot(t,X,type="n")
matlines(t,X,col="blue")
matplot(t,Y,type="n")
matlines(t,Y,col="red")

## a)

basis<-create.bspline.basis(c(10,110),nbasis=4, norder=4)
Y.fd<-smooth.basis(t,Y,basis)
X.fd<-smooth.basis(t,X,basis)
plot(Y.fd)
plot(X.fd)

## b)

X.list<-vector("list",2)
X.list[[1]]<-rep(1,18)
X.list[[2]]<-X.fd$fd

beta0.basis<-create.constant.basis(rangeval=c(10,110), names="const")
#beta0.basis<-create.bspline.basis(c(10,110),nbasis=10, norder=4)
beta1.basis<-create.constant.basis(rangeval=c(10,110), names="const")
#beta1.basis<-create.monomial.basis(rangeval=c(10,110), nbasis=2)

beta.list<-vector("list",2)
beta.list[[1]]<-beta0.basis
beta.list[[2]]<-beta1.basis

model<-fRegress(Y.fd$fd,X.list,beta.list)

beta0.fd<-model$betaestlist[[1]]
beta1.fd<-model$betaestlist[[2]]
plot(beta0.fd)
plot(beta1.fd)

Yhat.fd<-model$yhatfd
plot(Y.fd$fd)
lines(Yhat.fd, col="black", lwd=5, lty=1)
lines(mean(Y.fd$fd), col="red", lwd=7)

xstar<-X[,1]
xstar.fd<-smooth.basis(t,xstar,basis)
fit.fd<-beta0.fd$fd+(beta1.fd$fd*xstar.fd$fd)

lines(fit.fd,col="blue", lwd=7)



