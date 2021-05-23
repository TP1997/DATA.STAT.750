### Statistical Modeling 2
### Multivariate Linear Model
### 2021

############################
### Example 3.1


library(MASS)
data(Cars93)
data<-Cars93[,c(4:6,13,25)]
head(data)

attach(data)

Y<-cbind(Min.Price,Price,Max.Price)
model<-lm(Y~Horsepower+Weight)
summary(model)

## a)

B.hat<-coef(model)
B.hat

## b)

E<-residuals(model)
Sigma<-(t(E)%*%E)/model$df
Sigma

## c)

newdata<-data.frame(Horsepower=200,Weight=3500)
BP<-predict(model, newdata=newdata)
BP

X<-model.matrix(model)
xf<-t(cbind(1,newdata))
cov.error<-Sigma*as.numeric(1+t(xf)%*%solve(t(X)%*%X)%*%xf)

lower.bound<-BP-qnorm(0.9)*diag(cov.error)^0.5
upper.bound<-BP+qnorm(0.9)*diag(cov.error)^0.5
lower.bound
upper.bound

## d)

yf1<-c(25,29.50)
Bt<-t(coef(model))
B1t<-Bt[1:2,]
B2t<-Bt[3,]
Sigma11<-Sigma[1:2,1:2]
Sigma12<-Sigma[1:2,3]
Sigma21<-Sigma[3,1:2]
Sigma22<-Sigma[3,3]

BP<-B2t%*%xf+Sigma21%*%solve(Sigma11)%*%(yf1-B1t%*%xf)
BP
cov.error<-Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12

lower.bound<-BP-qnorm(0.9)*diag(cov.error)^0.5
upper.bound<-BP+qnorm(0.9)*diag(cov.error)^0.5
lower.bound
upper.bound



############################
### Example 3.2


## a)

Y<-cbind(Min.Price,Price,Max.Price)
model<-lm(Y~Horsepower+Weight)
summary(model)

model.H0<-lm(Y~Horsepower)
summary(model.H0)

anova(model.H0, model)

library(car)
Anova(model)

E<-residuals(model)
EH<-residuals(model.H0)

S<-t(E)%*%E
SH<-t(EH)%*%EH

V<-dim(Y)[2]-sum(diag(solve(SH%*%solve(S))))  ## Pillai
V

n<-dim(data)[1]	## sample size
r<-3  		## rank(X)
d<-1  		## rank(X2)
q<-dim(coef(model))[2]	## number of columns in B
s<-min(d,q)

F<-(V/(s-V)*((n-r-q+s)/(abs(q-d)+s)))
F

df1<-s*(abs(q-d) + s)
df2<-s*(n - r - q + s)
pvalue<-1-pf(F, df1=df1, df2=df2)
pvalue

## b)

model.H1<-lm(cbind(Min.Price, Price)~Horsepower+Weight)
model.H0<-lm(cbind(Min.Price, Price)~Horsepower)

anova(model.H0, model.H1)

## c)

model.H1<-lm(cbind(Price,Max.Price)~Horsepower+Weight+Min.Price)
summary(model.H1)

model.H0<-lm(cbind(Price,Max.Price)~Min.Price)
summary(model.H0)

anova(model.H0, model.H1)


############################
### Example 3.3

data<-read.table("knee2.txt", sep="\t",header=TRUE)
attach(data)

Y<-cbind(R1, R2, R3, R4)


plot(c(1,4),c(1,5), type="n")
for(i in 1:length(N)){

lines(jitter(1:4),jitter(Y[i,]), col="red")

}

model<-lm(Y~factor(Th))
summary(model)

newdata<-data.frame(Th=c(1,2))

pred<-predict(model, newdata=newdata)

plot(c(1,4),c(1,5), type="n", axes=FALSE, xlab="", ylab="")
points(1:4, pred[1,], col="red", pch=19)
lines(1:4, pred[1,], col="red")
points(1:4, pred[2,], col="blue", pch=19)
lines(1:4, pred[2,], col="blue")
axis(1, 1:4)
axis(2, 1:5)


## a)

Y<-cbind(R1, R2, R3, R4)
model.H0<-lm(Y~1)
anova(model.H0, model)

predH0<-predict(model.H0, newdata=newdata)
points(1:4, predH0[2,], col="black", pch=19)
lines(1:4, predH0[2,], col="black")

## b)

L<-t(t((1/4)*rep(1,4)))
YL<-Y%*%L
head(Y)
head(YL)

model<-lm(YL~factor(Th))
summary(model)

model.H0<-lm(YL~1)
anova(model.H0, model)


## c)

L1<-c(1,-1,0,0)
L2<-c(0,1,-1,0)
L3<-c(0,0,1,-1)
L<-cbind(L1,L2,L3)

YL<-Y%*%L
head(Y)
head(YL)

model<-lm(YL~factor(Th))
summary(model)

library(MASS)

K<-t(c(1,1))
T<-diag(2)-ginv(K)%*%K
T<-round(T,1)

X<-model.matrix(model)
X.star<-X%*%T

model.H0<-lm(YL~X.star-1)
summary(model.H0)
anova(model.H0,model)

## d)

L1<-c(1,-1,0,0)
L2<-c(0,1,-1,0)
L3<-c(0,0,1,-1)
L<-cbind(L1,L2,L3)

YL<-Y%*%L
head(Y)
head(YL)

model<-lm(YL~factor(Th))
summary(model)

model.H0<-lm(YL~1)
anova(model.H0,model)

## testing the same by general hypotheses testing procedure

K<-t(c(0,1))
T<-diag(2)-ginv(K)%*%K

X<-model.matrix(model)
X.star<-X%*%T

model.H0<-lm(YL~X.star-1)
anova(model.H0,model)




