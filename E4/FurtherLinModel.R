### Statistical Modeling 2
### Further Longinitudial Models
### 2021

############################
### Example 4.1

data<-read.table("knee2.txt", sep="\t",header=TRUE)
attach(data)

Y<-cbind(R1, R2, R3, R4)
t<-c(0,3,7,10)
T<-cbind(rep(1,4),t,t^2)

### a)

Y.star<-Y%*%T%*%solve(t(T)%*%T)
model<-lm(Y.star~factor(Th)-1)
summary(model)
Theta<-coef(model)
Theta

x1<-t(cbind(1,0))
x2<-t(cbind(0,1))
mu1<-T%*%t(Theta)%*%x1
mu2<-T%*%t(Theta)%*%x2

plot(c(0,10),c(1,5), type="n", axes=FALSE, xlab="", ylab="")
points(t, mu1, col="red", pch=19)
lines(t, mu1, col="red")
points(t, mu2, col="blue", pch=19)
lines(t, mu2, col="blue")
axis(1, t)
axis(2, 1:5)

t.c<-seq(0,10,0.1)
T.c<-cbind(1,t.c,t.c^2)
mu2.continuous<-T.c%*%t(Theta)%*%x2
lines(t.c, mu2.continuous, col="blue", lwd=10)


### b)

xf<-t(cbind(0,1))
BP<-T%*%t(Theta)%*%xf

X<-model.matrix(model)
E<-Y-X%*%Theta%*%t(T)
Sigma<-t(E)%*%E/model$df

PT<-T%*%solve(t(T)%*%T)%*%t(T)
cov.error<-Sigma+(as.numeric(t(xf)%*%solve(t(X)%*%X)%*%xf)*PT%*%Sigma%*%PT)

lower.bound<-BP-qnorm(0.9)*diag(cov.error)^0.5
upper.bound<-BP+qnorm(0.9)*diag(cov.error)^0.5
lower.bound
upper.bound

cbind(t,lower.bound, BP,upper.bound)


## c)

Y.star<-Y%*%T%*%solve(t(T)%*%T)
model<-lm(Y.star~factor(Th)-1)
model.H0<-lm(Y.star~1)
Theta.H0<-coef(model.H0)
Theta.H0
mu.H0<-T%*%t(Theta.H0)
lines(t, mu.H0, col="black", lwd=10)

anova(model.H0, model)

## d)

T1<-T[,1:2]
Z1<-T1%x%X
Z<-T%x%X
y<-as.vector(Y)

model.H0<-lm(y~Z1-1)
model.H1<-lm(y~Z-1)

MSE.H0<-mean(residuals(model.H0)^2)
MSE.H1<-mean(residuals(model.H1)^2)
MSE.H0
MSE.H1

AIC(model.H0)
AIC(model.H1)

anova(model.H0, model.H1)   # p-value only approximate

Y.star1<-Y%*%T1%*%solve(t(T1)%*%T1)
model.1<-lm(Y.star1~factor(Th)-1)
Theta.1<-coef(model.1)
Theta.1
x1<-t(cbind(1,0))
x2<-t(cbind(0,1))
mu1.1<-T1%*%t(Theta.1)%*%x1
mu2.1<-T1%*%t(Theta.1)%*%x2


plot(c(0,10),c(1,5), type="n", axes=FALSE, xlab="", ylab="")
points(t, mu1.1, col="red", pch=19)
lines(t, mu1.1, col="red")
points(t, mu2.1, col="blue", pch=19)
lines(t, mu2.1, col="blue")
axis(1, t)
axis(2, 1:5)


############################
### Example 4.2

library(lme4)
library(longitudinal)
data(tcell)

tcell.34

id<-unlist(strsplit(rownames(tcell.34), "-"))[seq(2,680,2)]
time<-rep(attr(tcell.34,"time"),each=34)
Y<-cbind(tcell.34)
y<-as.vector(Y)
variable<-rep(colnames(Y), each=dim(Y)[1])
T<-rep(time,dim(Y)[2])
ID<-rep(id, dim(Y)[2])
data<-data.frame(ID, T, variable, y)
dd<-data[data$variable%in%colnames(Y)[c(2:4)],]

### a)

model<-lmer(y~-1+variable+variable:T+(-1+variable|ID), data=dd)
summary(model)
fixef(model)
B<-matrix(fixef(model), nrow=2,ncol=3, by=1)

newdata<-data.frame(variable="CCNG1", T=attr(tcell.34,"time"))
mu1<-predict(model, newdata=newdata, re.form=NA)
plot(time, Y[,2])
lines(newdata$T, mu1, lwd=3)

newdata<-data.frame(variable="TRAF5", T=attr(tcell.34,"time"))
mu2<-predict(model, newdata=newdata, re.form=NA)
plot(time, Y[,3])
lines(newdata$T, mu2, lwd=3)

newdata<-data.frame(variable="CLU", T=attr(tcell.34,"time"))
mu3<-predict(model, newdata=newdata, re.form=NA)
plot(time, Y[,4])
lines(newdata$T, mu3, lwd=3)


### b)

n<-dim(Y)[1]
X<-model.matrix(model)
head(X)

E<-matrix(residuals(model), nrow=dim(Y)[1], ncol=3, by=2)
Sigma<-t(E)%*%E/(n-2)
Sigma

diag(diag(Sigma)^{-0.5})%*%Sigma%*%diag(diag(Sigma)^{-0.5})

### c)

newdata<-data.frame(variable=c("CCNG1","TRAF5","CLU"), T=80, ID=34)
muf<-predict(model, newdata=newdata)
predict(model, newdata=newdata, re.form=NA)
ranef(model)

lower.bound<-qnorm(0.1, mean=muf, sd=(diag(Sigma)^(0.5)))
upper.bound<-qnorm(0.9, mean=muf, sd=(diag(Sigma)^(0.5)))

rbind(lower.bound,muf,upper.bound)


### d)

y1f<-c(17,18)
Sigma11<-Sigma[1:2,1:2]
Sigma12<-Sigma[1:2,3]
Sigma21<-Sigma[3,1:2]
Sigma22<-Sigma[3,3]

BP<-muf[3]+Sigma21%*%solve(Sigma11)%*%(y1f-muf[1:2])

cov.error<-Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12

lower.bound<-BP-qnorm(0.9)*diag(cov.error)^0.5
upper.bound<-BP+qnorm(0.9)*diag(cov.error)^0.5
lower.bound
upper.bound



############################
### Example GEE-Method

library(geepack)
data(respiratory)
respiratory
?respiratory

model<-geeglm(outcome ~ factor(center) + treat + age + baseline, data=respiratory, id=id, 
family=binomial(), corstr="unstructured")        
summary(model)

fitted(model, type="response")
newdata<-data.frame(center=2,treat="A",sex="M",age=31, baseline=1)
predict(model, newdata=newdata, type="response")

model$geese$alpha  # working correlation

model.H0<-geeglm(outcome ~  age + baseline, data=respiratory, id=id, family=binomial(), corstr="unstructured")        
anova(model.H0, model)





