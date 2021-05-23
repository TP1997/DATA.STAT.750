### Statistical Modeling 2
### Generalized Linear Mixed Models
### 2021

############################
### Example 2.1


data<-read.table("retinal.txt", sep="\t", dec=".", header=TRUE)
attach(data)

plot(Time, Gas)


plot(Time, Gas, type="n")
for(i in 1:length(levels(factor(ID)))){

osa<-data[ID==i,]

if(osa$Level[1]=="15%"){
lines(osa$Time, osa$Gas, lwd=2, col="green")
}
if(osa$Level[1]=="20%"){
lines(osa$Time, osa$Gas, lwd=2, col="blue")
}
if(osa$Level[1]=="25%"){
lines(osa$Time, osa$Gas, lwd=2, col="red")
}}

## a)

library(lme4)

model.simple<-lmer(Gas~Level*Time+(1|ID), data=data)
summary(model.simple)
fixef(model.simple)
ranef(model.simple)

pred.data<-data.frame(data, pred=predict(model.simple, newdata=data, re.form=NA))
pred.data<-data.frame(data, pred=predict(model.simple, newdata=data, type="response"))

plot(Time, Gas, type="n")
for(i in 1:length(levels(factor(ID)))){
osa<-pred.data[ID==i,]
if(osa$Level[1]=="15%"){
lines(osa$Time, osa$pred, lwd=2, col="green")
}
if(osa$Level[1]=="20%"){
lines(osa$Time, osa$pred, lwd=2, col="blue")
}
if(osa$Level[1]=="25%"){
lines(osa$Time, osa$pred, lwd=2, col="red")
}}

##

model<-lmer(Gas~Level*Time+(1+Time|ID), data=data)
summary(model)

pred.data<-data.frame(data, pred=predict(model, newdata=data))
plot(Time, Gas, type="n")
for(i in 1:length(levels(factor(ID)))){
osa<-pred.data[ID==i,]
if(osa$Level[1]=="15%"){
lines(osa$Time, osa$pred, lwd=2, col="green")
}
if(osa$Level[1]=="20%"){
lines(osa$Time, osa$pred, lwd=2, col="blue")
}
if(osa$Level[1]=="25%"){
lines(osa$Time, osa$pred, lwd=2, col="red")
}}


model.log<-glmer(Gas~Level*Time+(1+Time|ID),family=gaussian(link="log"), data=data)
summary(model.log)

model.inverse<-glmer(Gas~Level*Time+(1+Time|ID),family=gaussian(link="inverse"), data=data)
summary(model.inverse)

pred.data<-data.frame(data, pred=predict(model.inverse, newdata=data, type="response"))
plot(Time, Gas, type="n")
for(i in 1:length(levels(factor(ID)))){
osa<-pred.data[ID==i,]
if(osa$Level[1]=="15%"){
lines(osa$Time, osa$pred, lwd=2, col="green")
}
if(osa$Level[1]=="20%"){
lines(osa$Time, osa$pred, lwd=2, col="blue")
}
if(osa$Level[1]=="25%"){
lines(osa$Time, osa$pred, lwd=2, col="red")
}}

AIC(model)
AIC(model.log)
AIC(model.inverse)


model.Gamma<-glmer(Gas~Level*Time+(1+Time|ID),family=Gamma(link="log"), data=data)
summary(model.Gamma)

model.IG<-glmer(Gas~Level*Time+(1+Time|ID),family=inverse.gaussian(link="log"), data=data)
summary(model.IG)



pred.data<-data.frame(data, pred=predict(model.Gamma, newdata=data, type="response"))
plot(Time, Gas, type="n")
for(i in 1:length(levels(factor(ID)))){
osa<-pred.data[ID==i,]
if(osa$Level[1]=="15%"){
lines(osa$Time, osa$pred, lwd=2, col="green")
}
if(osa$Level[1]=="20%"){
lines(osa$Time, osa$pred, lwd=2, col="blue")
}
if(osa$Level[1]=="25%"){
lines(osa$Time, osa$pred, lwd=2, col="red")
}}


## b)

newdata<-data.frame(ID=12, Level="15%", Time=20)
mu.hat<-predict(model.Gamma, newdata=newdata, type="response", re.form=NA)
mu.hat

mu.pred<-predict(model.Gamma, newdata=newdata, type="response")
mu.pred

## normal distribution

mu.hat<-predict(model.log, newdata=newdata, type="response", re.form=NA)
mu.hat

mu.pred<-predict(model.log, newdata=newdata, type="response")
mu.pred



## c)

model.H0<-glmer(Gas~Time+(1+Time|ID),family=Gamma(link="log"), data=data)
anova(model.H0, model.Gamma)


model.H0<-glmer(Gas~Time+(1+Time|ID),family=gaussian(link="log"), data=data)
anova(model.H0, model.log)



## d)

VarCorr(model.Gamma)
print(VarCorr(model.Gamma), comp="Variance")

phi<-sigma(model.Gamma)^2

ST<-getME(model.Gamma, "ST")$ID
S<-diag(diag(ST))
T<-ST-S+diag(2)
F<-T%*%S%*%S%*%t(T)
cov.bi<-phi*F
cov.bi

## normal distribution

sigma2<-sigma(model.log)^2

ST<-getME(model.log, "ST")$ID
S<-diag(diag(ST))
T<-ST-S+diag(2)
F<-T%*%S%*%S%*%t(T)
cov.bi<-sigma2*F
cov.bi



## e)

### Parametric bootstrap method

model.Gamma<-glmer(Gas~Level*Time+(1+Time|ID),family=Gamma(link="log"), data=data)

n<-dim(data)[1]
newdata<-data.frame(ID=12, Level="15%", Time=20)

mupred<-predict(model.Gamma, newdata=data, type="response")
varerror<-(sigma(model.Gamma)^2)*mupred^2

a.hat<-(mupred^2)/varerror
s.hat<-varerror/mupred

yf.star<-numeric()

for(b in 1:100){

yb<-rgamma(n, shape=a.hat,scale=s.hat)
model.b<-glmer(yb~Level*Time+(1+Time|ID),family=Gamma(link="log"), data=data)
mupred.star<-predict(model.b, newdata=newdata, type="response")
varerror.star<-(sigma(model.b)^2)*mupred.star^2
a.star<-(mupred.star^2)/varerror.star
s.star<-varerror.star/mupred.star
yf.star[b]<-rgamma(1, shape=a.star,scale=s.star)

}

lower.bound<-quantile(yf.star, 0.1)
upper.bound<-quantile(yf.star, 0.9)
lower.bound
upper.bound

## normal distribution 

n<-dim(data)[1]
newdata<-data.frame(ID=12, Level="15%", Time=20)

mupred<-predict(model.log, newdata=data, type="response")
varerror<-(sigma(model.log)^2)
yf.star<-numeric()

for(b in 1:100){

yb<-abs(rnorm(n, mean=mupred,sd=sqrt(varerror)))
model.b<-glmer(yb~Level*Time+(1+Time|ID),family=gaussian(link="log"), data=data)
mupred.star<-predict(model.b, newdata=newdata, type="response")
varerror.star<-(sigma(model.b)^2)
yf.star[b]<-rnorm(1, mean=mupred.star,sd=sqrt(varerror.star))

}

lower.bound<-quantile(yf.star, 0.1)
upper.bound<-quantile(yf.star, 0.9)
lower.bound
upper.bound


### Estimated quantiles method

newdata<-data.frame(ID=12, Level="15%", Time=20)
mupred<-predict(model.Gamma, newdata=newdata, type="response")
varerror<-(sigma(model.Gamma)^2)*mupred^2

a.hat<-(mupred^2)/varerror
s.hat<-varerror/mupred

lower<-qgamma(0.1, shape=a.hat,scale=s.hat)
upper<-qgamma(0.9, shape=a.hat,scale=s.hat)
lower
upper

## normal distribution 

newdata<-data.frame(ID=12, Level="15%", Time=20)
mupred<-predict(model.log, newdata=newdata, type="response")
varerror<-(sigma(model.log)^2)

lower<-qnorm(0.1, mean=mupred, sd=sqrt(varerror))
upper<-qnorm(0.9, mean=mupred, sd=sqrt(varerror))
lower
upper


############################
### Example 2.2

data<-read.table("ratescancer.txt", sep="\t", dec=".", header=TRUE)
attach(data)

interaction.plot(age, city, (cases/pop))

library(lme4)
model<-glmer(cases~offset(log(pop))+age+(1|city), data=data, family=poisson(link="log"))
summary(model)

model.nb<-glmer.nb(cases~offset(log(pop))+age+(1|city), data=data)
summary(model.nb)
getME(model.nb, "glmer.nb.theta")

predict(model, newdata=data, type="response")

### a)

newdata<-data.frame(city="Kolding", age="70-74", pop=535)

muhat<-predict(model, newdata=newdata, type="response", re.form=NA)
ratio.hat<-muhat/newdata$pop
ratio.hat

mupred<-predict(model, newdata=newdata, type="response")
ratio.pred<-mupred/newdata$pop
ratio.pred

ranef(model)
ratio.hat*exp(ranef(model)$city[3,1])

## plotting with whole data

muhat<-predict(model, newdata=data, type="response", re.form=NA)
ratio.hat<-muhat/data$pop
mupred<-predict(model, newdata=data, type="response")
ratio.pred<-mupred/data$pop

par(mfrow=c(2,2))
interaction.plot(age, city, (cases/pop))
interaction.plot(age, city, ratio.hat)
interaction.plot(age, city, ratio.pred)

### b)

model.H0<-glmer(cases~offset(log(pop))+(1|city), data=data, family=poisson(link="log"))
summary(model.H0)
anova(model.H0, model) # Not reliable

k1<-c(0,1,0,0,0,0)
k2<-c(0,0,1,0,0,0)
k3<-c(0,0,0,1,0,0)
k4<-c(0,0,0,0,1,0)
k5<-c(0,0,0,0,0,1)

K<-cbind(k1,k2,k3,k4,k5)
beta<-t(t(fixef(model)))

Wald<-t(beta)%*%K%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%beta
Wald
p.value<-1-pchisq(as.numeric(Wald), df=5)
p.value


### c)

### Parametric bootstrap method

n<-dim(data)[1]
newdata<-data.frame(city="Kolding", age="70-74", pop=535)
mupred.new<-predict(model, newdata=newdata, type="response")
ratio.new<-mupred.new/newdata$pop
ratio.new

mupred<-predict(model, newdata=data, type="response")
yf.star<-numeric()

for(b in 1:1000){

yb<-rpois(n, lambda=mupred)
data.b<-data.frame(yb,data)
model.b<-glmer(yb~offset(log(pop))+age+(1|city), data=data.b, family=poisson(link="log"))
mupred.star<-predict(model.b, newdata=newdata, type="response")
yf.star[b]<-rpois(1, lambda=mupred.star)

}

lower.bound<-quantile(yf.star, 0.1)/newdata$pop
upper.bound<-quantile(yf.star, 0.9)/newdata$pop
lower.bound
upper.bound


### Estimated quantiles method

newdata<-data.frame(city="Kolding", age="70-74", pop=535)
mupred.new<-predict(model, newdata=newdata, type="response")

lowerRatio.bound<-qpois(0.1, lambda=mupred.new)/newdata$pop
upperRatio.bound<-qpois(0.9, lambda=mupred.new)/newdata$pop
lowerRatio.bound
upperRatio.bound

#################


############################
### Example 2.3

data<-read.table("kneepain.txt", sep="\t", dec=".", header=TRUE)
attach(data)

plot(data$T,jitter(data$Y, amount=0.05))
id<-as.numeric(levels(factor(data$N)))

par(mfrow=c(2,2))
plot(data$T,jitter(data$Y, amount=0.05), type="n")
lines(data$T[data$N==id[1]],jitter(data$Y[data$N==id[1]], amount=0.05), col="red")
plot(data$T,jitter(data$Y, amount=0.05), type="n")
lines(data$T[data$N==id[2]],jitter(data$Y[data$N==id[2]], amount=0.05), col="red")
plot(data$T,jitter(data$Y, amount=0.05), type="n")
lines(data$T[data$N==id[4]],jitter(data$Y[data$N==id[4]], amount=0.05), col="blue")
plot(data$T,jitter(data$Y, amount=0.05), type="n")
lines(data$T[data$N==id[5]],jitter(data$Y[data$N==id[5]], amount=0.05), col="blue")

library(lme4)


model.1<-glmer(Y~T+factor(Th)+factor(Sex)+(1|N), data=data, family=binomial(link="logit"))
summary(model.1)

model.2<-glmer(Y~T+factor(Th)+factor(Sex)+(1+T|N), data=data, family=binomial(link="logit"))
summary(model.2)

model.22<-glmer(Y~T+factor(Th)+factor(Sex)+(1|N)+(T-1|N), data=data, family=binomial(link="logit"))
summary(model)

### a)

## M1

newdata<-data.frame(N=127, Th=2, Sex=0, T=11)
mupred<-predict(model.1, newdata=newdata, type="response", se.fit=TRUE)
mupred

newdata<-data.frame(N=1, Th=1, Sex=1, T=11)
mupred<-predict(model.1, newdata=newdata, type="response")
mupred

## M2

newdata<-data.frame(N=127, Th=2, Sex=0, T=11)
mupred<-predict(model.2, newdata=newdata, type="response", se.fit=TRUE)
mupred

newdata<-data.frame(N=1, Th=1, Sex=1, T=11)
mupred<-predict(model.2, newdata=newdata, type="response")
mupred

### b)

model.H0<-glmer(Y~T+factor(Sex)+(1|N), data=data, family=binomial(link="logit"))
summary(model.H0)
anova(model.H0,model.1)

model.H0<-glmer(Y~T+factor(Sex)+(1+T|N), data=data, family=binomial(link="logit"))
summary(model.H0)
anova(model.H0,model.2)


### c)

ST<-getME(model.2, "ST")$N
S<-diag(diag(ST))
T<-ST-S+diag(2)
F<-T%*%S%*%S%*%t(T)
cov.bi<-F
cov.bi
1.912347/sqrt(0.8643005*4.231250)

### d)

n<-dim(data)[1]
newdata<-data.frame(Th=1, Sex=1, T=3)
muhat.new<-predict(model.1, newdata=newdata, type="response", re.form=NA)
YS.pred<-100*muhat.new

mupred<-predict(model.1, newdata=data, type="response")

e.b<-numeric()

for(b in 1:1000){

yb<-numeric()
for(i in 1:n){

yb[i]<-sample(0:1,1,prob=c(1-mupred[i],mupred[i]))

}

model.b<-glmer(yb~T+factor(Th)+factor(Sex)+(1|N), data=data, family=binomial(link="logit"))
muhat.star<-predict(model.b, newdata=newdata, type="response", re.form=NA)
YS.predB<-100*muhat.star

yf.b<-sample(0:1,100,prob=c(1-muhat.star,muhat.star), replace=TRUE)

e.b[b]<-sum(yf.b)-YS.predB

}

var.error<-var(e.b)
var.error

z<-qnorm(c(0.1), lower.tail=FALSE)
lower.bound<-YS.pred-z*sqrt(var.error)
upper.bound<-YS.pred+z*sqrt(var.error)
lower.bound
upper.bound




