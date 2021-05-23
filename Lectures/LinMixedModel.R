### Statistical Modeling 2
### Linear Mixed Models
### 2021

############################
### Example 1.1

data<-read.table("pollutiondata.txt",header=TRUE,sep="\t")
data<-as.matrix(data)
edata<-embed(data,2)
data<-data.frame(edata[,c(1,5,6)])
colnames(data)<-c("Mortality","TemperatureT1","ParticulatesT1")

training<-data[-507,]
newdata<-data[507,]

plot(training)

modelLM<-lm(Mortality~TemperatureT1+ParticulatesT1,data=training)
summary(modelLM)

### We are using residuals to investigate the ARMA process 

e<-residuals(modelLM)

plot(e,type="n")
lines(e)

acf(e)
pacf(e)

ar(e,demean=FALSE)

## a) 

library(nlme)

modelGLS<-gls(Mortality~TemperatureT1+ParticulatesT1,data=training,correlation=corARMA(p=4,q=0))
summary(modelGLS)

paraARMA<-coef(modelGLS$model$corStruct,unconstrained=FALSE)
paraARMA

## b) 

modelH1<-gls(Mortality~TemperatureT1+ParticulatesT1,data=training,correlation=corARMA(p=4,q=0), method="ML")
summary(modelH1)

modelH0<-gls(Mortality~TemperatureT1,data=training,correlation=corARMA(p=4,q=0), method="ML")
summary(modelH0)

anova(modelGLS)
anova(modelH1)
anova(modelH0, modelH1)

## c)

nn<-dim(data)[1]
T<-data.frame(t=1:nn)
cs<-corARMA(paraARMA,form = ~ 1,p=4,q=0)
cs<-Initialize(cs, data = T)

Sigma<-corMatrix(cs)
V<-Sigma[1:(nn-1),1:(nn-1)]
w<-Sigma[1:(nn-1),nn]

blup<-predict(modelGLS,newdata=newdata)+t(w)%*%solve(V)%*%residuals(modelGLS)
blup

xf<-t(t(c(1,73.33,57.58)))
beta<-coef(modelGLS)
X<-model.matrix(~TemperatureT1+ParticulatesT1, data=training)
y<-training$Mortality

blup<-t(xf)%*%beta+t(w)%*%solve(V)%*%(y-X%*%beta)
blup

## d)

sigma2<-summary(modelGLS)$sigma^2
X<-model.matrix(~TemperatureT1+ParticulatesT1, data=training)
vf<-Sigma[nn,nn]
xf<-t(t(c(1,73.33,57.58)))
var.error<-sigma2*(vf-t(w)%*%solve(V)%*%w+(t(xf)-t(w)%*%solve(V)%*%X)%*%solve(t(X)%*%solve(V)%*%X)%*%(xf-t(X)%*%solve(V)%*%w))

lowerbound<-blup-qnorm(0.9)*sqrt(var.error)
lowerbound

upperbound<-blup+qnorm(0.9)*sqrt(var.error)
upperbound


############################
### Example 1.4


data<-read.table("saltgrass.txt", sep="\t", dec=".", header=TRUE)
attach(data)

interaction.plot(Salinity, Location,Ashcontent)

library(lme4)

model<-lmer(Ashcontent~factor(Salinity)+(1|Location), data=data)
summary(model)
fixef(model)
ranef(model)

## a)

VarCorr(model)
print(VarCorr(model), comp="Variance")
as.numeric(VarCorr(model))

## b)

newdata<-expand.grid(Location=1:16, Salinity=1:4)
mu.hat<-predict(model, newdata=newdata,re.form=NA)
interaction.plot(newdata$Salinity, newdata$Location,mu.hat)

mu.pred<-predict(model, newdata=newdata)
interaction.plot(newdata$Salinity, newdata$Location,mu.pred)

#### answer: 4=50 dS/m

## c)

ranef(model)
ranef(model)$Location[[1]][12]


## d)

model.H1<-lmer(Ashcontent~factor(Salinity)+(1|Location), data=data, REML=FALSE)
summary(model.H1)
model.H0<-lmer(Ashcontent~1+(1|Location), data=data, REML=FALSE)
summary(model.H0)
anova(model.H0, model.H1)
anova(model.H0, model.H1)$Chisq[2]

## Testing Location

library(lmerTest)
ranova(model)

## e)

newdata<-expand.grid(Location=1:16, Salinity=1:4)
mu.pred<-predict(model, newdata=newdata)
mu.pred
data.frame(newdata, mu.pred)[21,]

## f)

beta<-t(t(fixef(model)))
b<-as.matrix(ranef(model)$Location)

X<-getME(model, "X")
xf<-t(t(c(1,1,0,0)))

Z<-getME(model, "Z")
zf<-t(t(rep(0,16)))
zf[5]<-1
zf<-as.matrix(zf)

blup<-t(xf)%*%beta+t(zf)%*%b

newdata<-expand.grid(Location=5, Salinity=2)
predict(model, newdata=newdata)

VarCorr(model)
attr(VarCorr(model)$Location,"stddev")^2

sigma2<-sigma(model)^2
sigma2

ST<-getME(model, "ST")$Location
S<-diag(ST)
T<-ST-S+diag(1)
F<-T%*%S%*%S%*%t(T)
cov.bh<-sigma2*F
cov.bh
G<-kronecker(diag(16), F)
cov.b<-sigma2*G 
cov.b

V<-diag(dim(X)[1])+Z%*%G%*%t(Z)
cov.y<-sigma2*V
cov.y

y<-getME(model, "y")

beta.hat<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%y
beta.hat
b.hat<-G%*%t(Z)%*%solve(V)%*%(y-X%*%beta.hat)
b.hat

### Prediction interval - Asymptotic method

w<-Z%*%G%*%zf
vf<-1+t(zf)%*%G%*%zf

var.error<-sigma2*(vf-t(w)%*%solve(V)%*%w+(t(xf)-t(w)%*%solve(V)%*%X)%*%solve(t(X)%*%solve(V)%*%X)%*%(xf-t(X)%*%solve(V)%*%w))

lowerbound<-blup-qnorm(0.9)*sqrt(var.error)
lowerbound

upperbound<-blup+qnorm(0.9)*sqrt(var.error)
upperbound

### Prediction interval - Quantiles method

lowerbound<-qnorm(0.1, mean=blup, sd=sqrt(sigma2))
lowerbound

upperbound<-qnorm(0.9, mean=blup, sd=sqrt(sigma2))
upperbound


############################
### Example 1.5

data<-read.table("fertilizer.txt", sep="\t", dec=".", header=TRUE)
attach(data)

interaction.plot(week, plant,root)

model<-lmer(root~fertilizer*week+(1+week|plant), data=data)
summary(model)

model<-lmer(root~fertilizer*week+(1|plant)+(week-1|plant), data=data)
summary(model)


## a)

ranef(model)
ranef(model)$plant[[2]]

## b)

VarCorr(model)

sigma2<-sigma(model)^2
sigma2

F<-diag(c(as.numeric(getME(model, "ST")[1])^2,as.numeric(getME(model, "ST")[2])^2))
cov.bi<-sigma2*F
cov.bi

## c)

model.H0<-lmer(root~week+(1|plant)+(week-1|plant), data=data, REML=FALSE)
model.H1<-lmer(root~fertilizer*week+(1|plant)+(week-1|plant), data=data, REML=FALSE)
anova(model.H0, model.H1)
anova(model.H0, model.H1)$Chisq[2]


## d)

model<-lmer(root~fertilizer*week+(1|plant)+(week-1|plant), data=data)
summary(model)

beta<-t(t(fixef(model)))
b<-as.numeric(stack(ranef(model)$plant)[,1])

X<-getME(model, "X")
xf<-t(t(c(1,0,11,0)))

Z<-getME(model, "Z")
zf<-t(t(rep(0,24)))
zf[1]<-1
zf[13]<-11
zf<-as.matrix(zf)

blup<-t(xf)%*%beta+t(zf)%*%b

sigma2<-sigma(model)^2
sigma2

G1<-kronecker(diag(12), F[1,1])
G2<-kronecker(diag(12), F[2,2])
zero<-matrix(0, ncol=12, nrow=12)
G<-rbind(cbind(G1,zero),cbind(zero,G2))
cov.b<-sigma2*G 
cov.b

V<-diag(dim(X)[1])+Z%*%G%*%t(Z)
cov.y<-sigma2*V
cov.y

### Prediction interval - Asymptotic method

w<-Z%*%G%*%zf
vf<-1+t(zf)%*%G%*%zf

var.error<-sigma2*(vf-t(w)%*%solve(V)%*%w+(t(xf)-t(w)%*%solve(V)%*%X)%*%solve(t(X)%*%solve(V)%*%X)%*%(xf-t(X)%*%solve(V)%*%w))

lowerbound<-blup-qnorm(0.9)*sqrt(var.error)
lowerbound

upperbound<-blup+qnorm(0.9)*sqrt(var.error)
upperbound


############################
### Example 1.6


library(faraway)
data(psid)
psid
attach(psid)

interaction.plot(year, person, income)

model<-lmer(income~year+sex+(year|person), data=psid)
summary(model)


## a)

VarCorr(model)

sigma2<-sigma(model)^2
sigma2

ST<-getME(model, "ST")$person
S<-diag(diag(ST))
T<-ST-S+diag(2)
F<-T%*%S%*%S%*%t(T)
cov.bi<-sigma2*F
cov.bi

## b)

X<-getME(model, "X")
xf<-t(t(c(1,84,1)))

Z<-getME(model, "Z")
zf<-t(t(rep(0,170)))
zf[1]<-1
zf[2]<-84
zf<-as.matrix(zf)

beta<-t(t(fixef(model)))
b<-as.numeric(stack(data.frame(t(ranef(model)$person)))[,1])

G<-kronecker(diag(length(levels(factor(person)))), F)
cov.b<-sigma2*G 
cov.b

V<-diag(dim(X)[1])+Z%*%G%*%t(Z)
cov.y<-sigma2*V
cov.y

blup<-t(xf)%*%beta+t(zf)%*%b
blup

### Prediction interval - Asymptotic method

w<-Z%*%G%*%zf
vf<-1+t(zf)%*%G%*%zf

var.error<-sigma2*(vf-t(w)%*%solve(V)%*%w+(t(xf)-t(w)%*%solve(V)%*%X)%*%solve(t(X)%*%solve(V)%*%X)%*%(xf-t(X)%*%solve(V)%*%w))

lowerbound<-blup-qnorm(0.9)*sqrt(var.error)
lowerbound

upperbound<-blup+qnorm(0.9)*sqrt(var.error)
upperbound




