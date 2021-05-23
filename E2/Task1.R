setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("aids.txt", sep="\t", dec=".", head=TRUE)
attach(data)

plot(time, cd4, type="n")
for(i in 1:length(levels(factor(person)))){
  
  osa<-data[person==levels(factor(person))[i],]
  print(i)
  if(osa$drugs[1]==0){
    lines(osa$time, osa$cd4, lwd=2, col="green")
  }
  if(osa$drugs[1]==1){
    lines(osa$time, osa$cd4, lwd=2, col="blue")
  }
}

library(lme4)
# a)
# Select log link since cd4 is always positive (and other link functions caused some errors)
model.poi = glmer(cd4~time*factor(drugs)+age+(1+time|person), family=poisson(link='log'), data=data)
summary(model.poi)
model.gamma = glmer(cd4~time*factor(drugs)+age+(1+time|person), family=Gamma(link='log'), data=data)
summary(model.gamma)
model.ig = glmer(cd4~time*factor(drugs)+age+(1+time|person), family=inverse.gaussian(link='log'), data=data)
summary(model.ig)
model.gaus = glmer(cd4~time*factor(drugs)+age+(1+time|person), family=gaussian(link='log'), data=data)
summary(model.gaus)

par(mfrow=c(2,2))
plot(fitted(model.poi), residuals(model.poi, type="pearson"))
plot(fitted(model.gamma), residuals(model.gamma, type="pearson"))
plot(fitted(model.ig), residuals(model.ig, type="pearson"))
plot(fitted(model.gaus), residuals(model.gaus, type="pearson"))
par(mfrow=c(1,1))

par(mfrow=c(2,2))
qqnorm(residuals(model.poi, type='pearson'))
qqnorm(residuals(model.gamma, type='pearson'))
qqnorm(residuals(model.ig, type='pearson'))
qqnorm(residuals(model.gaus, type='pearson'))
par(mfrow=c(1,1))

sum(residuals(model.poi, type='response')^2)
sum(residuals(model.gamma, type='response')^2)
sum(residuals(model.ig, type='response')^2)
sum(residuals(model.gaus, type='response')^2)

sum(residuals(model.poi, type='pearson')^2)
sum(residuals(model.gamma, type='pearson')^2)
sum(residuals(model.ig, type='pearson')^2)
sum(residuals(model.gaus, type='pearson')^2)

sum(residuals(model.poi, type='deviance')^2)
sum(residuals(model.gamma, type='deviance')^2)
sum(residuals(model.ig, type='deviance')^2)
sum(residuals(model.gaus, type='deviance')^2)


new.data = data[,c(2,3,4,5)]
gt.y = data[,1]
poi.predy = predict(model.poi, newdata=new.data, type="response")
gamma.predy = predict(model.gamma, newdata=new.data, type="response")
ig.predy = predict(model.ig, newdata=new.data, type="response")
gaus.predy = predict(model.gaus, newdata=new.data, type="response")

par(mfrow=c(2,2))
plot(gt.y, fitted(model.poi))
abline(0,1)
title('Poisson model')
plot(gt.y, fitted(model.gamma))
abline(0,1)
title('Gamma model')
plot(gt.y, fitted(model.ig))
abline(0,1)
title('IG model')
plot(gt.y, fitted(model.gaus))
abline(0,1)
title('Gaussian model')

par(mfrow=c(3,2))
hist(gt.y)
hist(poi.predy)
hist(gamma.predy)
hist(ig.predy)
hist(gaus.predy)
par(mfrow=c(1,1))

# Shapiro test
shapiro.test(residuals(model.poi, type="pearson"))
shapiro.test(residuals(model.gamma, type="pearson"))
shapiro.test(residuals(model.ig, type="pearson"))
shapiro.test(residuals(model.gaus, type="pearson"))

# Select Gamma model
model.gamma = glmer(cd4~time*factor(drugs)+age+(1+time|person), family=Gamma(link='log'), data=data)
summary(model.gamma)

# b)
new.data = data.frame(time=2.02, drugs=1, age=13.72)
mu.hat = predict(model.gamma, newdata=new.data, re.form=NA, type='response')
mu.hat

# c)
### Parametric bootstrap method
new.data = data.frame(time=2.02, drugs=1, age=13.72, person=10396)
mu.pred = predict(model.gamma, newdata=new.data, type='response')
var.error = (sigma(model.gamma)^2)*mu.pred^2
a.hat = (mu.pred^2)/var.error
s.hat = var.error/mu.pred

yf.star = numeric()
n = dim(data)[1]
for(b in 1:1000){
  yb = rgamma(n, shape=a.hat, scale=s.hat)
  model.b = glmer(yb~time*factor(drugs)+age+(1+time|person), family=Gamma(link='log'), data=data)
  mu.pred.star = predict(model.b, newdata=new.data, type='response')
  var.error.star = (sigma(model.b)^2)*mu.pred.star^2
  a.star = (mu.pred.star^2)/var.error.star
  s.star = var.error.star/mu.pred.star
  yf.star[b] = rgamma(1, shape=a.star, scale=s.star)
}

lb = quantile(yf.star, 0.1)
ub = quantile(yf.star, 0.9)
lb
ub
