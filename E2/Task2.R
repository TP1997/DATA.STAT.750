setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("locust.txt", sep="\t", dec=".", head=TRUE)
attach(data)

library(lme4)
# a)
model.bin = glmer(move~time+sex+feed+(1+time|id), family=binomial(link="logit"), data=data)
summary(model.bin)
model.bin2 = glmer(move~time+sex+feed+(1|id)+(time-1|id), family=binomial(link="logit"), data=data)
summary(model.bin2)

beta1.hat = summary(model.bin)$coef[2]
beta1.hat

# b)
new.data = data.frame(id=24, sex=0, time=1.35, feed=0)
mu.pred = predict(model.bin, newdata=new.data, type='response')
mu.pred

# c)
model.H0 = glmer(move~time+sex+(1+time|id), family=binomial(link="logit"), data=data)
summary(model.H0)
anova(model.H0, model.bin)
p.val = anova(model.H0, model.bin)$"Pr(>Chisq)"[2]
p.val

# d)
ST = getME(model.bin, 'ST')$id
S = diag(diag(ST))
T = ST-S+diag(2)
F = T%*%S%*%S%*%t(T)
cov.bi = F
cov.bi
# -0.64*sqrt(1.299)*sqrt(2.385)

# e)
new.data = data.frame(sex=0, time=1.35, feed=0)
mu.hat.new = predict(model.bin, newdata=new.data, type="response", re.form=NA)
YS.pred = 100*mu.hat.new

mu.pred = predict(model.bin, newdata=data, type="response")

e.b = numeric()
n = dim(data)[1]
for(b in 1:200){
  yb = numeric()
  for(i in 1:n){
    yb[i]<-sample(0:1, 1, prob=c(1-mu.pred[i],mu.pred[i]))
  }
  model.b = glmer(yb~time+sex+feed+(1+time|id), family=binomial(link="logit"), data=data)
  mu.hat.star = predict(model.b, newdata=new.data, type="response", re.form=NA)
  YS.predB = 100*mu.hat.star
  
  yf.b = sample(0:1, 100, prob=c(1-mu.hat.star,mu.hat.star), replace=TRUE)
  
  e.b[b] = sum(yf.b) - YS.predB
}

var.error = var(e.b)
var.error

z = qnorm(c(0.1), lower.tail=FALSE)
lb = YS.pred-z*sqrt(var.error)
ub = YS.pred+z*sqrt(var.error)
lb
ub