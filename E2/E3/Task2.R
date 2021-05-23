setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("excretory.txt", sep="\t", dec=".", head=TRUE)
attach(data)

# a)
Y = cbind(creatinine,chloride,chlorine)
model = lm(Y~gravity+factor(obesity), data=data)
summary(model)

new.data = data.frame(gravity=35, obesity='Group III')
BP = predict(model, newdata=new.data)
BP

E = residuals(model)
Sigma = (t(E)%*%E) / model$df.residual
X = model.matrix(model)
# Levels: Group I Group II Group III Group IV
# Levels: 1 2 3 4
xf = t(cbind(1, 35, 0, 1, 0))
cov.error = Sigma*as.numeric(1+t(xf)%*%solve(t(X)%*%X)%*%xf)
lb = BP-qnorm(0.9)*diag(cov.error)^0.5
ub = BP+qnorm(0.9)*diag(cov.error)^0.5

lb
BP
ub

# b)
model.H0 = lm(cbind(creatinine,chloride)~factor(obesity), data=data)
summary(model.H0)
model.H1 = lm(cbind(creatinine,chloride)~gravity+factor(obesity), data=data)
summary(model.H1)

anova(model.H0, model.H1)

# c)
model.H0 = lm(cbind(creatinine,chloride)~chlorine, data=data)
summary(model.H0)
model.H1 = lm(cbind(creatinine,chloride)~factor(obesity)+chlorine, data=data)
summary(model.H1)

anova(model.H0, model.H1)
