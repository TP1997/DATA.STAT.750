setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("growthheight.txt", sep="\t", dec=".", head=TRUE)
attach(data)

# a)
Y = cbind(Y16,Y17,Y18)
model = lm(Y~factor(gender), data=data)
summary(model)

B = coef(model)
K = t(c(1,1))
L1 = c(1,-1,0)
L2 = c(0,1,-1)
L<-cbind(L1,L2)

YL = Y%*%L
model.H1 = lm(YL~factor(gender))
summary(model.H1)

library(MASS)
T = diag(2) - ginv(K)%*%K
X = model.matrix(model)
X.star = X%*%T

model.H0 = lm(YL~X.star-1)
summary(model.H0)
anova(model.H0,model.H1)


