setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("growthheight.txt", sep="\t", dec=".", head=TRUE)
attach(data)

# a)
Y = cbind(Y10,Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18)
model = lm(Y~factor(gender), data=data)
summary(model)

B.hat = coef(model)
B.hat

# b)
E = residuals(model)
Sigma.hat = (t(E)%*%E) / model$df.residual
Sigma.hat

# c)
new.data = data.frame(gender = c('girl', 'boy'))
BP = predict(model, newdata=new.data)
BP

plot(c(10,18), c(130,190), type="n", axes=F, xlab="", ylab="")
points(10:18, BP[1,], col="red", pch=19)
lines(10:18, BP[1,], col="red")
points(10:18, BP[2,], col="blue", pch=19)
lines(10:18, BP[2,], col="blue")
axis(1, 10:18)
axis(2, round(seq(130,190, length.out=19), 0))
legend(10,190,legend=c("girl","boy"), col=c("red","blue"), lty=1:1, cex=0.8)

# d)
model.H0 = lm(Y~1, data=data)
summary(model.H0)

anova(model.H0, model)

# e)
new.data = data.frame(gender="girl")
xf = t(cbind(1, 1))
yf1 = c(148.7, 156.6, 164.7)
Bt = t(coef(model))
B1t = Bt[1:3,]
B2t = Bt[4:9,]
Sigma11 = Sigma.hat[1:3,1:3]
Sigma12 = Sigma.hat[1:3,4:9]
Sigma21 = Sigma.hat[4:9,1:3]
Sigma22 = Sigma.hat[4:9,4:9]

BP.yf2 = B2t%*%xf + Sigma21%*%solve(Sigma11)%*%(yf1 - B1t%*%xf)

cov.error = Sigma22 - Sigma21%*%solve(Sigma11)%*%Sigma12
lb = BP.yf2 - qnorm(0.9)*diag(cov.error)^0.5
ub = BP.yf2 + qnorm(0.9)*diag(cov.error)^0.5
lb
BP.yf2
ub

