setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("growthheight.txt", sep="\t", dec=".", head=TRUE)
attach(data)

Y = cbind(Y10,Y11,Y12,Y13,Y14,Y15,Y16,Y17,Y18)
t = 10:18
T = cbind(rep(1,9), t, t^2, t^3)
Y.star = Y%*%T%*%solve(t(T)%*%T)

# a)
model = lm(Y.star~factor(gender)-1, data=data)
summary(model)
Theta = coef(model)
Theta

# b)
t10 = 19
new.t = cbind(1,t10,t10^2,t10^3)
mu.10 = new.t%*%t(Theta)%*%t(cbind(0,1))
T%*%t(Theta)%*%t(cbind(0,1))
mu.10

# c)
BP = T%*%t(Theta)
BP
plot(c(10,18), c(130,190), type="n", axes=F, xlab="", ylab="")
points(10:18, BP[,1], col="blue", pch=19)
lines(10:18, BP[,1], col="blue")
points(10:18, BP[,2], col="red", pch=19)
lines(10:18, BP[,2], col="red")
axis(1, 10:18)
axis(2, round(seq(130,190, length.out=19), 0))
legend(10,190,legend=c("girl","boy"), col=c("red","blue"), lty=1:1, cex=0.8)

# d)
Y.star = Y%*%T%*%solve(t(T)%*%T)
model.H0 = lm(Y.star~1, data=data)
model.H1 = lm(Y.star~factor(gender)-1, data=data)

anova(model.H0, model.H1)

# e)
X = model.matrix(model)
T1 = T[,1:2]
Z1 = T1%x%X
Z = T%x%X
Y.vec = as.vector(Y)
model.H0 = lm(Y.vec~Z1-1, data=data)
model.H1 = lm(Y.vec~Z-1, data=data)

mean(residuals(model.H0)^2)
mean(residuals(model.H1)^2)

AIC(model.H0)
AIC(model.H1)

# f)
new.data = data.frame(gender='girl', Y10=148.7, Y11=156.6, Y12=164.7)

model = lm(Y.star~factor(gender)-1, data=data)
Theta = coef(model)
T = cbind(rep(1,9), t, t^2, t^3)
T1 = T[1:3,]
T2 = T[-(1:3),]
X = model.matrix(model)
E = Y-X%*%Theta%*%t(T)
Sigma = (t(E)%*%E)/model$df.residual
Sigma11 = Sigma[1:3,1:3]
Sigma12 = Sigma[1:3,4:9]
Sigma21 = Sigma[4:9,1:3]
Sigma22 = Sigma[4:9,4:9]

yf1 = t(cbind(148.7, 156.6, 164.4))
xf = t(cbind(0,1))
BP = T2%*%t(Theta)%*%xf + Sigma21%*%solve(Sigma11)%*%(yf1-T1%*%t(Theta)%*%xf)

cov.error = Sigma22-Sigma21%*%solve(Sigma11)%*%Sigma12
lb = BP-qnorm(0.9)*diag(cov.error)^0.5
ub = BP+qnorm(0.9)*diag(cov.error)^0.5
lb
BP
ub