setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("cmort.txt", sep="\t", dec=".", head=TRUE)
attach(data)
head(data)

# a)
model.M = lm(y~t+I(t^2)+I(t^3)+I(cos(2*pi*t/52))+I(sin(2*pi*t/52)), data=data)
summary(model.M)

library(nlme)
model.M = gls(y~t+I(t^2)+I(t^3)+I(cos(2*pi*t/52))+I(sin(2*pi*t/52)), 
                data=data,
                correlation=corARMA(p=1,q=0))
summary(model.M)

arma.phi = coef(model.M$model$corStruct,unconstrained=FALSE)
arma.phi

cat('Restricted MLE for phi =', arma.phi)

# b)
model.H0 = gls(y~t+I(t^2)+I(t^3), data=data, correlation=corARMA(p=1,q=0))
model.H1 = model.M
test = anova(model.H0, model.H1)
test$L.Ratio

cat('Likelihood ratio test statistic =', test$L.Ratio[2],
    '\np-value =', test$`p-value`[2],
    '\nModel H1 should be selected')

# c)
arma.par = coef(model.M$model$corStruct,unconstrained=FALSE)
nn = dim(data)[1]
T = data.frame(t=1:nn)
cs = corARMA(arma.par, form = ~1, p=1, q=0)
cs = Initialize(cs, data=T)

Sigma = corMatrix(cs)
V = Sigma[1:(nn),1:(nn)]
w = Sigma[1:(nn),nn]

new.data = data.frame('t'=509)
yf.bulp = predict(model.M, newdata=new.data)+t(w)%*%solve(V)%*%residuals(model.M)

cat('Empirical blup for new obs. y =', yf.bulp)

# d)
sigma2 = summary(model.M)$sigma^2
X = model.matrix(~t+I(t^2)+I(t^3)+I(cos(2*pi*t/52))+I(sin(2*pi*t/52)),data=data)
vf = Sigma[nn,nn]
xf = t(t(c(1,509,509^2,509^3,cos(2*pi*509/52),sin(2*pi*509/52))))
var.error = sigma2*(vf-t(w)%*%solve(V)%*%w+(t(xf)-t(w)%*%solve(V)%*%X)%*%solve(t(X)%*%solve(V)%*%X)%*%(xf-t(X)%*%solve(V)%*%w))




