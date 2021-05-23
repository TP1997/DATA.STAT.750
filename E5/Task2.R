library(fda)
library(FDboost)
data(emotion)
Y = t(emotion$EMG)
X1 = t(emotion$EEG)
t = emotion$t
matplot(t, Y, type="n", xlab="t", ylab="Y")
matlines(t, Y, lty=1, col="red")
matplot(t, X1, type="n", xlab="t", ylab="X1")
matlines(t, X1, lty=1, col="blue")
x2 = emotion$power
x2
x3 = emotion$game_outcome
x3
x4 = emotion$control
x4

matplot(t,Y[,1:2],type="n")
matlines(t,Y[,1:2],col="red")

matplot(t,X1[,1:2],type="n")
matlines(t,X1[,1:2],col="blue")
# a)
basis.x1 = create.fourier.basis(c(1,384), nbasis=80)
X1.fd = smooth.basis(t, X1, basis.x1)
plot(X1.fd)
title('X1.fd')
x1.coefs = X1.fd$fd$coefs

#basis.y = create.bspline.basis(c(1,384), nbasis=110, norder=4)
basis.y = create.fourier.basis(c(1,384), nbasis=150)
Y.fd = smooth.basis(t, Y, basis.y)
plot(Y.fd)
title('Y.fd')
y.coefs = Y.fd$fd$coefs

# b)
X.list = vector("list",2)
X.list[[1]] = rep(1,dim(X1)[2])
X.list[[2]] = X1.fd$fd

beta0.basis = create.fourier.basis(c(1,384), nbasis=50)#, norder=4)
#beta0.basis = create.bspline.basis(c(1,384), nbasis=20, norder=4)
#beta0.basis = create.constant.basis(rangeval=c(1,384), names='const')
#beta1.basis = create.fourier.basis(c(1,384), nbasis=170)#names='const')
#beta1.basis = create.bspline.basis(c(1,384), nbasis=50, norder=4)
beta1.basis = create.constant.basis(rangeval=c(1,384), names='const')

beta.list = vector('list', 2)
beta.list[[1]] = beta0.basis
beta.list[[2]] = beta1.basis

model = fRegress(Y.fd$fd, X.list, beta.list)
Yhat.fd = model$yhatfdobj
plot(Y.fd, type='n')
lines(Y.fd, col='red')
lines(Yhat.fd, col='black', lwd=3)

beta0.fd = model$betaestlist[[1]]
beta1.fd = model$betaestlist[[2]]
plot(beta0.fd)
title('beta0')
plot(beta1.fd)
title('beta1')

# c)
X.list = vector('list',4)
X.list[[1]] = rep(1,dim(Y)[2])
X.list[[2]] = as.numeric(x2=='high') #high=1, low=0
X.list[[3]] = as.numeric(x3=='gain')#gain=1, loss=0
X.list[[4]] = as.numeric(x4=='high') #high=1, low=0
  
beta0.basis = create.fourier.basis(c(1,384), nbasis=50)
beta2.basis = create.constant.basis(rangeval=c(1,384), names='const')
beta3.basis = create.constant.basis(rangeval=c(1,384), names='const')
beta4.basis = create.constant.basis(rangeval=c(1,384), names='const')

beta.list = vector('list',4)
beta.list[[1]] = beta0.basis
beta.list[[2]] = beta2.basis
beta.list[[3]] = beta3.basis
beta.list[[4]] = beta4.basis

model = fRegress(Y.fd$fd, X.list, beta.list)
Yhat.fd = model$yhatfdobj
plot(Y.fd, type='n')
lines(Y.fd, col='red')
lines(Yhat.fd, col='black', lwd=3)

beta0.fd = model$betaestlist[[1]]
beta2.fd = model$betaestlist[[2]]
beta3.fd = model$betaestlist[[3]]
beta4.fd = model$betaestlist[[4]]

time = seq(1,384,0.1)
fit.1 = eval.fd(time, beta0.fd$fd)+eval.fd(time, beta2.fd$fd)+eval.fd(time, beta3.fd$fd)+eval.fd(time, beta4.fd$fd)
fit.2 = eval.fd(time, beta0.fd$fd)

plot(Y.fd, type='n')
lines(Y.fd, col='red')
plot(time, fit.1, col='black', type='l', lwd=2)
lines(time, fit.2, col='blue', lwd=2)
title('y estimates, blue=(low, loss, low)')

# d)
#basis.x1 = create.bspline.basis(c(1,384), nbasis=50, norder=4)#create.fourier.basis(c(1,384), nbasis=80)
basis.x1 = create.fourier.basis(c(1,384), nbasis=80)
X1.fd = smooth.basis(t, X1, basis.x1)

X.list = vector('list',5)
X.list[[1]] = rep(1,dim(Y)[2])
X.list[[2]] = X1.fd$fd
X.list[[3]] = as.numeric(x2=='high') #high=1, low=0
X.list[[4]] = as.numeric(x3=='gain')#gain=1, loss=0
X.list[[5]] = as.numeric(x4=='high') #high=1, low=0

beta0.basis = create.fourier.basis(c(1,384), nbasis=80)
beta1.basis = create.constant.basis(rangeval=c(1,384), names='const')#create.bspline.basis(c(1,384), nbasis=50, norder=4)
#beta1.basis = create.fourier.basis(c(1,384), nbasis=80)
beta2.basis = create.constant.basis(rangeval=c(1,384), names='const')
#beta2.basis = create.fourier.basis(c(1,384), nbasis=80)
beta3.basis = create.constant.basis(rangeval=c(1,384), names='const')
#beta3.basis = create.fourier.basis(c(1,384), nbasis=80)
beta4.basis = create.constant.basis(rangeval=c(1,384), names='const')
#beta4.basis = create.fourier.basis(c(1,384), nbasis=80)

beta.list = vector('list',5)
beta.list[[1]] = beta0.basis
beta.list[[2]] = beta1.basis
beta.list[[3]] = beta2.basis
beta.list[[4]] = beta3.basis
beta.list[[5]] = beta4.basis

model = fRegress(Y.fd$fd, X.list, beta.list)
Yhat.fd = model$yhatfdobj
plot(Yhat.fd)
