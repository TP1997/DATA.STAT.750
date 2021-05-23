library(fda)
library(refund)
data(growth)
Y1 = growth$hgtf
Y2 = growth$hgtm
t = growth$age

par(mfrow=c(1,2))
matplot(t, Y1, pch=10, type='n', main="Height - Girls")
matlines(t, Y1, col="red")
matplot(t, Y2, pch=19, type='n', main="Height - Boys")
matlines(t, Y1, col="blue")

head(Y1[,1:5])
head(Y2[,1:5])

# a)
par(mfrow=c(1,1))
max(t)
min(t)
basis = create.bspline.basis(c(1,18), nbasis=15, norder=6)
plot(basis)
Y1.fd = smooth.basis(t, Y1, basis)
plot(Y1.fd)
title('Girl heights, nbasis=15, norder=6')

c = Y1.fd$fd$coefs
c

# b)
basis = create.bspline.basis(c(1,18), nbasis=20, norder=6)
plot(basis)
basisPar = fdPar(basis, Lfdobj=2, lambda=0.01)
Y1.fd = smooth.basis(t, Y1, basisPar)
plot(Y1.fd)
title('Girl heights, nbasis=20, norder=6, penalized')

# c)
Y2.fd = smooth.basis(t, Y2, basisPar)
plot(Y2.fd)
title('Boy heights, nbasis=20, norder=6, penalized')

Y1fd.mean = mean.fd(Y1.fd$fd)
Y2fd.mean = mean.fd(Y2.fd$fd)
plot(Y2fd.mean, type='n')
lines(Y2fd.mean, col='blue', lwd=2)
lines(Y1fd.mean, col='red', lwd=2)
title('Mean curves, blue=boys')

# Sample means at t=15:
# Boys:
eval.fd(15, Y2fd.mean)
# Girls:
eval.fd(15, Y1fd.mean)

# d)
time = seq(1,18,0.001)
Y1.dfit = eval.fd(time, Y1.fd$fd, 1)
plot(c(1,19), c(-5,30), type="n", main="Derivates")
for(j in 1:dim(Y1.dfit)[2]){
  lines(time, Y1.dfit[,j], col="orange", lwd=1)
}

# e)

time.f = seq(18,19,0.01)
basis.f = create.bspline.basis(c(1,19), nbasis=10, norder=3)
#basisPar.f = fdPar(basis.f, Lfdobj=2, lambda=0.01)
Y1f.fd = smooth.basis(t, Y1, basis.f)
plot(Y1f.fd)

Phi.f = eval.basis(time.f, basis.f)
C.hat = Y1f.fd$y2cMap%*%Y1
Pred = Phi.f%*%C.hat

y1f.pred = Pred[,1]
y1f.pred

plot(y1f.pred)
plot(c(1,19),c(0,180), type="n")
lines(Y1f.fd$fd[1], lwd=2)
lines(time.f, y1f.pred, col="brown", lwd=3)
title('Predicted height values for girl 1')

# f)
y1f = Y1[1:7,4]

mean.fd = mean.fd(Y1.fd$fd)
mu = eval.fd(t, mean.fd)
mu1 = mu[1:7]
mu2 = mu[-(1:7)]
ones = rbind(rep(1,dim(Y1)[2]))

V = Y1-(mu%*%ones)
Sigma = V%*%t(V)/(dim(Y1)[1]-1)

Sigma11 = Sigma[1:7,1:7]
Sigma12 = Sigma[1:7,-(1:7)]
Sigma21 = Sigma[-(1:7),1:7]
Sigma22 = Sigma[-(1:7),-(1:7)]

y2f.hat = mu2+Sigma21%*%solve(Sigma11)%*%(y1f-mu1)
y2f.hat

plot(t[-(1:7)], y2f.hat, type='l')
title('Predicted yf2 values')
