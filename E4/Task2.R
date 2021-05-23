setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("stageforest.txt", sep="\t", dec=".", head=TRUE)
attach(data)

library(lme4)

# a)
Y = cbind(dbhib.cm, height.m)
T = cbind(1,Age,Age^2)
id.cnt = cumsum(as.numeric(table(Tree.ID)))
Ys = as.vector(Y[1:5,])
Xs = diag(2)%x%T[1:5,]
Zs = list(diag(2)%x%T[1:5,1:2])
for (i in 1:(length(id.cnt)-1)) {
  s = id.cnt[i]+1
  e = id.cnt[i+1]
  print(s)
  Ys = append(Ys, as.vector(Y[s:e,]))
  Xs = rbind(Xs, diag(2)%x%T[s:e,])
  Zs = append(Zs, diag(2)%x%T[s:e,1:2])
}
Zs = bdiag(Zs)
dd
model = lmer(Ys~-1+Xs+(1+rep(Age,2)|rep(Tree.ID,2)), data=data)
















Y = cbind(dbhib.cm, height.m)
y = as.vector(Y)
X = cbind(1,Age,Age^2)
x = diag(3)%x%X #as.vector(X)
Z = cbind(Age)
z = as.vector(Z)
model = lmer(y~-1+x+(1+z|Tree.ID), data=data)


T = cbind(1,Age,Age^2)
Y1 = Y[:5,]
y2s = as.vector(Y1)
X1 = T[1:5,]
X1s = diag(2)%x%X1
Z1 = T[1:5,1:2]
Z1s = diag(2)%x%Z1

Y2 = Y[6:15,]
y2s = as.vector(Y2)
X2 = T[6:15,]
X2s = diag(2)%x%X2
Z2 = T[6:15,1:2]
Z2s = diag(2)%x%Z2

