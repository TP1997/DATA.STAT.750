setwd('/home/tuomas/R/Projects/DATA.STAT.750/datasets')
data = read.table("stageforest.txt", sep="\t", dec=".", head=TRUE)
attach(data)
head(data)


# a)
library(lme4)
model = lmer(height.m~dbhib.cm+I(dbhib.cm^2)+(1|Forest.ID), data=data)
summary(model)
sigma2.z = as.numeric(VarCorr(model))
cat('RMLE for sigma2.z =', sigma2.z)

# b)
fnames=c('Clark Fork','Clearwater','Coeur dAlene','Kaniksu','Nez Perce','Payette','St. Joe','Umatilla','Wallowa')
reffects = ranef(model)$Forest.ID
largest = which.max(reffects[,1])
cat('Forest area with largest random effect:\n',
    fnames[largest],
    '\n',reffects[largest,])

# c)
model.c = lmer(height.m~dbhib.cm+I(dbhib.cm^2)+(1+dbhib.cm|Tree.ID)+(1|Forest.ID), data=data)
summary(model.c)
#model.c = lmer(height.m~dbhib.cm+I(dbhib.cm^2)+(1|Tree.ID)+(dbhib.cm-1|Tree.ID)+(1|Forest.ID), data=data)
#summary(model.c)
fixef(model.c)
ranef(model.c)

new.data = expand.grid(Tree.ID=67, Forest.ID='Wallowa', dbhib.cm=35.5)
mu.pred = predict(model.c, newdata=new.data)

cat('ML prediction for mu =', mu.pred)

