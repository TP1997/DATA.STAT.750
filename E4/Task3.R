library(geepack)
data(ohio)
head(ohio)
tail(ohio)

model = geeglm(resp~smoke, family=binomial(link='logit'), id=id, data=ohio, corstr="unstructured")
summary(model)

new.data = data.frame(id=536, age=1.25, smoke=1)
mu.hat = predict(model, new.data, type='response')
mu.hat
