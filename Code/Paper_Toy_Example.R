set.seed(12345)
n=100
tm=rexp(n)
ev=factor(rbinom(n,2,0.3))
x=rnorm(n)
y=rbinom(n,1,0.5)
u=runif(n)
dset=cbind.data.frame(tm=tm,ev=ev,x=x,y=y,u=u)
fgm.res=crrf(Surv(tm,ev)~x+y+u,etype="1",dset,firth=TRUE,CI=TRUE)
fgm.res

fgm.test <- crrftest(Surv(tm,ev)~x+y+u,etype="1",dset,test=~x+y,firth=TRUE)
fgm.test

write.csv(dset, file="Y:/Anna/CRR/Code/ToyData.csv")
