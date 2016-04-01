# PIT application - simple poisson example with 1 parameter

#1) generate some data
Y=rpois(100,5)
U=rep(0,100)
sd=1
lam=10
n=99

for(i in 1:100){  #loop over leave-one-out
  Cur.Y=Y[-i]
  sum.y=sum(Cur.Y)
  for(iiter in 1:1000){  #MCMC
    prop=lam+rnorm(1,0,sd)
    MH=exp((log(prop)-log(lam))*sum.y+n*(lam-prop))
    if(runif(1)<MH)lam=prop
    temp.pred=rpois(1,lam)
    U[i]=U[i]+(temp.pred<Y[i])+0.5*(temp.pred==Y[i])
  }
}

U=U/1000
cat(mean(U))
hist(U)