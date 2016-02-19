### demonstrate conservative nature of posterior p-values when
# conditioning on observed data and parameters

n.sims=10000
n.pois=10
mcmc.iter=1000
mh.lambda=2.2

mcmc.pois<-function(Counts,mcmc.iter){
  MCMC=rep(0,mcmc.iter)
  MCMC[1]=mean(Counts)
  old.log.post=sum(dpois(Counts,MCMC[1],log=1))#-0.5*log(MCMC[1])#Jeffreys prior
  Accept=rep(0,mcmc.iter)
  for(iiter in 2:mcmc.iter){
    prop=MCMC[iiter-1]+rnorm(1,0,mh.lambda)
    MCMC[iiter]=MCMC[iiter-1]
    if(prop>0){
      new.log.post=sum(dpois(Counts,prop,log=1))#-0.5*log(prop) #Jeffreys prior
      if(runif(1)<exp(new.log.post-old.log.post)){
        old.log.post=new.log.post
        MCMC[iiter]=prop
        Accept=Accept+1
      }
    }
  }
  MCMC=list(MCMC=MCMC,Accept=Accept)
}

chisq<-function(Data,mu){
  sum((Data-mu)^2)/mu
}

Pval1=Pval2=Pval3=rep(0,n.sims)
Chisq1=Chisq2=rep(0,mcmc.iter)
for(isim in 1:n.sims){
  lam.exp=runif(1,1,10)
  Counts=rpois(n.pois,lam.exp)
  Res=mcmc.pois(Counts,mcmc.iter)
  #calculate posterior p-value using different discrepancy functions
  for(idata in 1:mcmc.iter){
    Cur.data1=rpois(n.pois,Res$MCMC[idata])
    Chisq1[idata]=sum((Cur.data1-Res$MCMC[idata])^2/Res$MCMC[idata]) #dividing by var(y|theta)=theta
    Chisq2[idata]=sum((Counts-Res$MCMC[idata])^2/Res$MCMC[idata])
  }
  Pval1[isim]=sum(Chisq2>Chisq1)/mcmc.iter
  lam.exp.est=mean(Res$MCMC)
  Cur.data2=matrix(rpois(mcmc.iter*n.pois,lam.exp.est),mcmc.iter,n.pois)
  Chisq1=apply(Cur.data2,1,"chisq",mu=lam.exp.est)
  chisq2=sum((Counts-lam.exp.est)^2)/lam.exp.est
  Pval2[isim]=sum(chisq2>Chisq1)/mcmc.iter
}

pdf(file="Pval_hist.pdf")
par(cex.lab=1.3)
hist(Pval1,breaks=20,xlab="P-value",ylab="Empirical density",main='',freq=FALSE)
dev.off()