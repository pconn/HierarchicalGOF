##### Run spatial negative binomial regression model using JAGS

#1) simulate some data
setwd('c:/users/paul.conn/git/hierarchicalGOF/hierarchicalGOF/inst')
library(MASS)
library(fields)
library(mvtnorm)

n_Y=250
Knots=expand.grid(c(0:5),c(0:5))
n_k=nrow(Knots)
Coords=matrix(runif(n_Y*2,1,4),n_Y,2)
theta_true=1 
tau_true=1
Dyy=rdist(Coords,Coords)
Dyk=rdist(Coords,Knots)
Dkk=rdist(Knots,Knots)
Cov_kk=Exp.cov(Knots,theta=theta_true,distMat=Dkk)/tau_true
Cov_yk=Exp.cov(Coords,Knots,theta=theta_true,distMat=Dyk)/tau_true
X_RE = Cov_yk %*% solve(Cov_kk)
Eta_k_true=matrix(rmvnorm(1,rep(0,n_k),Cov_kk),ncol=1)
Eta_true=X_RE%*%Eta_k_true
Sim_lambda=rgamma(n_Y,100,100/exp(3+Eta_true))
Y=rpois(n_Y,Sim_lambda)
hist(Y)

#2) jags model
neg.bin.model <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Lambda[i] ~ dgamma(psi,Par2[i])
    Par2[i] <- psi/Mu[i]
    log(Mu[i]) <- beta0 + Eta_y[i]
  }
  Eta_y[1:n_Y] <- X_RE[1:n_Y,1:n_k,theta]%*%Eta_k[1:n_k]
  Eta_k[1:n_k] ~ dmnorm(Zero_k[1:n_k],Sigma_inv[1:n_k,1:n_k])
  Sigma_inv[1:n_k,1:n_k] <- tau*Q_k[1:n_k,1:n_k,theta]
  # Priors
  psi ~ dgamma(1,0.01)
  beta0 ~ dnorm(0,0.01)
  tau ~ dgamma(1.0,0.01)
  theta ~ dcat(P)
  # Functions of interest:
  # Posterior predictions
  #for(i in 1:n_Y){
  #  Y_sim[i] ~ dpois(Lambda[i])
  #  Y_sim2[i] ~ dpois(Lambda_sim[i])
  #  Lambda_sim[i] ~ dpois 
  #}
}


library(rjags)
library(R2jags)
set.factory("bugs::Conjugate",FALSE,type="sampler")
Zero_k = rep(0,n_k)
max_range=3
n_incr=10
incr=max_range/n_incr
Cor_k=Q_k=array(0,dim=c(n_k,n_k,n_incr))
X_RE=array(0,dim=c(n_Y,n_k,n_incr))
for(i in 1:n_incr){
  Cor_k[,1:n_k,i]=Exp.cov(Knots,theta=incr*i,distMat=Dkk)
  Q_k[,1:n_k,i]=solve(Cor_k[,1:n_k,i])
  X_RE[,1:n_k,i]=Exp.cov(Coords,Knots,theta=incr*i,distMat=Dyk)%*%Q_k[,1:n_k,i]
}
P=rep(1/n_incr,n_incr)
jags_data = list("n_Y","Y","Zero_k","Q_k","X_RE","n_k","P")
jags_params = c("beta0","psi","Lambda","tau","theta","Eta_k")
jags_save =c("beta0","psi","tau","theta")
jags.inits = function(){
  list("beta0"=rnorm(1),"psi"=rgamma(1,1,1),"Lambda"=rnorm(n_Y,3.0,0.5),"tau"=runif(1,1,10),"theta"=round(runif(1,0.5,max_range+0.5)),"Eta_k"=rnorm(n_k))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,n.iter=10000,
                model.file=neg.bin.model)

fit_mcmc <- as.mcmc(jags_fit)
#plot(fit_mcmc)
                #parameters.to.save=jags_save,
                #n.chains=3,
                #n.burnin=100,
                #n.thin=1,
                #DIC=FALSE)
plot(jags_fit$BUGSoutput$sims.matrix[,"theta"])
plot(jags_fit$BUGSoutput$sims.matrix[,"beta0"])
plot(jags_fit$BUGSoutput$sims.matrix[,"tau"])
plot(jags_fit$BUGSoutput$sims.matrix[,"psi"])


