pois.mcmc<-function(Data){ #simple Poisson mcmc metropolis-hastings algorithm w improper uniform prior
  Theta.mcmc=rep(0,1100)
  Theta.mcmc[1]=mean(Data)
  mh.inc=mean(Data/7)
  Accept=0
  for(iiter in 2:1100){
    prop=Theta.mcmc[iiter-1]+rnorm(1,0,mh.inc)
    if(prop>0 & runif(1)<exp(sum(dpois(Data,prop,log=1))-sum(dpois(Data,Theta.mcmc[iiter-1],log=1)))){
      Theta.mcmc[iiter]=prop
      Accept=Accept+1
    }
    else Theta.mcmc[iiter]=Theta.mcmc[iiter-1]
  } 
  cat(paste0("Acceptance rate = ",Accept/1100))
  return(Theta.mcmc[101:1100])
}

#2 nonspatial Poisson model used in spatial regression simulations
pois.no.spat <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Y_sim[i] ~ dpois(Lambda[i]) #posterior predictions
    log(Lambda[i]) <- XB[i]
    log_lik_pred[i] <- Y_sim[i]*log(Lambda[i])-Lambda[i]-logfact(Y_sim[i])
    log_lik[i] <- Y[i]*log(Lambda[i])-Lambda[i]-logfact(Y[i])
  }
  XB[1:n_Y] <- X[,] %*% Beta[1:n_B]
  # Priors
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  dev_pred <- -2*sum(log_lik_pred)
  dev_y <- -2*sum(log_lik)
}

pois.overd.no.spat <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Y_sim1[i] ~ dpois(Lambda[i]) #posterior predictions
    Y_sim2[i] ~ dpois(Lam_pred[i])
    Lambda[i] <- exp(Nu0[i])
    Lam_pred[i] <- exp(Nu_pred[i])
    Nu0[i] ~ dnorm(Mu[i],tau_iid)
    Nu_pred[i] ~ dnorm(Mu[i],tau_iid)
    log_lik_sim1[i] <- Y_sim2[i]*log(Lam_pred[i])-Lam_pred[i]-logfact(Y_sim2[i]) #poisson part
    log_lik1[i] <- Y[i]*log(Lambda[i])-Lambda[i]-logfact(Y[i])
    log_lik_sim2[i] <- -0.5*tau_iid*(Nu_pred[i]-Mu[i])^2  #gaussian part
    log_lik2[i] <- -0.5*tau_iid*(Nu0[i]-Mu[i])^2
  }
  Mu[1:n_Y] <- X[,] %*% Beta[1:n_B]
  # Priors
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  tau_iid~dgamma(1.0,0.01)
  dev_pred1 <- -2*sum(log_lik_sim1)
  dev_pred2 <- -2*sum(log_lik_sim2)
  dev_y1 <- -2*sum(log_lik1)  
  dev_y2 <- -2*sum(log_lik2) 
}

gauss.no.center <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dnorm(b0,tau_comb)
    Y_sim[i] ~ dnorm(b0,tau_comb) #posterior predictions
    #Z[i] <- b0 + X[i,]%*%Beta[1:n_B]
    log_lik[i]=(Y[i]-b0)^2
    log_lik_sim[i]=(Y_sim[i]-b0)^2
  }
  dev_y <- -n_Y*log(tau_comb)+sum(log_lik)*tau_comb
  dev_sim <- -n_Y*log(tau_comb)+sum(log_lik_sim)*tau_comb
  tau_comb<- tau_p * tau_iid / (tau_p + tau_iid)
  b0 ~ dnorm(0,0.01)
  #tau_iid~dgamma(1.0,0.01)
  tau_p~dgamma(1.0,0.01)
}

gauss.center <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dnorm(Z[i],tau_iid)
    Y_sim[i] ~ dnorm(Z[i],tau_iid) #posterior predictions
    Y_rep[i] ~ dnorm(b0,tau_comb) #posterior predictions
    Z[i] ~ dnorm(b0,tau_p)
    log_lik1[i]=(Y[i]-b0)^2
    log_lik2[i]=(Y[i]-Z[i])^2
    log_lik_sim1[i]=(Y_sim[i]-b0)^2
    log_lik_sim2[i]=(Y_sim[i]-Z[i])^2
    log_lik_rep[i]=(Y_rep[i]-b0)^2
  }
  dev_y1 <- n_Y*log(tau_comb)+sum(log_lik1)*tau_comb
  dev_y2 <- n_Y*log(tau_iid)+sum(log_lik2)*tau_iid
  dev_sim1 <- n_Y*log(tau_comb)+sum(log_lik_sim1)*tau_comb
  dev_sim2 <- n_Y*log(tau_iid)+sum(log_lik_sim2)*tau_iid
  dev_rep <- n_Y*log(tau_comb)+sum(log_lik_rep)*tau_comb
  tau_comb <- tau_iid * tau_p / (tau_iid + tau_p)
  # Priors
  #for(i in 1:n_B){
  #  Beta[i] ~ dnorm(0,tau_p)
  #}
  b0 ~ dnorm(0,0.01)
  #tau_iid~dgamma(1.0,0.01)
  tau_p~dgamma(1.0,0.01)
}


#cleaner version of poisson-normal mixture JAGS model for cross validation 
pois.overd.no.spat.cv <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Y_sim[i] ~ dpois(Lam_sim[i])
    Lambda[i] <- exp(Nu0[i])
    Lam_sim[i] <- exp(Nu_sim[i])
    Nu0[i] ~ dnorm(Mu[i],tau_iid)
    Nu_sim[i] ~ dnorm(Mu[i],tau_iid)
  }
  for (i in 1:n_k){
    Y_pred[i] ~ dpois(Lam_pred[i])
    Lam_pred[i] <- exp(Nu_pred[i])
    Nu_pred[i] ~ dnorm(Mu_pred[i],tau_iid)
  }
  Mu[1:n_Y] <- X[,] %*% Beta[1:n_B]
  Mu_pred[1:n_k] <- X_pred[,] %*% Beta[1:n_B]
  # Priors
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  tau_iid~dgamma(1.0,0.01)
}

### poisson model with spatially autocorrelated REs (blocky predictive process) AND IID Gaussian REs
pois.spat.model <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Y_sim1[i] ~ dpois(Lambda[i]) #posterior predictions
    Y_sim2[i] ~ dpois(Lam_pred[i])
    Lambda[i] <- exp(Nu0[i])
    Lam_pred[i] <- exp(Nu_pred[i])
    Nu0[i] ~ dnorm(Mu[i],tau_iid)
    Nu_pred[i] ~ dnorm(Mu[i],tau_iid)
    log_lik_sim1[i] <- Y_sim2[i]*log(Lam_pred[i])-Lam_pred[i]-logfact(Y_sim2[i]) #poisson part
    log_lik1[i] <- Y[i]*log(Lambda[i])-Lambda[i]-logfact(Y[i])
    log_lik_sim2[i] <- -0.5*tau_iid*(Nu_pred[i]-Mu[i])^2  #gaussian part
    log_lik2[i] <- -0.5*tau_iid*(Nu0[i]-Mu[i])^2
  }
  Mu[1:n_Y] <- X[,] %*% Beta[1:n_B] + Eta_y[1:n_Y]
  Eta_y[1:n_Y] <- X_RE[1:n_Y,1:n_k,theta]%*%Eta_k[1:n_k]
  Eta_k[1:n_k] ~ dmnorm(Zero_k[1:n_k],Sigma_inv[1:n_k,1:n_k])
  Sigma_inv[1:n_k,1:n_k] <- tau_sp*Q_k[1:n_k,1:n_k,theta]
  
  # Priors
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  tau_iid~dgamma(1.0,0.01)
  tau_sp~dgamma(1.0,0.01)
  theta~dcat(P)
  dev_pred1 <- -2*sum(log_lik_sim1)
  dev_pred2 <- -2*sum(log_lik_sim2)
  dev_y1 <- -2*sum(log_lik1)  
  dev_y2 <- -2*sum(log_lik2) 
}


# Inputs:
#   theta - Vector (or array) containing theta1, theta2, and theta3.
#             Theta1 = range.
#             Theta2 = smoothness.
#             Theta3 = nugget.
#   D     - Matrix (or array) of distances between sites. Typically D is an
#           (nxn) matrix with zeros along the main diagonal.
cor.matern <- function(theta,D) {  #from Hoeting et al. 2005 (Ecological Applications)
  if (is.vector(D)) names(D) <- NULL
  if (is.matrix(D)) dimnames(D) <- list(NULL,NULL)
  th1 <- theta[1]                             # Range.
  th2 <- theta[2]                             # Smoothness.
  th3 <- ifelse(length(theta)==3,theta[3],0)  # Nugget effect.
  u <- 2*D*sqrt(th2)/th1
  u <- ifelse(u>0,(1-th3)*u^th2/(2^(th2-1)*gamma(th2))*besselK(u,th2),1)
  return(u)
}

library(MASS)
library(fields)
library(mvtnorm)

n_Y=100
Knots=expand.grid(c(0:5),c(0:5))
n_k=nrow(Knots)
Coords=matrix(runif(n_Y*2,0,10),n_Y,2)
theta_true=1 
tau_true=1
#Dyy=rdist(Coords,Coords)
Dyk=rdist(Coords,Knots)
Dkk=rdist(Knots,Knots)
#Cov_kk=cor.matern(theta_true,Dkk)/tau_true
#Q_k = solve(cor.matern(theta_true,Dkk)/tau_true)
#X_RE=(cor.matern(theta_true,Dyk)/tau_true)%*% Q_k
Cov_kk=Exp.cov(Knots,theta=theta_true,distMat=Dkk)/tau_true
Cov_yk=Exp.cov(Coords,Knots,theta=theta_true,distMat=Dyk)/tau_true
X_RE = Cov_yk %*% solve(Cov_kk)
Eta_k_true=matrix(rmvnorm(1,rep(0,n_k),Cov_kk),ncol=1)
Eta_true=X_RE%*%Eta_k_true
#Cov_yy=cor.matern(theta_true,Dyy)/tau_true
#Eta_true=rmvnorm(1,rep(0,n_Y),Cov_yy)
Sim_lambda=matrix(exp(3+Eta_true)*rgamma(n_Y,1,rate=1),ncol=1)
#Sim.lambda=matrix(rgamma(100,1,1/exp(3)) ,ncol=1)
crap=rpois(n_Y,Sim_lambda)
hist(crap)

#knots is n_knots x 2 matrix giving easting, northing locations of knots
#coords in n_Y x 2 matrix giving easting, northing locations of data
negbin.reg.mcmc <- function(Y,X=NULL,Control,Prior=NULL,Coords=NULL,Knots=NULL){
  DEBUG=FALSE
  n_Y=length(Y)
  SPAT=FALSE
  Eta_y=rep(0,n_Y)
  if(is.null(X))X=matrix(1,nrow=n_Y,ncol=1) #intercept only model if X is omitted
  n_beta=ncol(X)
  MCMC=list(Beta=matrix(0,n_beta,Control$iter),psi=rep(0,Control$iter))
  Accept=list(Beta=rep(0,n_beta),psi=0)
  if(is.null(Control$inits$Beta)){
    Control$inits$Beta=rep(0,n_beta)
    Control$inits$Beta[1]=c(log(mean(Y))+rnorm(1))
    if(n_beta>1)Control$inits$Beta[2:n_beta]=rnorm(n_beta-2+1)
  }
  if(is.null(Control$inits$psi))Control$inits$psi=runif(1,1,10)
  if(is.null(Control$inits$Lambda))Control$inits$Lambda=matrix(exp(log(mean(Y))+rnorm(n_Y,0,2)),ncol=1)
  MCMC$Beta[,1]=Control$inits$Beta
  MCMC$psi[1]=Control$inits$psi
  Lambda=Control$inits$Lambda
  if(is.null(Control$MH_beta))Control$MH_beta=0.1
  if(is.null(Control$MH_psi))Control$MH_psi=0.1
  if(is.null(Control$burnin))Control$burnin=0
  if(is.null(Prior$a))Prior$a=1
  if(is.null(Prior$b))Prior$b=0.01
  if(DEBUG)MCMC$Beta[,1]=3
  if(DEBUG)MCMC$psi[1]=1
  Beta=MCMC$Beta[,1]

  if(is.null(Knots)==FALSE){
    SPAT=TRUE
    if(is.null(Coords))stop("coords must be provided whenever knots are provided")
    n_k=nrow(Knots)
    MCMC$theta=MCMC$tau=rep(0,Control$iter)
    max_range=max(Dkk)/3
    if(is.null(Control$inits$theta)==TRUE)Control$inits$theta=c(runif(1,0,max_range))
    MCMC$theta[1]=Control$inits$theta
    MCMC$tau[1]=runif(1,0.5,20)  
    MCMC$tau[1]=1
    Dyk=rdist(Coords,Knots)
    Dkk=rdist(Knots,Knots)
    Cor_k=Q_k=X_RE=vector("list",100)
    incr=max_range/100
    range_index=round(runif(1,0.5,100.5))
    for(i in 1:100){  #compute exponential correlation function, inverses, etc., over a set of discretized range values
      Cor_k[[i]]=Exp.cov(Knots,theta=incr*i,distMat=Dkk)
      Q_k[[i]]=solve(Cor_k[[i]])
      X_RE[[i]]=Exp.cov(Coords,Knots,theta=incr*i,distMat=Dyk)%*%Q_k[[i]]
    }
    #Q_k = solve(Sigma_k)
    #X_RE=(Exp.cov(Coords,Knots,theta=MCMC$Theta[1,1],distMat=Dyk)/MCMC$Theta[2,1])%*%Q_k
    #Sigma_k=cor.matern(MCMC$Theta[c(1:2),1],Dkk)/MCMC$Theta[3,1]
    #Q_k = solve(cor.matern(MCMC$Theta[c(1:2),1],Dkk)/MCMC$Theta[3,1])
    #X_RE=(cor.matern(MCMC$Theta[c(1:2),1],Dyk)/MCMC$Theta[3,1])%*% Q_k
    Eta_k= matrix(rmvnorm(1,rep(0,n_k),Cor_k[[range_index]]/MCMC$tau[1]),ncol=1)
    Eta_y= X_RE[[range_index]] %*% Eta_k
    Accept$theta=Accept$tau=0
    Accept$Eta=rep(0,n_k)
    if(is.null(Control$MH_Eta))Control$MH_Eta=rep(1.5,n_Y)
    if(is.null(Control$MH_theta))Control$MH_theta=20
    if(is.null(Control$MH_tau))Control$MH_tau=.3
  }
  Time=rep(0,5)
  names(Time)=c("Beta","psi and Lambda","Eta","theta","tau")
  for(iiter in 2:Control$iter){
    if(iiter%%1000==0)cat(paste("Iteration ",iiter," out of ",Control$iter,"\n"))
    MCMC$Beta[,iiter]=MCMC$Beta[,iiter-1]
    MCMC$psi[iiter]=MCMC$psi[iiter-1]
    #update beta
    cur.time=proc.time()
    if(!DEBUG){
    for(ibeta in 1:n_beta){
      beta_prop=MCMC$Beta[ibeta,iiter]+rnorm(1,0,Control$MH_beta)
      Beta_tmp=Beta
      Beta_tmp[ibeta]=beta_prop
      log_FC_new=-MCMC$psi[iiter]*(sum(X[,ibeta])*beta_prop+crossprod(Lambda,exp(-X%*%Beta_tmp-Eta_y))) #full conditional
      log_FC_old=-MCMC$psi[iiter]*(sum(X[,ibeta])*Beta[ibeta]+crossprod(Lambda,exp(-X%*%Beta-Eta_y)))
      if(runif(1)<exp(log_FC_new-log_FC_old)){
        Beta[ibeta]=beta_prop
        MCMC$Beta[ibeta,iiter]=beta_prop
        Accept$Beta[ibeta]=Accept$Beta[ibeta]+1
      }
    }
    }
    XB=X%*%Beta
    Mu=exp(XB+Eta_y)
    Time[1]=Time[1]+(proc.time()-cur.time)[3]
    cur.time=proc.time()
    
    #log.FC.old=-n.Y*lgamma(cur.psi)+cur.psi*(n.Y*log(cur.psi)-sum(XB+log(Lambda)+Lambda*exp(-XB)))+(Prior$a-1)*log(cur.psi)-Prior$b*cur.psi #from Lord and Park, Appendix D; doesn't appear to work
    #log.FC.old=sum(dgamma(Lambda,cur.psi,rate=cur.psi/exp(XB),log=TRUE))+(Prior$a-1)*log(cur.psi)-Prior$b*cur.psi  #works but likely slower than direct

    #update psi
    if(!DEBUG){
    log_FC_old=-n_Y*lgamma(MCMC$psi[iiter])+MCMC$psi[iiter]*sum(log(MCMC$psi[iiter])-XB-Eta_y-Lambda/Mu)+sum((MCMC$psi[iiter]-1)*log(Lambda))+(Prior$a-1)*log(MCMC$psi[iiter])-Prior$b*MCMC$psi[iiter]
    cur_psi=MCMC$psi[iiter]+rnorm(1,0,Control$MH_psi)
    if(cur_psi>0){
      log_FC_new=-n_Y*lgamma(cur_psi)+cur_psi*sum(log(cur_psi)-XB-Eta_y-Lambda/Mu)+sum((cur_psi-1)*log(Lambda))+(Prior$a-1)*log(cur_psi)-Prior$b*cur_psi
      if(runif(1)<exp(log_FC_new-log_FC_old)){
        Accept$psi=Accept$psi+1
        MCMC$psi[iiter]=cur_psi
      }
    }
    }
    #update lambda
    #Lambda=Sim.lambda
    Lambda=matrix(rgamma(n_Y,Y+MCMC$psi[iiter],rate=1+MCMC$psi[iiter]*exp(-XB-Eta_y)),ncol=1)
    Time[2]=Time[2]+(proc.time()-cur.time)[3]
    cur.time=proc.time()
    
    
    #Eta_k = Eta_k_true
    if(SPAT){#update spatial random effects and autocorrelation parameters
      MCMC$theta[iiter]=MCMC$theta[iiter-1]
      MCMC$tau[iiter]=MCMC$tau[iiter-1]
      if(!DEBUG){
      XE=X_RE[[range_index]]%*%Eta_k
      #log_FC_old = sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu,log=TRUE))- 0.5*(crossprod(Eta_k,Q_k[[range_index]]) %*% Eta_k)*MCMC$tau[iiter]
      log_FC_old = sum(-MCMC$psi[iiter]*(Lambda/Mu+XE))- 0.5*(crossprod(Eta_k,Q_k[[range_index]]) %*% Eta_k)*MCMC$tau[iiter]
      Mu_tmp=Mu
      Eta_tmp=Eta_k
      Rnorm=rnorm(n_Y,0,Control$MH_Eta)
      ## update eta_k
      if(iiter%%10 == 2){
      for(ieta in 1:n_k){
        Eta_tmp [ieta] = Eta_tmp[ieta] + Rnorm[ieta]
        XE_prop = X_RE[[range_index]]%*%Eta_tmp
        Mu_tmp= exp(XB + XE_prop)
        #log_FC_new = sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu_tmp,log=TRUE))- 0.5*(crossprod(Eta_tmp,Q_k[[range_index]]) %*% Eta_tmp)*MCMC$tau[iiter]
        log_FC_new = sum(-MCMC$psi[iiter]*(Lambda/Mu_tmp+XE_prop))- 0.5*(crossprod(Eta_tmp,Q_k[[range_index]]) %*% Eta_tmp)*MCMC$tau[iiter]
        if(runif(1)<exp(log_FC_new-log_FC_old)){
          Mu=Mu_tmp
          Eta_k[ieta]=Eta_tmp[ieta]
          log_FC_old=log_FC_new
          Accept$Eta[ieta]=Accept$Eta[ieta]+1
        }
        else{
          Mu_tmp=Mu
          Eta_tmp[ieta]=Eta_k[ieta]
        }
      }
      }
      }
      
      if(DEBUG)Eta_k=Eta_k_true
      Eta_y = X_RE[[range_index]]%*%Eta_k
      Time[3]=Time[3]+(proc.time()-cur.time)[3]
      cur.time=proc.time()
      
      
      ### update theta (exponential covariance range parameter) 
      range_index_prop=range_index+round(rnorm(1,0,Control$MH_theta))
      if(range_index_prop>0 & range_index_prop<101){
        theta_prop=incr*range_index_prop
        log_FC_old=-0.5*( log(det(Cor_k[[range_index]])) + crossprod(Eta_k,Q_k[[range_index]]) %*% Eta_k * MCMC$tau[iiter])+sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu,log=TRUE))
        Mu_tmp=exp(XB + X_RE[[range_index_prop]]%*%Eta_k)
        log_FC_new=-0.5*( log(det(Cor_k[[range_index_prop]])) + crossprod(Eta_k,Q_k[[range_index_prop]]) %*% Eta_k * MCMC$tau[iiter])+sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu_tmp,log=TRUE))
        if(runif(1)<exp(log_FC_new-log_FC_old)){
          range_index=range_index_prop
          MCMC$theta[iiter]=theta_prop
          Mu=Mu_tmp
          Accept$theta=Accept$theta+1
          Eta_y = X_RE[[range_index]]%*%Eta_k
        }        
      }
      Time[4]=Time[4]+(proc.time()-cur.time)[3]
      cur.time=proc.time()
      
      
      #update precision (tau) for spatial random effects
      #note don't have to include component for lambda because tau cancels out of X_RE
      #comp1=det(Cor_k[[range_index]])
      comp2=crossprod(Eta_k,Q_k[[range_index]]) %*% Eta_k
      log_FC_old=-0.5*( -n_k*log(MCMC$tau[iiter]) + comp2*MCMC$tau[iiter])+dgamma(MCMC$tau[iiter],Prior$a,Prior$b,log=TRUE)
      tau_prop=MCMC$tau[iiter]+rnorm(1,0,Control$MH_tau)
      if(tau_prop>0){
        log_FC_new=-0.5*(-n_k*log(tau_prop) + comp2*tau_prop)+dgamma(tau_prop,Prior$a,Prior$b,log=TRUE)
        if(runif(1)<exp(log_FC_new-log_FC_old)){
          MCMC$tau[iiter]=tau_prop
          Accept$tau=Accept$tau+1
        }
      }
      Time[5]=Time[5]+(proc.time()-cur.time)[3]
      cur.time=proc.time()
      
      
      #log_FC_old=dgamma(MCMC$Theta[2,iiter],1,0.01,log=TRUE)-0.5*( log(det(Sigma_k)/) + crossprod(Eta_k,Q_k) %*% Eta_k)+sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu,log=TRUE))
      #Sigma_tmp=Exp.cov(Knots,theta=Theta_prop[1],distMat=Dkk)/Theta_prop[2]
      #Q_tmp = solve(Sigma_tmp)
      #X_RE_tmp=(Exp.cov(Coords,Knots,theta=Theta_prop[1],distMat=Dyk)/Theta_prop[2])%*%Q_tmp
      #Mu_tmp= exp(XB + X_RE_tmp%*%Eta_k)
      #log_FC_new=dgamma(Theta_prop[2],1,0.01,log=TRUE)-0.5*(log(det(Sigma_tmp)) + crossprod(Eta_k,Q_tmp) %*% Eta_k)+sum(dgamma(Lambda,MCMC$psi[iiter],MCMC$psi[iiter]/Mu_tmp,log=TRUE))
      #if(runif(1)<exp(log_FC_new-log_FC_old)){
      #  Q_k=Q_tmp
      #  Sigma_k=Sigma_tmp
      #  X_RE = X_RE_tmp
      #  Mu=Mu_tmp
      #  MCMC$Theta[,iiter]=Theta_prop
      #  Accept$theta=Accept$theta+1
      #  Eta_y = X_RE%*%Eta_k
      #}
    }
  }
  return(list(MCMC=MCMC,Accept=Accept,Lambda=Lambda,Eta_k=Eta_k,Time=Time))
}

#set.seed(11111)
#poop=negbin.reg.mcmc(Y=crap,X=matrix(1,nrow=length(crap),ncol=1),Control=list(iter=20000),Coords=Coords,Knots=Knots)
#acf(poop$MCMC$Beta[1,])
#autocorr_time=1+2*sum(acf(poop$MCMC$Beta[1,],lag.max=800)$acf) #take lag.max= point where autocorr = 0 
#eff_n = length(poop$MCMC$Beta[1,])/autocorr_time

norm.mcmc<-function(Data){
  Theta.mcmc=matrix(0,2,1100)
  Theta.mcmc[,1]=c(mean(Data),1/var(Data))
  mu.prior.sd=100
  mu.prior.var=mu.prior.sd^2
  mu.prior.prec=1/mu.prior.var
  n=length(Data)
  y.bar=mean(Data)
  alpha0=1
  beta0=0.01
  for(iiter in 2:1100){
    #use straight Gibbs
    #mu update
    mu=mu.prior.sd^2*n*y.bar/(Theta.mcmc[2,iiter-1]+n*mu.prior.var)
    mu.var=1/(mu.prior.prec+n*Theta.mcmc[2,iiter-1])
    Theta.mcmc[1,iiter]=rnorm(1,mu,sqrt(mu.var))
    #tau (precision) update
    Theta.mcmc[2,iiter]=rgamma(1,0.5*n+alpha0,0.5*sum((Data-Theta.mcmc[1,iiter])^2))
  }
  Theta.mcmc
}
