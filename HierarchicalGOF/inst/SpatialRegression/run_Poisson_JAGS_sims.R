##### Run spatial poisson regression model simulations using JAGS


####   set simulation specifications ###########
n_sims = 1 #number of simulation replicates
n_Y=200  #number of sample locations
Knots=expand.grid(c(0:8),c(0:8))/2  #knots for predictive process spatial model
n_k=nrow(Knots)
theta_true=2 #controls exponenial correlation decay
tau_eta_true=1  #precision of spatial autocorrelation process
tau_iid_true=5  #iid precision for extra-Poisson error
incr=(max(Knots)-min(Knots))/99
Beta_true=c(2,0.75)  #intercept, slope for relationship of log abundance to simulated covariate
#MCMC options
n_iter = 2000
n_burnin = 500  
n_thin = 1
n_chains=3  
n_mcmc=n_chains*(n_iter-n_burnin)/n_thin #total number of recorded MCMC iterations

##### Initialize data structures for storing simulation results ######
Colnames=c("pois0_postp_chisq","pois0_postp_ft","pois0_postp_dev","pois0_postp_tail","pois0_postp_0","pois0_moran","pois0_pivot","pois0_sampled",
           "poisIID_postp_chisq","poisIID_postp_ft","poisIID_postp_dev1","poisIID_postp_dev2","poisIID_postp_tail","poisIID_postp_0","poisIID_moran","poisIID_pivot1","poisIID_pivot2","poisIID_sampled1","poisIID_sampled2","poisIID_cv","poisIID_cv_outliers",
           "poisSpat_postp_chisq","poisSpat_postp_ft","poisSpat_postp_ft_mixed","poisSpat_postp_dev1","poisSpat_postp_dev2","poisSpat_postp_tail","poisSpat_postp_0","poisSpat_moran","poisSpat_pivot1","poisSpat_pivot2","poisSpat_sampled1","poisSpat_sampled2")
Results = data.frame(matrix(0,n_sims,length(Colnames)))
colnames(Results) = Colnames

#### Initialize some structures needed for spatially autocorrelated Poisson estimation model
max_range=4
n_incr=10
incr=max_range/n_incr
Cor_k=Q_k=array(0,dim=c(n_k,n_k,n_incr))
Q_k_partial = array(0,dim=c(n_k,n_k-1,n_k-1,n_incr)) #for partial posterior prediction
Cor_partial = array(0,dim=c(n_k,n_k-1,n_incr)) 
X_RE=array(0,dim=c(n_Y,n_k,n_incr))
Tmp=c(1:n_k)
P=rep(incr,n_incr)*c(1:n_incr)

var.fcn.pois = mean.fcn.pois = function(x){
  return(x)
}
mean.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(Theta[1:n_par])
}
var.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(1/Theta[n_par])
}

mean.fcn.pois<-var.fcn.pois<-function(Theta){
  return(Theta)
}

my.pnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(pnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}
my.dnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(dnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}

my.ppois<-function(Y,Theta)ppois(Y,Theta)
my.dpois<-function(Y,Theta)dpois(Y,Theta)


#####  BEGIN SIMULATIONS #####
cur_time = proc.time()
for(isim in 1:n_sims){
  set.seed(10000+isim)  #so random number seed known for each simulation
  cat(paste('Simulation ',isim,' out of ',n_sims,'\n'))
  if(isim ==2){
    elapsed=as.numeric(proc.time()-cur_time)[3]
    cat(paste("Time for 1 simulation: ",elapsed/3600," hours \n"))
    cat(paste("Estimated time remaining: ",elapsed/3600*(n_sims-1)," hours \n"))
  }
  #determine sample locations
  Coords=matrix(runif(n_Y*2,0.5,3.5),n_Y,2)  
  
  #distances between data points and knots, etc.
  Dyy=rdist(Coords,Coords)
  Dyk=rdist(Coords,Knots)
  Dkk=rdist(Knots,Knots)
  
  #covariance matrices
  Cov_kk=Exp.cov(Knots,theta=theta_true,distMat=Dkk)/tau_eta_true
  Cov_yk=Exp.cov(Coords,Knots,theta=theta_true,distMat=Dyk)/tau_eta_true
  #Cov_pk=Exp.cov(Coords_plot,Knots,theta=theta_true,distMat=Dpk)/tau_eta_true
  
  #spatial random effects "design matrix"
  X_RE_true = Cov_yk %*% solve(Cov_kk)
  
  #generate spatially autocorrelated covariate
  Matern_points = rMatClust(kappa=12,r=.25,mu=100)  #generate some points from a matern process to construct a hypothetical covariate
  Matern_xy=cbind(Matern_points$x,Matern_points$y)*max(Knots[,1])
  Matern_kde=kde2d(Matern_xy[,1],Matern_xy[,2],100,h=0.5,lims=c(0,max(Knots[,1]),0,max(Knots[,2])))
  Wts=1/(rdist(Coords,expand.grid(x=Matern_kde$x,y=Matern_kde$y)))^2
  Wts=Wts*1/rowSums(Wts)
  Pred_cov =  Wts %*% matrix(Matern_kde$z,ncol=1) #hypothetical habitat covariate
  denom=mean(Pred_cov)
  Pred_cov=Pred_cov/denom  #stan Wts %*% as.vector(Matern_kde$z)dardize to have a mean on 1.0

  #design matrix for fixed effects
  X=matrix(1,nrow=n_Y,ncol=2)  
  X[,2]=Pred_cov  # linear term only

  #generate spatial random effects
  Eta_k_true=matrix(rmvnorm(1,rep(0,n_k),Cov_kk),ncol=1)
  Eta_true=X_RE_true%*%Eta_k_true
  
  #simualte observed counts
  Sim_lambda=exp(rnorm(n_Y,X%*%Beta_true+Eta_true,sqrt(1/tau_iid_true)))
  Y=rpois(n_Y,Sim_lambda)
  Y_mat=matrix(Y,n_Y,n_mcmc)  #needed for pivot function

  #conduct MCMC analysis of counts using a Poisson regression model without random effects
  n_B=ncol(X)
  jags_data = list("n_Y","n_B","Y","X")
  jags_params = c("Beta","Y_sim")
  jags_save =c("Beta","Y_sim","Lambda","dev_pred","deviance","dev_y")
  jags.inits = function(){
    list("Beta"=rnorm(n_B),"Y_sim"=rgamma(n_Y,10,1))
  }
  jags_fit = jags(data=jags_data,
                  inits=jags.inits,
                  jags_params,
                  n.iter=n_iter,
                  model.file=pois.no.spat,
                  DIC=FALSE,
                  parameters.to.save=jags_save,
                  n.chains=n_chains,
                  n.burnin=n_burnin,
                  n.thin=n_thin,
                  working.directory=getwd())
  
  # calculate goodness-of-fit metrics on basic Poisson regressino model
  MoranI=rep(0,n_mcmc)
  #select appropriate columns form jags MCMC output matrix
  Cols_Y_sim=grep("Y_sim",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Lambda=grep("Lambda",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Beta=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
  Dists_inv=1/Dyy   #(Dyy<.38)*1
  diag(Dists_inv)=0
  p_omni=p_tail=p_0=p_ft=0
  y_quant95 = quantile(Y,0.95)
  y_0 = sum(Y==0)
  for(i in 1:n_mcmc){
    MoranI[i]=Moran.I(Y-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda],Dists_inv)$observed
    chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])
    chisq_y=sum((Y-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])  
    ft=sum((sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim])-sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda]))^2)
    ft_y=sum((sqrt(Y)-sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda]))^2)
    p_omni = p_omni+(chisq_y<chisq)
    p_ft = p_ft + (ft_y<ft)
    p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim],0.95))
    p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]==0)) + 0.5
  }
  Results[isim,"pois0_postp_chisq"]=p_omni/n_mcmc
  Results[isim,"pois0_postp_tail"] = p_tail/n_mcmc
  Results[isim,"pois0_postp_0"] = p_0/n_mcmc
  Results[isim,"pois0_postp_ft"] = p_ft/n_mcmc
  Results[isim,"pois0_postp_dev"] = sum(jags_fit$BUGSoutput$sims.matrix[,"deviance"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred"])/n_mcmc
  Results[isim,"pois0_moran"]=sum(MoranI<0)/n_mcmc
  pivot_out = chisq.cdf.test(Y=Y_mat,Theta=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Lambda]),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)
  Results[isim,"pois0_sampled"]=sample(pivot_out$pval_omni,1)
  Results[isim,"pois0_pivot"]=median(pivot_out$pval_omni)
    
  #Run poisson model with IID Gaussian random effects
  jags_data = list("n_Y","n_B","Y","X")
  jags_params = c("Beta","tau_iid","Nu0","Nu_pred","Y_sim1","Y_sim2")
  jags_save =c("Beta","Y_sim1","Y_sim2","Lambda","Lam_pred","Nu0","Mu","deviance","dev_pred1","dev_y1","dev_pred2","dev_y2","tau_iid")
  jags.inits = function(){
    #list("Beta"=rnorm(n_B),"psi"=rgamma(1,1,1),"Nu"=rgamma(n_Y,1,1),"Nu_pred"=rgamma(n_Y,1,1))
    list("Beta"=rnorm(n_B),"tau_iid"=rgamma(1,1,.1),"Nu0"=log(Y+0.1)+rnorm(n_Y),"Nu_pred"=log(Y+0.1)+rnorm(n_Y),"Y_sim1"=exp(log(Y+0.1)+rnorm(n_Y)),"Y_sim2"=exp(log(Y+0.1)+rnorm(n_Y)))
  }
  jags_fit = jags(data=jags_data,
                  inits=jags.inits,
                  jags_params,
                  n.iter=n_iter,
                  model.file=pois.overd.no.spat,
                  DIC=FALSE,
                  parameters.to.save=jags_save,
                  n.chains=n_chains,
                  n.burnin=n_burnin,
                  n.thin=n_thin,
                  working.directory=getwd())
  
  ##### calculate Moran's I for residuals, Chi-square omnibus test, Targeted skewness check
  MoranI=rep(0,n_mcmc)
  #select appropriate columns form jags MCMC output matrix
  Cols_Y_sim1=grep("Y_sim1",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Y_sim2=grep("Y_sim2",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Lambda=grep("Lambda",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Lam_Pred=grep("Lam_pred",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Mu=grep("Mu",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Nu=grep("Nu0",colnames(jags_fit$BUGSoutput$sims.matrix)) 
  Cols_Beta=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
  Dists_inv=1/Dyy
  diag(Dists_inv)=0
  p_omni=p_tail=p_0=p_ft=0
  for(i in 1:n_mcmc){
    Mu_tmp = exp(jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])*exp(0.5/jags_fit$BUGSoutput$sims.matrix[i,"tau_iid"]) #includes lognormal bias correction
    MoranI[i]=Moran.I(Y-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lam_Pred],Dists_inv)$observed
    chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]-Mu_tmp)^2/Mu_tmp)
    chisq_y=sum((Y-Mu_tmp)^2/Mu_tmp)  
    ft=sum((sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2])-sqrt(Mu_tmp))^2)
    ft_y=sum((sqrt(Y)-sqrt(Mu_tmp))^2)
    p_ft = p_ft + (ft_y<ft)
    p_omni = p_omni+(chisq_y<chisq)
    p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2],0.95))
    p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]==0)) 
  }
  Results[isim,"poisIID_postp_chisq"]=p_omni/n_mcmc
  Results[isim,"poisIID_postp_tail"] = p_tail/n_mcmc
  Results[isim,"poisIID_postp_0"] = p_0/n_mcmc
  Results[isim,"poisIID_moran"]=sum(MoranI<0)/n_mcmc
  Results[isim,"poisIID_postp_dev1"] = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y1"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred1"])/n_mcmc
  Results[isim,"poisIID_postp_dev2"] = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y2"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred2"])/n_mcmc
  Results[isim,"poisIID_postp_ft"] = p_ft/n_mcmc
  #test poisson part of poisson-normal mixture
  pivot_out=chisq.cdf.test(Y=Y_mat,Theta=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Lambda]),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)
  Results[isim,"poisIID_sampled1"]=sample(pivot_out$pval_omni,1)
  Results[isim,"poisIID_pivot1"]=median(pivot_out$pval_omni)
    #test normal part of poisson-normal mixture using pivotal discrepancy
  pivot_out=chisq.cdf.test(Y=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu]),
                           Theta=t(cbind(jags_fit$BUGSoutput$sims.matrix[,Cols_Mu],
                                         jags_fit$BUGSoutput$sims.matrix[,"tau_iid"])),
                           mean.fcn=mean.fcn.normal,var.fcn=var.fcn.normal,
                           cdf.fcn=my.pnorm,pmf.fcn=my.dnorm,K=5,L=5,DISCRETE=FALSE)
  Results[isim,"poisIID_sampled2"]=sample(pivot_out$pval_omni,1)
  Results[isim,"poisIID_pivot2"]=median(pivot_out$pval_omni)
  #conduct k-fold cross validation (40 folds of n=5)
  n_cv_k = 5
  n_cv_reps = n_Y/n_cv_k
  n_cv_Y = n_Y - n_cv_k
  Samples = matrix(0,n_Y/n_cv_k,n_cv_k)
  Cur_y = c(1:n_Y)
  for(i in 1:n_cv_reps){  #generate test groupings for k-fold cross validation
    Samples[i,]=sample(Cur_y,n_cv_k)
    Cur_y = Cur_y[-which(Cur_y %in% Samples[i,])]
  }
  D_cv = rep(0,n_Y)
  for(i in 1:n_cv_reps){
    which_Y = c(1:n_Y)[-Samples[i,]]
    Y_cv = Y[which_Y]
    X_cv = X[which_Y,]
    X_pred = X[Samples[i,],]
    jags_data = list("n_Y"=n_cv_Y,"n_B"=n_B,"n_k"=n_cv_k,"Y"=Y_cv,"X"=X_cv,"X_pred"=X_pred)
    jags_params = c("Beta","tau_iid","Nu0","Nu_sim","Nu_pred","Y_sim","Y_pred")
    jags_save =c("Lam_pred")
    jags.inits = function(){
      list("Beta"=rnorm(n_B),"tau_iid"=rgamma(1,1,.1),"Nu0"=log(Y_cv+0.1)+rnorm(n_cv_Y),"Nu_sim"=log(Y_cv+0.1)+rnorm(n_cv_Y),"Nu_pred"=log(mean(Y_cv)+0.1)+rnorm(n_cv_k),"Y_sim"=exp(log(Y_cv+0.1)+rnorm(n_cv_Y)),"Y_pred"=exp(log(mean(Y_cv))+rnorm(n_cv_k)))
    }
    jags_fit = jags(data=jags_data,
                    inits=jags.inits,
                    jags_params,
                    n.iter=n_iter,
                    model.file=pois.overd.no.spat.cv,
                    DIC=FALSE,
                    parameters.to.save=jags_save,
                    n.chains=n_chains,
                    n.burnin=n_burnin,
                    n.thin=n_thin,
                    working.directory=getwd())
    for(k in 1:n_cv_k){
      D_cv[(i-1)*n_cv_k+k] = sum((jags_fit$BUGSoutput$sims.matrix[,k]<Y[Samples[i,k]])+0.5*(jags_fit$BUGSoutput$sims.matrix[,k]==Y[Samples[i,k]]))
    }
  }
  D_cv = D_cv/n_mcmc
  Chisq_bounds = c(0:10)/10
  Chisq_bounds[1]=-1
  Chisq_bounds[11]=1.1
  D_cat=tabulate(as.numeric(cut(D_cv,Chisq_bounds)),nbins=10)
  T_cv = sum((D_cat - n_Y/10)^2)*10/n_Y #chi-square value w/ 10 bins 
  Results[isim,"poisIID_cv"]=1-pchisq(T_cv,df=9)  
  Results[isim,"poisIID_cv_outliers"]=sum(D_cv<0.001)
  
  #run Poisson model w/ IID Gaussian REs AND spatially autocorrelated REs (blocky predictive process formulation)  
  Zero_k = rep(0,n_k)
  #discretize covariance matrices and related quantities so all matrix inverse
  #calculations can be done up front and provided as input to JAGS
  for(i in 1:n_incr){
    Cor_k[,1:n_k,i]=Exp.cov(Knots,theta=incr*i,distMat=Dkk)
    Q_k[,1:n_k,i]=solve(Cor_k[,1:n_k,i])
    X_RE[,1:n_k,i]=Exp.cov(Coords,Knots,theta=incr*i,distMat=Dyk)%*%Q_k[,1:n_k,i]
    for(k in 1:n_k){
      Q_k_partial[k,,,i]=solve(Cor_k[Tmp[-k],Tmp[-k],i])
      Cor_partial[k,,i] = Cor_k[k,Tmp[-k],i]
    }
  }
  jags_data = list("n_Y","n_B","Y","X","Zero_k","Q_k","X_RE","n_k","P")
  jags_params = c("Beta","Nu0","Nu_pred","tau_sp","tau_iid","theta","Eta_k","Y_sim1","Y_sim2")
  jags_save =c("Beta","tau_sp","tau_iid","theta","Y_sim1","Y_sim2","Mu","Nu0","Nu_pred","Eta_k","dev_pred1","dev_y1","dev_pred2","dev_y2")
  jags.inits = function(){
    list("Beta"=rnorm(n_B),"Nu0"=rnorm(n_Y,3.0,0.5),"Nu_pred"=rnorm(n_Y,3.0,0.5),"tau_iid"=runif(1,1,10),"tau_sp"=runif(1,1,10),"theta"=round(runif(1,0.5,max_range+0.5)),"Eta_k"=rnorm(n_k))
  }
  jags_fit = jags(data=jags_data,
                  inits=jags.inits,
                  jags_params,
                  n.iter=n_iter,
                  model.file=pois.spat.model,
                  DIC=FALSE,
                  parameters.to.save=jags_save,
                  n.chains=n_chains,
                  n.burnin=n_burnin,
                  n.thin=n_thin,
                  working.directory=getwd())
  
  ##### calculate GOF tests for spatially autocorrelated model
  MoranI=rep(0,n_mcmc)
  #select appropriate columns form jags MCMC output matrix
  Cols_Y_sim1=grep("Y_sim1",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Y_sim2=grep("Y_sim2",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Mu=grep("Mu",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Nu=grep("Nu0",colnames(jags_fit$BUGSoutput$sims.matrix)) 
  Cols_Nu_pred = grep("Nu_pred",colnames(jags_fit$BUGSoutput$sims.matrix)) 
  Cols_Beta=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
  Cols_Eta_k = grep("Eta_k",colnames(jags_fit$BUGSoutput$sims.matrix))
  Lambda=exp(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu])
  Lam_pred = exp(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu_pred])
  Dists_inv=1/Dyy
  diag(Dists_inv)=0
  p_omni=p_tail=p_0=p_ft=0
  for(i in 1:n_mcmc){
    Mu_tmp = exp(jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])*exp(0.5/jags_fit$BUGSoutput$sims.matrix[i,"tau_iid"]) #includes lognormal bias correction
    MoranI[i]=Moran.I(Y-Lam_pred[i,],Dists_inv)$observed
    chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]-Mu_tmp)^2/Mu_tmp)
    chisq_y=sum((Y-Mu_tmp)^2/Mu_tmp)  
    ft=sum((sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2])-sqrt(Mu_tmp))^2)
    ft_y=sum((sqrt(Y)-sqrt(Mu_tmp))^2)
    p_ft = p_ft + (ft_y<ft)
    p_omni = p_omni+(chisq_y<chisq)
    p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2],0.95))
    p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]==0)) 
  }
  Results[isim,"poisSpat_postp_chisq"]=p_omni/n_mcmc
  Results[isim,"poisSpat_postp_tail"] = p_tail/n_mcmc
  Results[isim,"poisSpat_postp_0"] = p_0/n_mcmc
  Results[isim,"poisSpat_moran"]=sum(MoranI<0)/n_mcmc
  Results[isim,"poisSpat_postp_dev1"] = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y1"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred1"])/n_mcmc
  Results[isim,"poisSpat_postp_dev2"] = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y2"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred2"])/n_mcmc
  Results[isim,"poisSpat_postp_ft"] = p_ft/n_mcmc
  #test poisson part of poisson-normal mixture
  pivot_out=chisq.cdf.test(Y=Y_mat,Theta=t(Lambda),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)
  Results[isim,"poisSpat_sampled1"]=sample(pivot_out$pval_omni,1)
  Results[isim,"poisSpat_pivot1"]=median(pivot_out$pval_omni)
  #test Gaussian part of mixture
  pivot_out=chisq.cdf.test(Y=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu]),
                           Theta=t(cbind(jags_fit$BUGSoutput$sims.matrix[,Cols_Mu],
                                         jags_fit$BUGSoutput$sims.matrix[,"tau_iid"])),
                           mean.fcn=mean.fcn.normal,var.fcn=var.fcn.normal,
                           cdf.fcn=my.pnorm,pmf.fcn=my.dnorm,K=5,L=5,DISCRETE=FALSE)
  Results[isim,"poisSpat_sampled2"]=sample(pivot_out$pval_omni,1)
  Results[isim,"poisSpat_pivot2"]=median(pivot_out$pval_omni)
  ### post hoc mixed predictive p-value by resimulating spatial random effects
  p_omni_mixed = 0
  Eta_tmp=rep(0,n_k)
  for(i in 1:n_mcmc){
    theta=jags_fit$BUGSoutput$sims.matrix[i,"theta"]
    Beta=jags_fit$BUGSoutput$sims.matrix[i,Cols_Beta]
    for(ik in 1:n_k){
      Tmp_cor = matrix(Cor_partial[ik,,theta],nrow=1)
      Eta_tmp[ik]=rnorm(1,Tmp_cor%*%Q_k_partial[ik,,,theta]%*%jags_fit$BUGSoutput$sims.matrix[i,Cols_Eta_k][-ik],
                        sqrt((1-Tmp_cor%*%Q_k_partial[ik,,,theta]%*%t(Tmp_cor))/jags_fit$BUGSoutput$sims.matrix[i,"tau_sp"]))
    }
    Mu_tmp=X%*%Beta+X_RE[,,theta]%*%Eta_tmp
    Yrep = rpois(n_Y,exp(rnorm(n_Y,Mu_tmp,jags_fit$BUGSoutput$sims.matrix[i,"tau_iid"]^-0.5))) 
    Mu=X%*%Beta
    Mu_tmp= exp(Mu + 0.5*(1/jags_fit$BUGSoutput$sims.matrix[i,"tau_iid"]+1/jags_fit$BUGSoutput$sims.matrix[i,"tau_sp"])) #includes lognormal bias correction
    ft=sum((sqrt(Yrep)-sqrt(Mu_tmp))^2)
    ft_y=sum((sqrt(Y)-sqrt(Mu_tmp))^2)  
    p_omni_mixed = p_omni_mixed + (ft_y<ft)
  }
  Results[isim,"poisSpat_postp_ft_mixed"]=p_omni_mixed/n_mcmc
  save(Results,file="Result.Rda")
}


##### plot Results
load('Result.Rda')
apply(Results,2,'mean')

Plot_data <- data.frame(matrix(NA,22000,4))
colnames(Plot_data)=c("Sim.num","Est.mod","Pval.type","Pval")
Plot_data[,"Sim.num"]=rep(c(1:1000),22)
  
j=1
for(i in 1:33){
  Places=(j-1)*1000+c(1:1000)
  Str = unlist(strsplit(colnames(Results)[i],'_'))
  n_str = length(Str)
  Plot_data[Places,"Pval"]=Results[,i]
  Plot_data[Places,"Sim.num"]=c(1:1000)
  if(Str[2]=="postp" & Str[n_str]=="chisq")Plot_data[Places,"Pval.type"]='PP.ChiSq'
  if(Str[2]=="postp" & Str[n_str]=="ft")Plot_data[Places,"Pval.type"]='PP.FT'
  if(Str[2]=="postp" & Str[n_str]=="dev1")Plot_data[Places,"Pval.type"]='PP.Dev.Pois'
  if(Str[2]=="postp" & Str[n_str]=="dev2")Plot_data[Places,"Pval.type"]='PP.Dev.Gauss'
  if(Str[2]=="postp" & Str[n_str]=="tail")Plot_data[Places,"Pval.type"]='PP.Tail'
  #if(Str[2]=="postp" & Str[n_str]=="0")Plot_data[Places,"Pval.type"]='PP.ZeroInfl'
  if(Str[n_str]=="moran")Plot_data[Places,"Pval.type"]='Moran'
  if(Str[n_str]=="pivot1")Plot_data[Places,"Pval.type"]='Pivot.Pois'
  if(Str[n_str]=="pivot2")Plot_data[Places,"Pval.type"]='Pivot.Gauss'
  if(Str[n_str]=="sampled1")Plot_data[Places,"Pval.type"]='Sampled.Pois'
  if(Str[n_str]=="sampled2")Plot_data[Places,"Pval.type"]='Sampled.Gauss'
  if(Str[n_str]=="cv")Plot_data[Places,"Pval.type"]='Cross.val'
  if(Str[n_str]=="mixed")Plot_data[Places,"Pval.type"]='Mixed.FT'
  
  if(Str[1]=="poisIID" & Str[n_str]!="outliers" & Str[n_str]!='0'){
    Plot_data[Places,"Est.mod"] = "PoisMix"
    j=j+1
  }
  if(Str[1]=="poisSpat" & Str[n_str]!='0'){
    Plot_data[Places,"Est.mod"] = "PoisMixSp"
    j=j+1
  }
}

fmt <- function(){
  f <- function(x) {
    f= as.character(round(x,3))
    #if(x<.1 | x>0.9)f= as.character(round(x,1))
    #else f = as.character(round(x,3))
  }
  f
}

Plot_data[,"Pval.type"]=factor(Plot_data[,"Pval.type"])
Plot = ggplot(Plot_data)+geom_freqpoly(position="identity",size=1.1,aes(x=Pval,linetype=factor(Est.mod)),binwidth=0.1,alpha=.5)+facet_wrap(~Pval.type,scales='free_y')
Plot = Plot + theme(text=element_text(size=14),axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position='bottom')+xlab('P-value')+ylab('Relative frequency')+labs(linetype="Estimation model")
Plot = Plot + guides(linetype=guide_legend(keywidth=3,override.aes = list(alpha = 1)))+ scale_x_continuous(labels = fmt(),limits=c(0,1))
Plot

pdf('SpatRegSimPvals.pdf')
Plot
dev.off()
