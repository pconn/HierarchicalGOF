##### Run spatial negative binomial regression model using JAGS

setwd('c:/users/paul.conn/git/hierarchicalGOF/hierarchicalGOF/inst')
library(MASS)
library(fields)
library(mvtnorm)
library(spatstat)
library(rjags)
library(R2jags)
library(ape)
library(ggplot2)
library(RColorBrewer)
set.factory("bugs::Conjugate",FALSE,type="sampler")

#' plot a map of estimated abundance or related quantity.  
#' @param N quantity to plot
#' @param Coords data.frame holding x and y coordinates of each prediction cell
#' @param highlight [optional] A vector of cells to highlight.  Default = NULL
#' @param myPalette A prespecified palette from RColorBrewer (default is YlOrBr(9) )
#' @param leg.title Title of the legend
#' @return ggplot2 map of predictions on a grid
#' @export 
#' @keywords 
#' @author Paul Conn \email{paul.conn@@noaa.gov}
plot.prediction.map<-function(N,Coords,highlight=NULL,myPalette=NULL,leg.title="Abundance"){
  #require(rgeos)
  #require(ggplot2)
  #library(RColorBrewer)
  if(is.null(myPalette))myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  if(is.null(highlight)==FALSE){
    midpoints=Coords[highlight,]
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N
  Cur.df=cbind(Coords,Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  Grid.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+Grid.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
  }
  p1
}

#######    1) simulate some data   ########
set.seed(12345)
n_Y=200  #number of sample locations
Knots=expand.grid(c(0:8),c(0:8))/2
n_k=nrow(Knots)
Coords=matrix(runif(n_Y*2,0.5,3.5),n_Y,2)
incr=(max(Knots)-min(Knots))/99
Coords_plot=expand.grid(x=c(0:99)*incr,y=c(c(0:99)*incr))  #everything with _plot after it is for a 100x100 surface tessalation
#construct SpatialPolygonsDataFrame for use with plotting function
Dpk=rdist(Coords_plot,Knots)
theta_true=2 
tau_true=1
Dyy=rdist(Coords,Coords)
Dyk=rdist(Coords,Knots)
Dkk=rdist(Knots,Knots)
Cov_kk=Exp.cov(Knots,theta=theta_true,distMat=Dkk)/tau_true
Cov_yk=Exp.cov(Coords,Knots,theta=theta_true,distMat=Dyk)/tau_true
Cov_pk=Exp.cov(Coords_plot,Knots,theta=theta_true,distMat=Dpk)/tau_true
X_RE = Cov_yk %*% solve(Cov_kk)
X_RE_plot = Cov_pk %*% solve(Cov_kk)
Matern_points = rMatClust(kappa=12,r=.25,mu=100)  #generate some points from a matern process to construct a hypothetical covariate
Matern_xy=cbind(Matern_points$x,Matern_points$y)*max(Knots[,1])
Matern_kde=kde2d(Matern_xy[,1],Matern_xy[,2],100,h=0.5,lims=c(0,max(Knots[,1]),0,max(Knots[,2])))
Wts=1/(rdist(Coords,expand.grid(x=Matern_kde$x,y=Matern_kde$y)))^2
Wts=Wts*1/rowSums(Wts)
Pred_cov =  Wts %*% matrix(Matern_kde$z,ncol=1) #hypothetical habitat covariate
denom=mean(Pred_cov)
Pred_cov=Pred_cov/denom  #stan Wts %*% as.vector(Matern_kde$z)dardize to have a mean on 1.0
Pred_cov_plot=as.vector(Matern_kde$z)/denom #have same standardization as for covariates used for fitting 

#plot.prediction.map(N=as.vector(Matern_kde[[3]]),Coords=data.frame(Coords_plot),myPalette=colorRampPalette(brewer.pal(9, "Greens")),highlight=NULL,leg.title="Covariate")


#pdf("Sim_data_covs.pdf")
cov_map = plot.prediction.map(N=as.vector(Matern_kde[[3]]),Coords=data.frame(Coords_plot),myPalette=colorRampPalette(brewer.pal(9, "Greens")),highlight=NULL,leg.title="Covariate")
#dev.off()


X=matrix(1,nrow=n_Y,ncol=2)  #design matrix for fixed effects
X[,2]=Pred_cov  # linear term only; quadratic appears to result in a multimodal posterior (more evident with larger data set)
X_plot=matrix(1,nrow=10000,ncol=2) 
X_plot[,2]=Pred_cov_plot
Beta_true=c(2,0.75)
Eta_k_true=matrix(rmvnorm(1,rep(0,n_k),Cov_kk),ncol=1)
Eta_true=X_RE%*%Eta_k_true
Eta_plot = X_RE_plot%*%Eta_k_true
Mu_plot = exp(X_plot %*% Beta_true+X_RE_plot%*%Eta_k_true)
Fixed_effect = X_plot%*%Beta_true
Sim_lambda=rgamma(n_Y,1,1/exp(X%*%Beta_true+Eta_true))
Y=rpois(n_Y,Sim_lambda)
hist(Y)
re_map = plot.prediction.map(N=as.vector(Eta_plot),Coords=data.frame(Coords_plot),myPalette=colorRampPalette(brewer.pal(9, "Blues")),highlight=NULL,leg.title="RE Value")
Mu_map = plot.prediction.map(N=as.vector(Mu_plot),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Expectation")
Pts_df = data.frame("Easting"=as.vector(Coords[,1]),"Northing"=as.vector(Coords[,2]),"Abundance"=rep(0,n_Y))
Mu_map=Mu_map + geom_point(data=Pts_df,aes(x=Easting,y=Northing))
Knots_df = data.frame("Easting"=as.vector(Knots[,1]),"Northing"=as.vector(Knots[,2]),"Abundance"=rep(0,n_k))
Mu_map=Mu_map + geom_point(data=Knots_df,aes(x=Easting,y=Northing),size=5,shape=1)
Mu_map

#MCMC options
n_iter = 2000
n_burnin = 1000  
n_thin = 1
n_chains=3  
n_mcmc=n_chains*(n_iter-n_burnin)/n_thin #total number of recorded MCMC iterations

#2) jags models
neg.bin.no.spat <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dnegbin(p[i],r)
    Y_sim[i] ~ dnegbin(p[i],r) #posterior predictions
    p[i] <- r/(r+Lambda[i])
    log(Lambda[i]) <- XB[i]
    log_dev_pred_i[i]<- loggam(Y_sim[i]+r)-loggam(r)-loggam(Y_sim[i]+1) + r*log(p[i]) + Y_sim[i]*log(1-p[i]) #for deviance-based omnibus Bayesian p-value
  }
  XB[1:n_Y] <- X[,] %*% Beta[1:n_B]
  # Priors
  r ~ dgamma(1,0.01)
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  dev_pred <- -2*sum(log_dev_pred_i)
}
n_B=ncol(X)
n_p=nrow(X_plot)
jags_data = list("n_Y","n_B","Y","X")
#jags_params = c("Beta","psi","Lambda","Y_sim")
jags_params = c("r","Beta","Y_sim")
jags_save =c("r","Beta","Y_sim","Lambda","deviance","dev_pred")
jags.inits = function(){
  list("Beta"=rnorm(n_B),"r"=rgamma(1,1,1))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,
                n.iter=n_iter,
                model.file=neg.bin.no.spat,
                DIC=FALSE,
                parameters.to.save=jags_save,
                n.chains=n_chains,
                n.burnin=n_burnin,
                n.thin=n_thin,
                working.directory=getwd())
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[1]"])
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[2]"])
plot(jags_fit$BUGSoutput$sims.matrix[,"r"])

n_samples=1000
Lambda_pred=rep(0,10000)
Sampled=sample(c(1:n_mcmc),n_samples)
Beta_cols=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
for(isamp in 1:n_samples){
  Lambda_pred = Lambda_pred+exp(X_plot %*% jags_fit$BUGSoutput$sims.matrix[Sampled[isamp],Beta_cols])
}
Lambda_pred=Lambda_pred/n_samples
Mu_negbin_noRE_map = plot.prediction.map(N=as.vector(Lambda_pred),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Estimate")

##### calculate Moran's I for residuals, Chi-square omnibus test, Targeted skewness checks, Deviance omnibus test
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
  p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]==0)) 
}
p_omni=p_omni/n_mcmc
p_tail = p_tail/n_mcmc
p_0 = p_0/n_mcmc
p_ft = p_ft/n_mcmc
p_deviance = sum(jags_fit$BUGSoutput$sims.matrix[,"deviance"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred"])/n_mcmc
p_moran=sum(MoranI<0)/n_mcmc

neg.bin.no.spat.latent <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Y_sim[i] ~ dpois(Lambda_pred[i]) #posterior predictions
    Lambda[i] ~ dgamma(psi,Par2[i])
    Lambda_pred[i] ~ dgamma(psi,Par2[i])
    #Lambda[i] <- Mu[i]*Nu[i]
    #Nu[i] ~ dgamma(psi,psi)
    #Lambda_pred[i] <- Mu[i] * Nu_pred[i]
    #Nu_pred[i] ~ dgamma(psi,psi)
    log_dev_pred_i[i] <- Y_sim[i] * log(Lambda_pred[i])-Lambda_pred[i]-loggam(Y_sim[i]+1)
    log_dev_y_i[i] <- Y[i] * log(Lambda_pred[i])-Lambda_pred[i]-loggam(Y[i]+1)
    Par2[i] <- psi/Mu[i]
    log(Mu[i]) <- XB[i] 
  }
  XB[1:n_Y] <- X[,] %*% Beta[1:n_B]
  # Priors
  psi ~ dgamma(1,0.01)
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  dev_pred <- -2*sum(log_dev_pred_i)
  dev_y <- -2*sum(log_dev_y_i)
}
n_B=ncol(X)
jags_data = list("n_Y","n_B","Y","X")
jags_params = c("Beta","psi","Lambda","Y_sim")
#jags_params = c("Beta","psi","Nu","Y_sim","Nu_pred")
jags_save =c("Beta","psi","Y_sim","Lambda","Mu","deviance","dev_pred","dev_y")
jags.inits = function(){
  #list("Beta"=rnorm(n_B),"psi"=rgamma(1,1,1),"Nu"=rgamma(n_Y,1,1),"Nu_pred"=rgamma(n_Y,1,1))
  list("Beta"=rnorm(n_B),"psi"=rgamma(1,1,1),"Lambda"=exp(log(Y+0.1)+rnorm(n_Y)),"Lambda_pred"=exp(log(Y+0.1)+rnorm(n_Y)))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,
                n.iter=n_iter,
                model.file=neg.bin.no.spat.latent,
                DIC=FALSE,
                parameters.to.save=jags_save,
                n.chains=n_chains,
                n.burnin=n_burnin,
                n.thin=n_thin,
                working.directory=getwd())

##### calculate Moran's I for residuals, Chi-square omnibus test, Targeted skewness check
MoranI=rep(0,n_mcmc)
#select appropriate columns form jags MCMC output matrix
Cols_Y_sim=grep("Y_sim",colnames(jags_fit$BUGSoutput$sims.matrix))
Cols_Lambda=grep("Lambda",colnames(jags_fit$BUGSoutput$sims.matrix))
Cols_Mu=grep("Mu",colnames(jags_fit$BUGSoutput$sims.matrix))
Cols_Beta=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
Dists_inv=1/Dyy
diag(Dists_inv)=0
p_omni=p_tail=p_0=0
for(i in 1:n_mcmc){
  MoranI[i]=Moran.I(Y-jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu],Dists_inv)$observed
  #MoranI[i]=Moran.I(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-Y,Dists_inv)$observed
  chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])
  chisq_y=sum((Y-jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Mu])  
  #chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])
  #chisq_y=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda]-Y)^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])  
  p_omni = p_omni+(chisq_y<chisq)
  p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim],0.95))
  p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]==0)) 
}
p_omni=p_omni/n_mcmc
p_tail = p_tail/n_mcmc
p_0 = p_0/n_mcmc
p_moran=sum(MoranI<0)/n_mcmc


n_samples=1000
Lambda_pred=rep(0,10000)
Sampled=sample(c(1:n_mcmc),n_samples)
Beta_cols=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
for(isamp in 1:n_samples){
  Lambda_pred = Lambda_pred+exp(X_plot %*% jags_fit$BUGSoutput$sims.matrix[Sampled[isamp],Beta_cols])
}
Lambda_pred=Lambda_pred/n_samples
Mu_negbin_noRE_map2 = plot.prediction.map(N=as.vector(Lambda_pred),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Estimate")


neg.bin.model <- function(){
  for (i in 1:n_Y) {
    Y[i] ~ dpois(Lambda[i])
    Lambda[i] ~ dgamma(psi,Par2[i])
    Par2[i] <- psi/Mu[i]
    log(Mu[i]) <- XB[i] + Eta_y[i]
  }
  XB[1:n_Y] <- X[,] %*% Beta[1:n_B]
  Eta_y[1:n_Y] <- X_RE[1:n_Y,1:n_k,theta]%*%Eta_k[1:n_k]
  Eta_k[1:n_k] ~ dmnorm(Zero_k[1:n_k],Sigma_inv[1:n_k,1:n_k])
  Sigma_inv[1:n_k,1:n_k] <- tau*Q_k[1:n_k,1:n_k,theta]
  # Priors
  psi ~ dgamma(1,0.01)
  for(i in 1:n_B){
    Beta[i] ~ dnorm(0,0.01)
  }
  tau ~ dgamma(1.0,0.01)
  theta ~ dcat(P)
  # Functions of interest:
  # Posterior predictions
  for(i in 1:n_Y){
    Y_sim[i] ~ dpois(Lambda[i])
    Y_sim2[i] ~ dpois(Lambda_sim[i])
    Lambda_sim[i] <- exp(XB[i]+Eta_y_sim[i])
  }
  Eta_y_sim[1:n_Y] <- X_RE[1:n_Y,1:n_k,theta]%*%Eta_k_sim[1:n_k]
  Eta_k_sim[1:n_k] ~ dmnorm(Zero_k[1:n_k],Sigma_inv[1:n_k,1:n_k])
}

Zero_k = rep(0,n_k)
max_range=4
n_incr=10
incr=max_range/n_incr
Cor_k=Q_k=array(0,dim=c(n_k,n_k,n_incr))
X_RE=array(0,dim=c(n_Y,n_k,n_incr))
for(i in 1:n_incr){
  Cor_k[,1:n_k,i]=Exp.cov(Knots,theta=incr*i,distMat=Dkk)
  Q_k[,1:n_k,i]=solve(Cor_k[,1:n_k,i])
  X_RE[,1:n_k,i]=Exp.cov(Coords,Knots,theta=incr*i,distMat=Dyk)%*%Q_k[,1:n_k,i]
}
P=rep(incr,n_incr)*c(1:n_incr)
n_B=ncol(X)
jags_data = list("n_Y","n_B","Y","X","Zero_k","Q_k","X_RE","n_k","P")
jags_params = c("Beta","psi","Lambda","tau","theta","Eta_k")
jags_save =c("Beta","psi","tau","theta")
jags.inits = function(){
  list("Beta"=rnorm(n_B),"psi"=rgamma(1,1,1),"Lambda"=rnorm(n_Y,3.0,0.5),"tau"=runif(1,1,10),"theta"=round(runif(1,0.5,max_range+0.5)),"Eta_k"=rnorm(n_k))
}
jags_fit = jags(data=jags_data,
                inits=jags.inits,
                jags_params,n.iter=2000,
                model.file=neg.bin.model)

fit_mcmc <- as.mcmc(jags_fit)
#plot(fit_mcmc)
                #parameters.to.save=jags_save,
                #n.chains=3,
                #n.burnin=100,
                #n.thin=1,
                #DIC=FALSE)
plot(jags_fit$BUGSoutput$sims.matrix[,"theta"])
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[1]"])
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[2]"])
plot(jags_fit$BUGSoutput$sims.matrix[,"tau"])
plot(jags_fit$BUGSoutput$sims.matrix[,"psi"])


