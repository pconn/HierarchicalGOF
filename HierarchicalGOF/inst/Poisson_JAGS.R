##### Run spatial negative binomial regression model using JAGS
### development code - to be repackaged for full simulation analysis in "run_count_simulations.R"

### I can't get spatial code in the negative binomial linear predictor without
# some weirdness happening (different chains coverging to different - and wrong - values, etc.)
# switching to overdispersed Poisson....


setwd('c:/users/paul.conn/git/hierarchicalGOF/hierarchicalGOF/inst')
#Sys.setenv(JAGS_HOME='C:/program files/JAGS/JAGS-4.0.0/x64/bin') 
library(MASS)
library(fields)
library(mvtnorm)
library(spatstat)
library(rjags)
library(R2jags)
library(ape)
library(ggplot2)
library(RColorBrewer)
#set.factory("bugs::Conjugate",FALSE,type="sampler")
source('../R/mcmc_functions.R')
source('../R/pivot_functions.R')

#' plot a map of estimated abundance or related quantity on a grid.  
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
tau_eta_true=1
tau_iid_true=5
Dyy=rdist(Coords,Coords)
Dyk=rdist(Coords,Knots)
Dkk=rdist(Knots,Knots)
Cov_kk=Exp.cov(Knots,theta=theta_true,distMat=Dkk)/tau_eta_true
Cov_yk=Exp.cov(Coords,Knots,theta=theta_true,distMat=Dyk)/tau_eta_true
Cov_pk=Exp.cov(Coords_plot,Knots,theta=theta_true,distMat=Dpk)/tau_eta_true
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
cov_map = plot.prediction.map(N=as.vector(Matern_kde[[3]]),Coords=data.frame(Coords_plot),myPalette=colorRampPalette(brewer.pal(9, "Greens")),highlight=NULL,leg.title="Value")
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
Sim_lambda=exp(rnorm(n_Y,X%*%Beta_true+Eta_true,sqrt(1/tau_iid_true)))
Y=rpois(n_Y,Sim_lambda)
hist(Y)
title_theme=theme(plot.title = element_text(hjust = 0,size=12))
re_map = plot.prediction.map(N=as.vector(Eta_plot),Coords=data.frame(Coords_plot),myPalette=colorRampPalette(brewer.pal(9, "Blues")),highlight=NULL,leg.title="Value")
Mu_map = plot.prediction.map(N=as.vector(Mu_plot),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Value")
Pts_df = data.frame("Easting"=as.vector(Coords[,1]),"Northing"=as.vector(Coords[,2]),"Abundance"=rep(0,n_Y))
Mu_map=Mu_map + geom_point(data=Pts_df,aes(x=Easting,y=Northing),size=0.5)
Knots_df = data.frame("Easting"=as.vector(Knots[,1]),"Northing"=as.vector(Knots[,2]),"Abundance"=rep(0,n_k))
Mu_map=Mu_map + geom_point(data=Knots_df,aes(x=Easting,y=Northing),size=1,shape=1)
Mu_map

colnames(Pts_df)[3]="Count"
Count_df = Pts_df
Count_df[,3]=Y
myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
Grid.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
Count_map = ggplot(Count_df)+geom_point(aes(x=Easting,y=Northing,colour=Count),size=1.5)+Grid.theme+scale_color_gradientn(colours=myPalette(100),name="Value")

library(gridExtra)
pdf(file="sim_count_maps.pdf")
grid.arrange(arrangeGrob(cov_map+ggtitle("A. Covariate")+title_theme,re_map+ggtitle("B. Spatial random effects")+title_theme,Mu_map+ggtitle("C. Expected abundance")+title_theme,Count_map+ggtitle("D. Simulated count")+title_theme,nrow=2))
dev.off()


#MCMC options
n_iter = 2000
n_burnin = 500  
n_thin = 1
n_chains=3  
n_mcmc=n_chains*(n_iter-n_burnin)/n_thin #total number of recorded MCMC iterations

#2) jags models
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
n_B=ncol(X)
n_p=nrow(X_plot)
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
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[1]"])
plot(jags_fit$BUGSoutput$sims.matrix[,"Beta[2]"])

n_samples=1000
Lambda_pred=rep(0,10000)
Sampled=sample(c(1:n_mcmc),n_samples)
Beta_cols=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
for(isamp in 1:n_samples){
  Lambda_pred = Lambda_pred+exp(X_plot %*% jags_fit$BUGSoutput$sims.matrix[Sampled[isamp],Beta_cols])
}
Lambda_pred=Lambda_pred/n_samples
Mu_pois_noRE_map = plot.prediction.map(N=as.vector(Lambda_pred),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Estimate")

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

Y_mat=matrix(Y,n_Y,n_mcmc)
mean.fcn.pois <- function(Theta){ #takes Mu [1:n_Y]
  return(Theta)
}
var.fcn.pois <- function(Theta){
  return(Theta)
}
Y_mat=matrix(Y,n_Y,n_iter)
p_pivot1=chisq.cdf.test(Y=Y_mat,Theta=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Lambda]),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)


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
n_B=ncol(X)
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
  #MoranI[i]=Moran.I(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-Y,Dists_inv)$observed
  chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]-Mu_tmp)^2/Mu_tmp)
  chisq_y=sum((Y-Mu_tmp)^2/Mu_tmp)  
  ft=sum((sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2])-sqrt(Mu_tmp))^2)
  ft_y=sum((sqrt(Y)-sqrt(Mu_tmp))^2)
  p_ft = p_ft + (ft_y<ft)
    #chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])
  #chisq_y=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda]-Y)^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])  
  p_omni = p_omni+(chisq_y<chisq)
  p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2],0.95))
  p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]==0)) 
}
p_omni=p_omni/n_mcmc
p_tail = p_tail/n_mcmc
p_0 = p_0/n_mcmc
p_moran=sum(MoranI<0)/n_mcmc
p_deviance1 = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y1"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred1"])/n_mcmc
p_deviance2 = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y2"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred2"])/n_mcmc
p_ft = p_ft/n_mcmc

Y_mat=matrix(Y,n_Y,n_mcmc)
#test poisson part of poisson-normal mixture
mean.fcn.pois <- function(Theta){ #takes Lambda
  return(Theta)
}
var.fcn.pois <- function(Theta){ #takes Lambda
  return(Theta)
}
my.ppois<-function(Y,Theta)ppois(Y,Theta)
my.dpois<-function(Y,Theta)dpois(Y,Theta)
p_pivot2a=chisq.cdf.test(Y=Y_mat,Theta=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Lambda]),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)

#test normal part of poisson-normal mixture using pivotal discrepancy
mean.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(Theta[1:n_par])
}
var.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(1/Theta[n_par])
}
my.pnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(pnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}
my.dnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(dnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}
p_pivot2b=chisq.cdf.test(Y=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu]),
          Theta=t(cbind(jags_fit$BUGSoutput$sims.matrix[,Cols_Mu],
          jags_fit$BUGSoutput$sims.matrix[,"tau_iid"])),
          mean.fcn=mean.fcn.normal,var.fcn=var.fcn.normal,
          cdf.fcn=my.pnorm,pmf.fcn=my.dnorm,K=5,L=5,DISCRETE=FALSE)

#plot normal CDF values for a given MCMC iteration to 
tmp_iter=2000
CDF_vals = pnorm(jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Nu],jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Mu],1/sqrt(jags_fit$BUGSoutput$sims.matrix[100,"tau_iid"]))
plot(jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Mu],CDF_vals,xlab="Linear predictor value",ylab="Gaussian cumulative density function",cex.lab=1.3,cex=1.3)
abline(v=quantile(jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Mu],c(0.2,0.4,0.6,0.8)))
abline(h=c(0.2,0.4,0.6,0.8),lty=2)
pdf("pivot_CDF_plot.pdf")
plot(jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Mu],CDF_vals,xlab="Linear predictor value",ylab="Gaussian cumulative density function",cex.lab=1.3,cex=1.3,col="darkgray")
abline(v=quantile(jags_fit$BUGSoutput$sims.matrix[tmp_iter,Cols_Mu],c(0.2,0.4,0.6,0.8)))
abline(h=c(0.2,0.4,0.6,0.8),lty=2)
dev.off()

n_samples=1000
Lambda_pred=rep(0,10000)
Sampled=sample(c(1:n_mcmc),n_samples)
Beta_cols=grep("Beta",colnames(jags_fit$BUGSoutput$sims.matrix))
for(isamp in 1:n_samples){
  Lambda_pred = Lambda_pred+exp(X_plot %*% jags_fit$BUGSoutput$sims.matrix[Sampled[isamp],Beta_cols])
}
Lambda_pred=Lambda_pred/n_samples
Mu_pois_nospat_map2 = plot.prediction.map(N=as.vector(Lambda_pred),Coords=data.frame(Coords_plot),highlight=NULL,leg.title="Estimate")

#conduct leave-k-out cross validation
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
    #list("Beta"=rnorm(n_B),"psi"=rgamma(1,1,1),"Nu"=rgamma(n_Y,1,1),"Nu_pred"=rgamma(n_Y,1,1))
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
p_cv=1-pchisq(T_cv,df=9)  

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
  
#  Eta_y_sim[1:n_Y] <- X_RE[1:n_Y,1:n_k,theta]%*%Eta_k_sim[1:n_k]
#  for(k in 1:n_k){
#    Eta_k_sim[k] ~ dnorm(Mu_RE[k],Tau_k[k])
#    Mu_RE[k] <- 0
#    Tau_k[k] <- 1
    #Mu_RE[k] <- t(Cor_partial[k,,theta]) %*% Q_k_partial[k,,,theta] %*% Eta_k[Partial_index[k,1:(n_k-1)]]
    #poop[k] <- tau / (1 - t(Cor_partial[k,,theta])%*%(Cor_partial[k,,theta]))
#  



Zero_k = rep(0,n_k)
#discretize covariance matrices and related quantities so all matrix inverse
#calculations can be done up front and provided as input to JAGS
max_range=4
n_incr=10
incr=max_range/n_incr
Cor_k=Q_k=array(0,dim=c(n_k,n_k,n_incr))
Q_k_partial = array(0,dim=c(n_k,n_k-1,n_k-1,n_incr)) #for partial posterior prediction
Cor_partial = array(0,dim=c(n_k,n_k-1,n_incr)) 
X_RE=array(0,dim=c(n_Y,n_k,n_incr))
Tmp=c(1:n_k)
#Partial_index = matrix(0,n_k,n_k-1)  #for jags, but trying to do partial predictions outside of jags now
for(i in 1:n_incr){
  Cor_k[,1:n_k,i]=Exp.cov(Knots,theta=incr*i,distMat=Dkk)
  Q_k[,1:n_k,i]=solve(Cor_k[,1:n_k,i])
  X_RE[,1:n_k,i]=Exp.cov(Coords,Knots,theta=incr*i,distMat=Dyk)%*%Q_k[,1:n_k,i]
  for(k in 1:n_k){
    Q_k_partial[k,,,i]=solve(Cor_k[Tmp[-k],Tmp[-k],i])
    Cor_partial[k,,i] = Cor_k[k,Tmp[-k],i]
  }
}
#for(k in 1:n_k)Partial_index[k,1:(n_k-1)]=Tmp[-k]
P=rep(incr,n_incr)*c(1:n_incr)
n_B=ncol(X)
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


##### calculate Moran's I for residuals, Chi-square omnibus test, Targeted skewness check
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
  #MoranI[i]=Moran.I(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-Y,Dists_inv)$observed
  chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]-Mu_tmp)^2/Mu_tmp)
  chisq_y=sum((Y-Mu_tmp)^2/Mu_tmp)  
  ft=sum((sqrt(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2])-sqrt(Mu_tmp))^2)
  ft_y=sum((sqrt(Y)-sqrt(Mu_tmp))^2)
  p_ft = p_ft + (ft_y<ft)
  #chisq=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim]-jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])
  #chisq_y=sum((jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda]-Y)^2/jags_fit$BUGSoutput$sims.matrix[i,Cols_Lambda])  
  p_omni = p_omni+(chisq_y<chisq)
  p_tail = p_tail + (y_quant95<quantile(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2],0.95))
  p_0 = p_0 + (y_0<sum(jags_fit$BUGSoutput$sims.matrix[i,Cols_Y_sim2]==0)) 
}
p_omni=p_omni/n_mcmc
p_tail = p_tail/n_mcmc
p_0 = p_0/n_mcmc
p_moran=sum(MoranI<0)/n_mcmc
p_deviance1 = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y1"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred1"])/n_mcmc
p_deviance2 = sum(jags_fit$BUGSoutput$sims.matrix[,"dev_y2"]<jags_fit$BUGSoutput$sims.matrix[,"dev_pred2"])/n_mcmc
p_ft = p_ft/n_mcmc

Y_mat=matrix(Y,n_Y,n_mcmc)
#test poisson part of poisson-normal mixture
mean.fcn.pois <- function(Theta){ #takes Lambda
  return(Theta)
}
var.fcn.pois <- function(Theta){ #takes Lambda
  return(Theta)
}
my.ppois<-function(Y,Theta)ppois(Y,Theta)
my.dpois<-function(Y,Theta)dpois(Y,Theta)
p_pivot3a=chisq.cdf.test(Y=Y_mat,Theta=t(Lambda),mean.fcn=mean.fcn.pois,var.fcn=var.fcn.pois,cdf.fcn=ppois,pmf.fcn=dpois,K=5,L=5,DISCRETE=TRUE)

#test normal part of poisson-normal mixture using pivotal discrepancy
mean.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(Theta[1:n_par])
}
var.fcn.normal <- function(Theta){ #takes Mu[1:n],tau
  n_par=length(Theta)
  return(1/Theta[n_par])
}
my.pnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(pnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}
my.dnorm<-function(Y,Theta){
  n_par=length(Theta)
  return(dnorm(Y,Theta[1:n_par],sqrt(1/Theta[n_par])))
}
p_pivot3b=chisq.cdf.test(Y=t(jags_fit$BUGSoutput$sims.matrix[,Cols_Nu]),
                         Theta=t(cbind(jags_fit$BUGSoutput$sims.matrix[,Cols_Mu],
                                       jags_fit$BUGSoutput$sims.matrix[,"tau_iid"])),
                         mean.fcn=mean.fcn.normal,var.fcn=var.fcn.normal,
                         cdf.fcn=my.pnorm,pmf.fcn=my.dnorm,K=5,L=5,DISCRETE=FALSE)

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
p_omni_mixed=p_omni_mixed/n_mcmc


