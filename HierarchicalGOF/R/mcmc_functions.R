
#' JAGS function to conduct a Poisson GLM analysis with no spatial autocorrelation
#' or overdispersion 
#' @return RJags object that includes MCMC samples of Poisson regression parameters
#' @export
#' @keywords Jags, Poisson, MCMC
#' @author Paul B. Conn
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

#' JAGS function to conduct a Poisson GLM analysis with overdispersion
#' but no spatial autocorrelation
#' @return RJags object that includes MCMC samples of Poisson regression parameters and Gaussian precision for random effects on log scale
#' @export
#' @keywords Jags, Poisson, MCMC
#' @author Paul B. Conn
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


#' JAGS function to conduct a Poisson GLM analysis with overdispersion
#' but no spatial autocorrelation - this is a cleaner version for cross-validation
#' @return RJags object that includes MCMC samples of Poisson regression parameters and Gaussian precision for random effects on log scale
#' @export
#' @keywords Jags, Poisson, MCMC
#' @author Paul B. Conn
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

#' JAGS function to conduct a Poisson GLM analysis with overdispersion
#' and spatial autocorrelation 
#' @return RJags object that includes MCMC samples of Poisson regression parameters,
#' Gaussian precision for random effects on log scale, and Matern spatial autocorrelation parameters
#' @export
#' @keywords Jags, Poisson, MCMC
#' @author Paul B. Conn
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


#' Function to compute Matern correlation (from Hoeting et al. 2005 - Ecological Applications)
#' @param theta A vector containing range, smoothness, and nugget
#' @param D Matrix of distances between sites; typically a square matrix with zeros on main diagonal 
#' @return A matrix of correlations
#' @export
#' @keywords Matern, correlation
#' @author Andrew Merton
cor.matern <- function(theta,D) {  
  if (is.vector(D)) names(D) <- NULL
  if (is.matrix(D)) dimnames(D) <- list(NULL,NULL)
  th1 <- theta[1]                             # Range.
  th2 <- theta[2]                             # Smoothness.
  th3 <- ifelse(length(theta)==3,theta[3],0)  # Nugget effect.
  u <- 2*D*sqrt(th2)/th1
  u <- ifelse(u>0,(1-th3)*u^th2/(2^(th2-1)*gamma(th2))*besselK(u,th2),1)
  return(u)
}



