#' function to calculate Bayesian Goodness-of-fit diagnostics using chi-squared distribution on inverse CDF
#' Algorithm due to Yuan and Johnson 2012 Biometrics 68:156-164
#' @param Y An (Number of observations x Number of mcmc samples) matrix of "Data" (these can actually be parameters as well; these are whatever the "leaf" level of the hierarchy being tested are).
#'          Each column should give a sample of Y from the posterior
#' @param Theta A (n_params x n_samples) matrix giving draws from the posterior distribution of parameter values
#' @param mean.fcn  Name of function for producing E(Y|Theta) - should be a function of Theta
#' @param var.fcn Name of function for producing Var(Y|Theta) - should be a function of Theta
#' @param cdf.fcn Name of function for evaluating cumulative distribution or mass function for [Y|Theta] - should be a function of two arguments, Y and Theta.
#' @param pmf.fcn Name of function for evaluating the pdf [Y|Theta] (only needed if a discrete distribution) - should be a function of two arguments, Y and Theta
#' @param DISCRETE If TRUE (default is false), denotes a discrete response; in such cases a pmf.fcn needs to be provided
#' @param K Number of bins for mean response.  Note if there is no variation in mean response, set to 1 (default)
#' @param L Number of bins per mean response bin
#' @return A list object with the following elements:
#'       \itemize{
#'         \item{"pval.k"}{K-dimensional vector giving p-values for each partition of mean response}
#'         \item{"pval.tot}{Overall p-value associated with the test}
#'         \item{"gof.plot}{A plot of average residual difference from theoretical quantiles}
#'        }
#' #' @import mvtnorm RandomFields spatstat gridExtra
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
chisq.cdf.test<-function(Y,Theta,K=1,L=5,mean.fcn,var.fcn,cdf.fcn,pmf.fcn=NULL,DISCRETE=FALSE){
  if(DISCRETE & is.null(pmf.fcn)==TRUE)cat("ERROR: A pmf must be specified for discrete distributions")
  if(is.null(cdf.fcn)==TRUE)cat("ERROR: A cdf must be specified")
  
  B=c(0:L)/L
  n=nrow(Y)
  m=ncol(Y)
  O=array(0,dim=c(m,K,L))
  N=matrix(0,m,K)
  Quants=rep(0,K+1)
  Quants[1]=-1e+10
  Quants[K+1]=-Quants[1]
  T_k=matrix(0,m,K)
  T_omni=rep(0,m)
  for(i in 1:m){
    Mu=mean.fcn(Theta[,i])
    V=var.fcn(Theta[,i])
    #W=(Y[,i]-Mu)*V^(-0.5)
    if(K>1){  
      Quants[2:K]=quantile(Mu,c(1:(K-1))/K)
      R=cut(Mu,Quants)
      N[i,]=tabulate(R) #keep all factors for tabulate
      R=as.numeric(R)
    }
    else{ #no need to do binning for mean response if only one group / no variation in mean response
      R=rep(1,n)
      N[i,]=n
    }
    if(DISCRETE)CDF=cdf.fcn(Y[,i]-1,Theta[,i])+runif(n)*pmf.fcn(Y[,i],Theta[,i])
    else CDF=cdf.fcn(Y[,i],Theta[,i])
    CDF_cat=as.numeric(cut(CDF,B))
    for(j in 1:n){
      O[i,R[j],CDF_cat[j]]=O[i,R[j],CDF_cat[j]]+1
    }
    for(k in 1:K)T_k[i,k]=sum((O[i,k,]-N[i,k]/L)^2*L/N[i,k]) #chi-square calculation
    T_omni[i]=sum(T_k[i,])
  }
  pval_k=1-pchisq(T_k,df=L-1)  #apply(T_k,1,pchisq,df=L-1)
  pval_omni=1-pchisq(T_omni,df=K*(L-1))
  return(list(pval_k=pval_k,pval_omni=pval_omni))
}




