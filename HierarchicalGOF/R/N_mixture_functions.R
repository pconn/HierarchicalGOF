########################################################
### Functions for fitting N-mixture models
########################################################
###
### Written by: Perry Williams
###
### Modifications:
### 18 April: Wrote initial .R file
###

#' Title MCMC algorithm for fitting an N-mixture model with covariates in parallel.
#'
#' @param Y An n*J matrix of n sites and J repeat visits to each site
#' @param X Covariate matrix for abundance
#' @param W Covariate matrix for detection probability
#' @param n.iter Number of MCMC iterations
#' @param checkpoint The number of MCMC iterations that are run between each save of a file containing the MCMC output
#' @param thin The thinning rate
#' @param name Character name of file output is saved after every checkpoint
#' @param starting.values MCMC starting values
#'
#' @return MCMC output
#' @export
run.Nmixmcmc.parallel=function(Y,
                               X,
                               W,
                               n.iter,
                               checkpoint,
                               thin,
                               name,
                               starting.values){

    ## Pre-calculations
    y..=sum(Y,na.rm=TRUE)
    yi.=apply(Y,1,sum,na.rm=TRUE)
    n=dim(Y)[1]
    J=apply(!is.na(Y),1,sum)

    ## Starting values
    N=starting.values[[3]]
    JN=J*N
    alpha=starting.values[[1]]
    beta=starting.values[[2]]
    lambda=exp(X%*%beta)
    p=exp(W%*%alpha)/(1+exp(W%*%alpha))
    accept.alpha=numeric(length(alpha))
    accept.beta=numeric(length(beta))
    k=1
    alpha.tune=rep(0.001,length(alpha))
    beta.tune=rep(0.001,length(beta))

    ## Containers
    N.tot.save=matrix(NA,n.iter/thin)
    alpha.save=matrix(NA,n.iter/thin,length(alpha))
    beta.save=matrix(NA,n.iter/thin,length(beta))
    T.mcmc.chi2.save=numeric(n.iter/thin)
    T.data.chi2.save=numeric(n.iter/thin)
    beta.tune.save=matrix(NA,n.iter/thin,length(beta))
    alpha.tune.save=matrix(NA,n.iter/thin,length(alpha))

    ##  Begin Gibbs Loop
    for(k in 1:n.iter){
        if (k%%10000==0)
            cat(k,"")

        ## Sample beta
        for(i in 1:length(beta)){
            beta.star=beta
            beta.star[i]=rnorm(1,beta[i],beta.tune[i])
            lambda.star=exp(X%*%beta.star)
            mh1=sum(dpois(N,lambda.star,log=TRUE))+
                sum(dnorm(beta.star,0,10,log=TRUE))
            mh2=sum(dpois(N,lambda,log=TRUE))+
                sum(dnorm(beta,0,10,log=TRUE))
            mh=min(1,exp(mh1-mh2))
            if(mh>runif(1)){
                beta=beta.star
                lambda=lambda.star
                accept.beta[i]=accept.beta[i]+1
            }
        }

        ## Sample N
        N.star=sample(-1:1,n,replace=TRUE)+N
        mh1=apply(dbinom(Y,N.star,p,log=TRUE),1,sum)+
            dpois(N.star,lambda,log=TRUE)
        mh2=apply(dbinom(Y,N,p,log=TRUE),1,sum,na.rm=TRUE)+
            dpois(N,lambda,log=TRUE)
        mh=exp(mh1-mh2)
        mh[is.na(mh)]=0
        rcut=runif(length(mh))
        N=ifelse(mh>rcut,N.star,N)
        N.tot=sum(N,na.rm=TRUE)
        JN=J*N

        ## Sample alpha
        for(i in 1:length(alpha)){
            alpha.star=alpha
            alpha.star[i]=rnorm(1,alpha[i],alpha.tune[i])
            p.star=exp(W%*%alpha.star)/(1+exp(W%*%alpha.star))
            mh1=sum(dbinom(Y,N,p.star,log=TRUE))+
                sum(dnorm(alpha.star,0,1.5,log=TRUE))
            mh2=sum(dbinom(Y,N,p,log=TRUE))+
                sum(dnorm(alpha,0,1.5,log=TRUE))
            mh=min(1,exp(mh1-mh2))
            if(mh>runif(1)){
                alpha=alpha.star
                p=p.star
                accept.alpha[i]=accept.alpha[i]+1
            }
        }

        ## ppd (Chi-squared omnibus test statistic)
        y.ppd=matrix(rbinom(length(Y),N,p),dim(Y)[1],dim(Y)[2])
        y.ppd[is.na(Y)]=NA
        Expected.Y=matrix(N*p,dim(Y)[1],dim(Y)[2])
        T.mcmc.chi2=sum(((y.ppd-Expected.Y)^2)/Expected.Y,na.rm=TRUE)
        T.data.chi2=sum(((Y-Expected.Y)^2)/Expected.Y,na.rm=TRUE)

        ## Save samples
        if (k%%thin==0){
            alpha.save[k/thin,]=alpha
            beta.save[k/thin,]=beta
            N.tot.save[k/thin,]=N.tot
            T.mcmc.chi2.save[k/thin]=T.mcmc.chi2
            T.data.chi2.save[k/thin]=T.data.chi2
            beta.tune=ifelse(accept.beta/k<0.35,beta.tune*0.8,beta.tune)
            beta.tune=ifelse(accept.beta/k>0.45,beta.tune*1.2,beta.tune)
            alpha.tune=ifelse(accept.alpha/k<0.35,alpha.tune*0.8,alpha.tune)
            alpha.tune=ifelse(accept.alpha/k>0.45,alpha.tune*1.2,alpha.tune)
            beta.tune.save[k/thin,]=beta.tune
            alpha.tune.save[k/thin,]=alpha.tune
        }
        if (k%%checkpoint==0){
            cat("\n")
            out=cbind(
                alpha.save,
                beta.save,
                N.tot.save,
                T.mcmc.chi2.save,
                T.data.chi2.save,
                beta.tune.save,
                alpha.tune.save
            )
            colnames(out)=c("alpha0","alpha1","beta0","beta1","N.tot",
                            "T.mcmc.chi2.save","T.data.chi2.save","beta0.tune","beta1.tune",
                            "alpha0.tune","alpha1.tune"
            )

            save(out,file=name)
        }
    }
    return(out)
}


#' Title MCMC algorithm for fitting simple N-mixture model with varying lambda[i].
#'
#' @param Y An n*J matrix of n sites and J repeat visits to each site
#' @param q.p shape parameter for beta distribution prior on p
#' @param r.p shape parameter for beta distribution prior on p
#' @param alpha shape parameter for gamma distribution prior on lambda
#' @param beta rate parameter for gamma distribution prior on lambda
#' @param n.iter Number of MCMC iterations
#' @param checkpoint The number of MCMC iterations that are run between each save of a file containing the MCMC output
#' @param name Character name of file output is saved after every checkpoint
#' @param thin The thinning rate
#'
#' @return MCMC output
#' @export
Nmixmcmc=function(Y,q.p,r.p,alpha,beta,n.iter,checkpoint,name,thin){

    ## Packages and subroutines

    ## Dimensions
    n=dim(Y)[1]
    J=apply(!is.na(Y),1,sum)

    ## Containers
    p.save=matrix(NA,n.iter/thin,1)
    N.save=matrix(NA,n.iter/thin,n)
    lam.save=matrix(NA,n.iter/thin,n)
    N.total.save=numeric(n.iter/thin)
    T.mcmc.chi2.save=numeric(n.iter/thin)
    T.data.chi2.save=numeric(n.iter/thin)
    sev.save=numeric(n.iter/thin)

    ## Pre-calculations
    y..=sum(Y,na.rm=TRUE)
    yi.=apply(Y,1,sum,na.rm=TRUE)

    ## Starting values
    lam=apply(Y,1,mean,na.rm=TRUE)
    N=apply(Y,1,max,na.rm=TRUE)+1
    p=.8
    k=1

    ##  Begin Gibbs Loop
    for (k in 1:n.iter){
        if (k%%10000==0)
            cat(k,"")

        ## Sample N (primary)
        N.star=sample(-1:1,n,replace=TRUE)+N
        mh1=apply(dbinom(Y,N.star,p,log=TRUE),1,sum,na.rm=TRUE)+
            dpois(N.star,lam,log=TRUE)
        mh2=apply(dbinom(Y,N,p,log=TRUE),1,sum,na.rm=TRUE)+
            dpois(N,lam,log=TRUE)
        mh=exp(mh1-mh2)
        mh[is.na(mh)]=0
        rcut=runif(length(mh))
        N=ifelse(mh>rcut,N.star,N)
        JN=J*N

        ## Sample multiple lambdas
        lam=rgamma(n,shape=N+alpha,rate=beta+1)

        ##  Sample p
        p=rbeta(1,y..+q.p,sum(JN-yi.)+r.p)

        ## Chi-squared GOF
        y.ppd=matrix(rbinom(length(Y),N,p),dim(Y)[1],dim(Y)[2])
        y.ppd[is.na(Y)]=NA
        Expected.Y=matrix(N*p,dim(Y)[1],dim(Y)[2])
        T.mcmc.chi2=sum(((y.ppd-Expected.Y)^2)/Expected.Y,na.rm=TRUE)
        T.data.chi2=sum(((Y-Expected.Y)^2)/Expected.Y,na.rm=TRUE)

        ## Sum of empirical variance
        sev=sum(apply(y.ppd,1,var,na.rm=TRUE))

        ##  Save Samples
        if (k%%thin==0){
            p.save[k/thin]=p
            N.save[k/thin,]=N
            lam.save[k/thin,]=lam
            N.total.save[k/thin]=sum(JN/J,na.rm=TRUE)
            T.mcmc.chi2.save[k/thin]=T.mcmc.chi2
            T.data.chi2.save[k/thin]=T.data.chi2
            sev.save[k/thin]=sev
        }
        if (k%%checkpoint==0){
            ## cat("\n")
            out=list(
                p.save,
                N.save,
                lam.save,
                N.total.save,
                T.mcmc.chi2.save,
                T.data.chi2.save,
                Y,
                sev.save
            )
            save(out,file=name)
        }
    }
}



#' Title Wrapper for parallel processing
#'
#' @param models Vector of models to be run
#' @param Y.list List of n*J matrices of n sites and J repeat visits to each site
#' @param X Covariate matrix for abundance
#' @param W Covariate matrix for detection probability
#' @param n.iter Number of MCMC iterations
#' @param checkpoint The number of MCMC iterations that are run between each save of a file containing the MCMC output
#' @param thin The thinning rate
#' @param name.l List of character name of file output that is saved after every checkpoint
#' @param starting.values A list of MCMC starting values
#'
#' @return MCMC output
#' @export
run.chain.2pl.list=function(models,
                            Y.list,
                            X,
                            W,
                            n.iter,
                            checkpoint,
                            thin,
                            name.l,
                            starting.values
){
    chain.list=mclapply(models,
                        function(m){
                            ## Set the seed for this core
                            this.Y=Y.list[[m]]
                            this.start=starting.values.l[[m]]
                            this.name=name.l[[m]]
                            ## Run the chain on this core
                            this.chain=run.Nmixmcmc.parallel(this.Y,
                                                             X,
                                                             W,
                                                             n.iter,
                                                             checkpoint,
                                                             thin,
                                                             this.name,
                                                             this.start
                            )

                            ## Return the chain from this core
                            return(this.chain)
                        },
                        mc.cores=min(length(models),detectCores())
    )

    ## Save the initial random seed as the name of the chain
    names(chain.list)=models
    return(chain.list)
}






