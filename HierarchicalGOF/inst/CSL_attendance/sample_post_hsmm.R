# setwd("~/research/projects/sea_lion_analysis/ZC_attendance")
library(mvtnorm)
library(coda)
load("hsmm_fit.rda")
Rcpp::sourceCpp('attendance_hsmm.cpp')

mcmc_list = vector(mode="list",4)
names(mcmc_list)=names(fit_list)

mcmc_list[[1]] = mcmc_attendance(ddl, p_dm, m=c(30,30,10,10), model=2, 
                                 model_inits=fit_list[[1]]$par, 
                                 ln_prior=ln_prior, iter=55000, 
                                 Sig=solve(fit_list[[1]]$hessian)
                                 )
mcmc_list[[2]] = mcmc_attendance(ddl, p_dm_time, m=c(30,30,10,10), model=2, 
                                 model_inits=fit_list[[2]]$par, 
                                 ln_prior=ln_prior_time, iter=55000, 
                                 Sig=solve(fit_list[[2]]$hessian)
                                 )
mcmc_list[[3]] = mcmc_attendance(ddl, p_dm, m=rep(1,4), model=1, 
                                 model_inits=fit_list[[3]]$par, 
                                 ln_prior=ln_prior, iter=55000, 
                                 Sig=solve(fit_list[[3]]$hessian)
                                 )
mcmc_list[[4]] = mcmc_attendance(ddl, p_dm_time, m=rep(1,4), model=1, 
                                 model_inits=fit_list[[4]]$par, 
                                 ln_prior=ln_prior_time, iter=55000, 
                                 Sig=solve(fit_list[[4]]$hessian)
                                 )



save(list=ls(), file="post_hsmm.rda")




