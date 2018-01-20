load("post_hsmm.rda")
library(tidyverse)
library(cowplot)
library(coda)
Rcpp::sourceCpp('attendance_hsmm.cpp')

message("Making the plot")

obs_data = ddl$ehmat-1
obs_det = colSums(obs_data) %>% ifelse(.==0, NA, .)
n_occ = ncol(ddl$ehmat)

### Deviance functions

dev = function(par, data=NULL, model, m, p_mat){
  return(2*attend_likelihood(par=par, ddl=ddl, data=data, p_dm=p_mat, m=m, model=model))
}

### remove burn in (5000 iter)
for(i in 1:4){
  mcmc_list[[i]] = mcmc_list[[i]][-c(1:5000),]
}



plot(obs_det)
# par = colMeans(mcmc_list[[1]])

exp_det = attend_expected(fit_list[[1]]$par, ddl, p_dm, m=c(30,30,10,10), model=2) %>% 
  matrix(., ncol=n_occ, byrow=T) %>% colSums() %>% ifelse(.==0, NA, .)
lines(exp_det)

exp_det = attend_expected(fit_list[[2]]$par, ddl, p_dm_time, m=c(30,30,10,10), model=2) %>% 
  matrix(., ncol=n_occ, byrow=T) %>% colSums() %>% ifelse(.==0, NA, .)
lines(exp_det)


exp_det = attend_expected(fit_list[[3]]$par, ddl, p_dm, m=rep(1,4), model=1) %>% 
  matrix(., ncol=n_occ, byrow=T) %>% colSums() %>% ifelse(.==0, NA, .)
lines(exp_det)

exp_det = attend_expected(fit_list[[4]]$par, ddl, p_dm_time, m=rep(1,4), model=1) %>% 
  matrix(., ncol=n_occ, byrow=T) %>% colSums() %>% ifelse(.==0, NA, .)
lines(exp_det)

idx = seq(1, 50000, 10)
exp_det_hmm <- exp_det_hsmm <- matrix(NA, length(idx), 61)

for(i in 1:length(idx)){
  exp_det_hmm[i,] = attend_expected(mcmc_list[[3]][idx[i],], ddl, p_dm, m=rep(1,4), model=1) %>% 
    matrix(., ncol=n_occ, byrow=T) %>% colSums() 
  exp_det_hsmm[i,] = attend_expected(mcmc_list[[1]][idx[i],], ddl, p_dm, m=c(30,30,10,10), model=2) %>% 
    matrix(., ncol=n_occ, byrow=T) %>% colSums() 
}

hpd_exp_hmm = data.frame(HPDinterval(mcmc(exp_det_hmm), prob = .5), HPDinterval(mcmc(exp_det_hmm), prob = .9)) %>% 
  mutate_all(funs(ifelse(.==0, NA, .)))
hpd_exp_hsmm = data.frame(HPDinterval(mcmc(exp_det_hsmm), prob = .5), HPDinterval(mcmc(exp_det_hsmm), prob = .9)) %>% 
  mutate_all(funs(ifelse(.==0, NA, .)))


p1 = ggplot() + 
  xlab("Day") + ylab("Number of detected animals") +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=1:61), data=hpd_exp_hmm, fill="blue", alpha=0.5) + 
  geom_ribbon(aes(ymin=lower.1, ymax=upper.1, x=1:61), data=hpd_exp_hmm, fill="blue", alpha=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=1:61), data=hpd_exp_hsmm, fill="red", alpha=0.5) + 
  geom_ribbon(aes(ymin=lower.1, ymax=upper.1, x=1:61), data=hpd_exp_hsmm, fill="red", alpha=0.5) +
  geom_point(aes(x=1:61, y=obs_det))

save_plot("attendance_fit.pdf", p1, base_width=5.5)

