load("post_hsmm.rda")
library(magrittr)
Rcpp::sourceCpp('attendance_hsmm.cpp')

obs_data = ddl$ehmat-1
obs_det = colSums(obs_data)
n_occ = ncol(ddl$ehmat)

### Deviance functions

dev = function(par, data=NULL, model, m, p_mat){
  return(2*attend_likelihood(par=par, ddl=ddl, data=data, p_dm=p_mat, m=m, model=model))
}

### remove burn in (5000 iter)
for(i in 1:4){
  mcmc_list[[i]] = mcmc_list[[i]][-c(1:5000),]
}

### p-value storage

pval_tukey_list = vector(mode="list", length=length(mcmc_list))
names(pval_tukey_list) = names(mcmc_list)

for(i in 1:length(mcmc_list)){

  set.seed(111)
  idx = seq(1, nrow(mcmc_list[[i]]), 10)
  reps = length(idx)
  
  if(i<3){
    m=c(30,30,10,10)
    model=2
    if(i==2) p_mat = p_dm_time else p_mat = p_dm
  } else {
    m=rep(1,4)
    model=1
    if(i==4) p_mat = p_dm_time else p_mat = p_dm
  }
  
pval_tukey_list[[i]] = rep(NA, reps)
  T_tukey_obs <- T_tukey_rep <- rep(NA, reps)
  st = Sys.time()
  for(j in 1:reps){
    par = mcmc_list[[i]][idx[j],]
    nd = attend_sim_data(par, ddl, p_mat, m=m, model=model)$data %>% matrix(., ncol=n_occ, byrow=T)
    exp_det = attend_expected(par, ddl, p_mat, m=m, model=model) %>% matrix(., ncol=n_occ, byrow=T) %>% colSums()
    T_tukey_obs[j] = sum((sqrt(obs_det) - sqrt(exp_det))^2)
    rep_det = colSums(nd)
    T_tukey_rep[j] = sum((sqrt(rep_det) - sqrt(exp_det))^2)
    if(j ==10){
      tpi = difftime(Sys.time(), st, units="secs")/j
      message("Start time: ", st)
      message("TPI = ", tpi)
      message("Reps = ", reps)
      message(paste0("Time of completion: ", st + reps*tpi))
    }
    if(j >=100 & j%%100==0) cat(j, "(", i, ") ")
  }
  pval_tukey_list[[i]] = as.numeric(T_tukey_rep > T_tukey_obs)
  
}

save(list=ls(), file="gof_hsmm.rda")
