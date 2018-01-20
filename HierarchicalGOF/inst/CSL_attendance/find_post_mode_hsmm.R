
### Need hsmm package from github user Jlaake.
# Install with the devtools package
devtools::install_github("jlaake/hsmm/hsmm")

library(hsmm)
library(tidyverse)
library(lubridate)
library(numDeriv)
Rcpp::sourceCpp(file.path(dir, 'attendance_hsmm.cpp'))

data(attendance)

attendance %>% filter(year %in% c(2002)) %>% droplevels() -> attendance
ddl = process.data(attendance,model="ATTEND",strata.labels=c("P","B","S","L")) %>% 
  make.design.data()
# set p=0 for "S" at sea 
ddl$p$fix=ifelse(ddl$p$stratum=="S",0,NA)
# set p=0 for occasions in which not a single animal was seen; presumably no effort
zc_resight = ddl$p %>% select(brand, year, Time) %>% distinct() %>% 
  mutate(ch = as.vector(t(ddl$ehmat-1))) %>% group_by(year, Time) %>% 
  dplyr::summarize(effort = as.numeric(sum(ch)>0)) %>% ungroup()
ddl$p %>% left_join(zc_resight) %>% mutate(fix = ifelse(effort==0, 0, fix)) -> ddl$p

ddl$p %>% mutate(birth = ifelse(stratum %in% c("B","S","L"), 1, 0) %>% factor()) -> ddl$p

# p_dm = make_p_dm(ddl,~birth+id, contrasts.arg=list(id=diag(nlevels(ddl$p$id))))
p_dm_time = make_p_dm(ddl, ~ birth + time - 1)
p_dm = make_p_dm(ddl, ~ birth-1)

ln_prior = function(model_par, prior_par){
  dlogitbeta(model_par[5],1,9,T) + dlogitbeta(model_par[6],9,1,T)
}

ln_prior_time = function(model_par, prior_par){
  dlogitbeta(model_par[5],1,9,T) + dlogitbeta(model_par[6],9,1,T) + 
    sum(-abs(model_par[-c(1:6)])/2)
}

fit_list = vector(mode = "list", 0)
### Poisson DT fit
# No time effects
message("Fitting HSMM without time effects")
fit_list$post_mode_pois =fit_attendance(ddl, p_dm=p_dm, m=c(30,30,10,10), model=2, 
                               model_inits=c(2,2,1,0,-3,4),
                               ln_prior=ln_prior,
                               method="BFGS",
                               control=list(trace=6, REPORT=10),
                               hessian=T
)

# time effects
message("Fitting HSMM *with* time effects")
pois_inits = c(fit_list$post_mode_pois$par, rep(0, ncol(p_dm_time)-2))
fit_list$post_mode_pois_time =fit_attendance(ddl, p_dm=p_dm_time, m=c(30,30,10,10), model=2, 
                                    model_inits=pois_inits,
                                    ln_prior=ln_prior_time,
                                    method="BFGS",
                                    control=list(trace=6, REPORT=10),
                                    hessian=T
)

### HMM fit-- geometric DT distribution
message("Fitting HMM without time effects")
fit_list$post_mode_geom=fit_attendance(ddl, p_dm=p_dm, m=rep(1,4), model=1, 
                              model_inits=c(-2,-2,-1,0,-3,4),
                              ln_prior=ln_prior,
                              method="BFGS",
                              control=list(trace=6, REPORT=10),
                              hessian=T
)

#time effects
message("Fitting HMM *with* time effects")
geom_inits = c(fit_list$post_mode_geom$par, rep(0, ncol(p_dm_time)-2))
fit_list$post_mode_geom_time =fit_attendance(ddl, p_dm=p_dm_time, m=rep(1,4), model=1, 
                                    model_inits=geom_inits,
                                    ln_prior=ln_prior_time,
                                    method="BFGS",
                                    control=list(trace=6, REPORT=10),
                                    hessian=T
)

message('Saving output')
save(list=ls(), file=paste0("hsmm_fit", ".rda"))
