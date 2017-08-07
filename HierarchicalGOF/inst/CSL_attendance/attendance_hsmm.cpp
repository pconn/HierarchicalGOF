#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

arma::mat make_P_mat(const int& ch, const IntegerVector& m, const NumericVector p){
  int J = sum(m);
  arma::mat P_mat(J, J, fill::zeros);
  Rcpp::IntegerVector idx_end = Rcpp::cumsum(m);
  idx_end = idx_end - 1;
  Rcpp::IntegerVector idx_st = idx_end - m + 1;
  arma::vec pdiag(J, fill::zeros);
  for(int i=0; i<m.size(); i++) pdiag(span(idx_st(i),idx_end(i))) += R::dbinom(ch, 1, p(i), false);
  P_mat.diag() = pdiag;
  return P_mat;
}

Rcpp::List make_c_list(const NumericVector& par, const IntegerVector& m, const int& model){
  IntegerVector rvec;
  NumericVector cvec;
  NumericVector Fvec;
  NumericVector dvec;
  Rcpp::List clist(m.size());
  for(int i=0; i<m.size(); i++){
    rvec = seq_len(m(i));
    // add different models here
    if(model==1){
      dvec = Rcpp::dgeom(rvec-1, R::plogis(par(i), 0, 1, 1, 0));
      Fvec = Rcpp::pgeom(rvec-2, R::plogis(par(i), 0, 1, 1, 0));
    } else if(model==2){
      dvec = Rcpp::dpois(rvec-1, exp(par(i)));
      Fvec = Rcpp::ppois(rvec-2, exp(par(i)));
    }
    cvec = dvec / (1.0 - Fvec);
    clist[i] = Rcpp::ifelse(Rcpp::is_infinite(cvec), 1, cvec);
  }
  return clist;
}

arma::mat make_Gamma(const mat& A, const List& clist, const IntegerVector& m){
  int J = sum(m);
  IntegerVector idx_end = cumsum(m);
  idx_end = idx_end - 1;
  IntegerVector idx_st = idx_end - m + 1;
  arma::mat Gamma(J,J,fill::zeros);
  arma::vec cvec;
  for(int i=0; i<m.size(); i++){
    cvec = as<arma::vec>(clist[i]);
    for(int j=0; j<m.size(); j++){
      if(i==j){
        if(m(i)>1){
          Gamma( span(idx_st(i),idx_end(i)), span(idx_st(i),idx_end(i)) ) = arma::diagmat(1-cvec.head(m(i)-1), 1);
        }
        Gamma(idx_end(i), idx_end(i)) = 1 - cvec(m(i)-1);
      } else {
        Gamma(span(idx_st(i),idx_end(i)), idx_st(j)) = A(i,j) * cvec;
      }
    }
  }
  
  return Gamma;
}

// [[Rcpp::export]]
Rcpp::List attendance_hsmm(
    const NumericVector& dt_beta,
    const NumericMatrix& pmat,
    const int& model,
    const IntegerVector& m,
    const NumericVector& delta_prior,
    const IntegerVector& id,
    const IntegerVector& ch,
    const bool& debug
){
  int I=id.size();
  int k=4;
  int J = sum(m);
  double ll = 0.0;
  double u = 0.0;
  arma::mat P_mat(J,J);
  arma::mat phi(I, J, fill::zeros);
  
  // initial state prior (fixed)
  arma::rowvec delta = as<arma::rowvec>(delta_prior);
  // for(int j=1; j<=m(0); j++) delta(j-1) = 1.0/(j+1);
  // delta(0) = 1;
  delta = arma::normalise(delta, 1);
  
  // State transitions are fixed for attendance model //
  arma:mat A(k,k, fill::zeros);
  A(0,1) = 1; A(1,2) = 1; A(2,3) = 1; A(3,2) = 1;
  
  // make Gamma mat... move inside loop for time inhomogeneous HSMM
  Rcpp::List clist = make_c_list(dt_beta, m, model);
  mat Gamma = make_Gamma(A, clist, m);
  
  // Begin foreward alg. 
  for(int i=0; i<I; i++){
    P_mat = make_P_mat(ch(i), m, pmat(i,_));
    if( i==0 || id(i)!=id(i-1) ){
      phi.row(i) = delta * P_mat;
    } else{
      // Gamma = make_Gamma(A, clist);
      phi.row(i) = phi.row(i-1) * Gamma * P_mat;
    }
    u = arma::sum(phi.row(i));
    ll += log(u);
    
    phi.row(i) = arma::normalise(phi.row(i),1);
  }
  if(debug) Rcout << "-lnl = " << -ll << endl;
  return Rcpp::List::create(
    Rcpp::Named("loglik") = ll
  );
}

/*** R
attend_likelihood = function(par, ddl, data=NULL, p_dm, m, model, debug=FALSE){
  if(is.null(data)){
    ch = as.vector(t(ddl$ehmat-1))
  } else{
    if(!all(dim(ddl$ehmat)==dim(data))) stop("Dimensions of data and ddl$ehmat are not equal!")
    ch = as.vector(t(data))
  }
  id = rep(1:nrow(ddl$ehmat), each=ncol(ddl$ehmat)) - 1
  dt_beta = par[1:4]
  det_beta = par[-c(1:4)]
  pmat = plogis(p_dm %*% det_beta)
  pmat = matrix(ifelse(is.na(ddl$p$fix), pmat, ddl$p$fix), ncol=4, byrow=T)
  delta_prior = c(1, rep(0,sum(m)-1))
  return(-attendance_hsmm(dt_beta, pmat, model=model, m=m, delta_prior=delta_prior,
                          id=id, ch=ch, debug=debug)$loglik)
}

make_p_dm = function(ddl, pformula,...){
  p_dm = model.matrix(pformula, ddl$p,...)
  tmp = p_dm[is.na(ddl$p$fix),]
  idx = (apply(tmp, 2, var)==0) & (colMeans(tmp)!=1)
  return(p_dm[,!idx])
}

fit_attendance = function(ddl, p_dm, m, model, model_inits,
                          ln_prior=NULL, debug=FALSE, ...){
  if(is.null(ln_prior)){
    obj_foo = function(par){
      out = attend_likelihood(par, ddl=ddl, p_dm=p_dm, m=m, model=model)
      if(debug) cat("-lnl = ", out, "\n")
      return(out)
    }
  } else{
    obj_foo = function(par){
      out = attend_likelihood(par, ddl=ddl, p_dm=p_dm, m=m, model=model) - 
        ln_prior(par)
      if(debug){
        cat("-ln post = ", out, "\n")
      }
      return(out)
    }
  }
  if(missing(model_inits)) model_inits = rep(0, 4+ncol(p_dm))
  mle = optim(fn=obj_foo, par = model_inits, ...)
}

mcmc_attendance = function(ddl, p_dm, m, model, model_inits,
                           ln_prior=NULL, debug=FALSE, iter, Sig){
  inits = model_inits
  par_stor = matrix(NA, ncol=length(inits), nrow=iter)
  par_stor[1,] = inits
  tune = 2.4^2/length(inits)
  if(missing(Sig)) Sig = 3*diag(length(inits)) #solve(post_mode_pois$hessian)
  if(is.null(ln_prior)) ln_prior = function(par) return(0)
  jump = rep(0, iter)
  # Start MCMC
  st = Sys.time()
  # set.seed(111)
  for(ii in 2:iter){
    par_prop = par_stor[ii-1,] + rmvnorm(1, sigma=tune*Sig)
    mhr = exp(
      -attend_likelihood(par_prop, ddl=ddl, p_dm=p_dm, m=m, model=model) + 
        ln_prior(par_prop) +
        attend_likelihood(par_stor[ii-1,], ddl=ddl, p_dm=p_dm, m=m, model=model) - 
        ln_prior(par_stor[ii-1,])
    )
    if(runif(1,0,1)<mhr){
      par_stor[ii,] = par_prop
      jump[ii] = 1
    } else{
      par_stor[ii,] = par_stor[ii-1,]
    }
    # cat(mhr, "(", jump[ii],")")
    if(ii==10){
      tpi = (Sys.time() - st)/9
      message(paste0("Time of completion: ", st + iter*tpi))
    }
    if(ii %% 100 == 0 & ii>100){
      gamma = (100/ii)^0.25
      r = mean(jump[ii:(ii-100+1)])
      tune = exp(log(tune) + 2*(2^0.25)*gamma*(r-0.234))
      Sig_hat = var(par_stor[ii:(ii-100+1)])
      Sig = Sig + gamma*(Sig_hat- Sig)
      cat("iter", ii, "(accept", r, ") ")
    }
  }
  return(par_stor)
}
*/

/*** R
dlogitbeta = function(x, a, b, log=FALSE){
  if(!log){
    return(dbeta(plogis(x), a, b)*plogis(x)/(1+exp(x)))
  } else{
    return(dbeta(plogis(x), a, b, log=TRUE) + plogis(x, log=TRUE) - log((1+exp(x))))
  }
}
*/

// // [[Rcpp::export]]
// Rcpp::List attendance_decode(
//     const NumericVector& dtpar,
//     const NumericVector& detpar,
//     const int& model,
//     const IntegerVector& m,
//     const NumericVector& delta_prior,
//     const IntegerVector& id,
//     const NumericMatrix& fixmat,
//     const IntegerVector& ch
// ){
//   int I=id.size();
//   int k=4;
//   int J = sum(m);
//   double ll = 0.0;
//   double u = 0.0;
//   arma::mat P_mat(J,J);
//   arma::mat xi(I, J, fill::zeros);
//   IntegerVector state(I);
//   Rcpp::NumericVector p_par = Rcpp::plogis(detpar);
//   Rcpp::NumericVector p(k);
//   
//   // initial state prior (fixed)
//   arma::rowvec delta = as<arma::rowvec>(delta_prior);
//   delta = arma::normalise(delta, 1);
//   
//   // State transitions are fixed for attendance model //
//   arma:mat A(k,k, fill::zeros);
//   A(0,1) = 1; A(1,2) = 1; A(2,3) = 1; A(3,2) = 1;
//   
//   // make Gamma mat... move inside loop for time inhomogeneous HSMM
//   Rcpp::List clist = make_c_list(dtpar, m, model);
//   mat Gamma = make_Gamma(A, clist, m);
//   
//   // Begin foreward alg. 
//   for(int i=0; i<I; i++){
//     p = ifelse(is_na(fixmat(i,_)), p_par, fixmat(i,_));
//     P_mat = make_P_mat(ch(i), m, p);
//     if( i==0 || id(i)!=id(i-1) ){
//       xi.row(i) = delta * P_mat;
//     } else{
//       // Gamma = make_Gamma(A, clist);
//       for(int j=0; j<J; j++) xi(i,j) = ((xi.row(i-1) % Gamma.col(j).t()).max()) * P_mat(j,j);
//     }
//   }
//   
//   state(I-1) = xi.row(I-1).index_max();
//   for(int i=I-2; i>=0; i--){
//     if(id(i)!=id(i+1)){
//       state(i) = xi.row(i).index_max();
//     } else {
//       state(i) = (xi.row(i) % Gamma.col(state(i+1)).t()).index_max(); 
//     }
//   }
//   return Rcpp::List::create(
//     Rcpp::Named("state") = state
//   );
// }
// 
// /*** R
// attend_decode = function(par, ddl=NULL, data=NULL, m, model){
//   if(is.null(data)){
//     if(is.null(ddl)) stop("Data must be provided in 'ddl' or 'data' arguments")
//     data = ddl$ehmat-1
//   } 
//   ch = as.vector(t(data))
//   id = as.vector(t(row(data))) - 1
//   if(!is.null(ddl)){
//     fixmat = matrix(ddl$p$fix, ncol=4, byrow=T)
//   } else fixmat = matrix(NA, nrow=nrow(data), ncol=ncol(data))
//   dtpar = par[1:4]
//   detpar = c(par[5], rep(par[6],3))
//   delta_prior = as.numeric(c(1, rep(0,sum(m)-1)))
//   idx = attendance_decode(dtpar, detpar, model=model, m=m, delta_prior=delta_prior,
//                           id=id, fixmat=fixmat, ch=ch)$state + 1
//   return(
//     matrix(
//       rep(levels(ddl$p$stratum), m)[idx],
//       nrow=nrow(ddl$ehmat),
//       ncol=ncol(ddl$ehmat), byrow = T
//     )
//   )
// }
// */

// [[Rcpp::export]]
Rcpp::List attendance_expected(
    const NumericVector& dt_beta,
    const NumericMatrix& pmat,
    const int& model,
    const IntegerVector& m,
    const NumericVector& delta_prior,
    const IntegerVector& id
){
  int I=id.size();
  int J = sum(m);
  arma::rowvec expected(I);
  arma::rowvec state(J);
  arma::vec P_vec(J);
  
  arma::rowvec delta = as<arma::rowvec>(delta_prior);
  delta = arma::normalise(delta, 1);
  
  arma:mat A(4,4, fill::zeros);
  A(0,1) = 1; A(1,2) = 1; A(2,3) = 1; A(3,2) = 1;
  
  Rcpp::List clist = make_c_list(dt_beta, m, model);
  mat Gamma = make_Gamma(A, clist, m);
  
  for(int i=0; i<I; i++){
    P_vec = make_P_mat(1, m, pmat(i,_)).diag();
    if( i==0 || id(i)!=id(i-1) ){
      state = delta;
    } else{
      state = state * Gamma;
    }
    expected(i) = as_scalar(state * P_vec);
  }
  return Rcpp::List::create(
    Rcpp::Named("expected") = expected
  );
}

/*** R
attend_expected = function(par, ddl, p_dm, m, model){
  dt_beta = par[1:4]
  det_beta = par[-c(1:4)]
  pmat = plogis(p_dm %*% det_beta)
  pmat = matrix(ifelse(is.na(ddl$p$fix), pmat, ddl$p$fix), ncol=4, byrow=T)
  id = as.integer(matrix(ddl$p$id, ncol=4, byrow = T)[,1])
  delta_prior = c(1, rep(0,sum(m)-1))
  out = attendance_expected(dt_beta, pmat,model=model, m=m, delta_prior, id=id)$expected
  return(as.vector(out))
}
*/

int sample_du(arma::rowvec ppp){
  arma::rowvec cdf = cumsum(arma::normalise(ppp,1));
  double U = Rcpp::as<double>(Rcpp::runif(1));
  int out = 1;
  if(U<= cdf[0]) return(out);
  else
  {
    for(int i=1; i<ppp.n_elem; i++){ 
      if(U <= cdf[i]){
        out = i+1;
        return(out);
      }
    }
    return(out);
  }
}


// [[Rcpp::export]]
Rcpp::List attendance_sim_data(
    const NumericVector& dt_beta,
    const NumericMatrix& pmat,
    const int& model,
    const IntegerVector& m,
    const NumericVector& delta_prior,
    const IntegerVector& id
){
  int I=id.size();
  int J = sum(m);
  IntegerVector state_rel(I);
  NumericVector rep_data(I);  
  arma::rowvec state(J);
  arma::vec P_vec(J);
  IntegerVector state_vals = seq_len(4);
  
  arma::rowvec delta = as<arma::rowvec>(delta_prior);
  delta = arma::normalise(delta, 1);
  
  arma:mat A(4,4, fill::zeros);
  A(0,1) = 1; A(1,2) = 1; A(2,3) = 1; A(3,2) = 1;
  
  Rcpp::List clist = make_c_list(dt_beta, m, model);
  mat Gamma = make_Gamma(A, clist, m);
  for(int i=0; i<I; i++){
    P_vec = make_P_mat(1, m, pmat(i,_)).diag();
    if( i==0 || id(i)!=id(i-1) ){
      state_rel(i) = sample_du(delta);
    } else{
      state_rel(i) = sample_du(Gamma.row(state_rel(i-1)-1));
    }
    rep_data(i) = R::rbinom(1, P_vec(state_rel(i)-1));
  }
  return Rcpp::List::create(
    Rcpp::Named("state") = state_rel,
    Rcpp::Named("data") = rep_data
  );
}

/*** R
attend_sim_data = function(par, ddl, p_dm, m, model, state=FALSE){
  dt_beta = par[1:4]
  det_beta = par[-c(1:4)]
  pmat = plogis(p_dm %*% det_beta)
  pmat = matrix(ifelse(is.na(ddl$p$fix), pmat, ddl$p$fix), ncol=4, byrow=T)
  id = as.integer(matrix(ddl$p$id, ncol=4, byrow = T)[,1])
  delta_prior = c(1, rep(0,sum(m)-1))
  out =  attendance_sim_data(dt_beta, pmat,model=model, m=m, delta_prior, id=id)
}
*/

