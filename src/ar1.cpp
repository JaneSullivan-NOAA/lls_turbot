#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template <class Type>
Type square(Type x){return x * x;}

template<class Type>
bool isNA(Type x){return R_IsNA(asDouble(x));}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // model dimensions
  DATA_IVECTOR(model_yrs);
  
  // survey biomass (absolute biomass gives scale to the model)
  DATA_MATRIX(obs);
  DATA_MATRIX(obs_cv);
  
  // data for one-step-ahead (OSA) residuals
  // DATA_VECTOR(obsvec); // vector of all observations for OSA residuals
  // DATA_VECTOR_INDICATOR(keep, obsvec); // for OSA residuals
  // DATA_IMATRIX(keep_biomass_obs); // indices for biomass survey obs, can loop years/survey strata with keep(keep_biomass_obs(i,j))
  // DATA_IMATRIX(keep_cpue_obs);
  // derived quantities - model dimensions

  
  // parameter section
  
  // PARAMETER(dummy);  // dummy var for troubleshooting
  PARAMETER_VECTOR(log_PE); // process errors
  PARAMETER_VECTOR(logit_rho); // AR-1 parameter

  // random effects of predicted biomass
  PARAMETER_MATRIX(log_pred);
  
  int nyrs = model_yrs.size();
  int n_strata = obs.cols();
  
  vector<Type> PE(n_strata);
  // for(int j = 0; j < n_strata; j++) PE(j) = exp(log_PE(j));
  PE = exp(log_PE);
  
  vector<Type> rho(n_strata);
  // for(int j = 0; j < n_strata; j++) rho(j) = -1 + 2 / (1 + exp(-logit_rho(j)));
  rho = -1 + 2 / (1 + exp(-logit_rho));
  
  // negative log likelihood
  vector<Type> jnll(2); // ar1, biomass obs
  jnll.setZero();
  Type nll = 0;
  


  // model predictions
  
  // predicted biomass and cpue on the natural and log scale
  matrix<Type> pred(nyrs, n_strata);
  pred.setZero();
  
  // derived quantities - observations
  matrix<Type> log_obs(nyrs, n_strata);
  log_obs = log(obs.array());
  matrix<Type> log_obs_sd(nyrs, n_strata);
  for(int i = 0; i < nyrs; i++) {
    for(int j = 0; j < n_strata; j++) {
      log_obs_sd(i,j) = obs_cv(i,j) * obs_cv(i,j) + Type(1.0);
      log_obs_sd(i,j) = sqrt(log(log_obs_sd(i,j)));
    }
  }

  // add wide prior for first predicted biomass, but only when computing osa
  // residuals
  // if(CppAD::Variable(keep.sum())){
  //   Type huge = 10;
  //   for(int j = 0; j < n_strata; j++) {
  //     jnll -= dnorm(pred(0, j), Type(0), huge, true);
  //   }
  // }
  
  // random effects contribution to likelihood
  // for(int i = 1; i < nyrs; i++) {
  //   for(int j = 0; j < n_strata; j++) {
  //     jnll(0) -= dnorm(log_pred(i-1,j), log_pred(i,j), exp(log_PE(j)), 1);
  //   }
  // }
  
  // random effects contribution to likelihood
  for(int j = 0; j < n_strata; j++) {
    jnll(0) += SCALE(AR1(rho(j)), PE(j))(log_pred.col(j)); 
  }
  // 
  // SIMULATE {
  //   AR1(rho_re).simulate(re);
  //   re *= sig_re;
  // }
  // 
  // likelihood for observation error
    for(int i = 0; i < nyrs; i++) {
    for(int j = 0; j < n_strata; j++) {
      
      if(obs(i,j) > 0) {
        jnll(1) -= dnorm(log_obs(i,j), log_pred(i,j), log_obs_sd(i,j), 1);
        // jnll(1) -= keep(keep_obs(i,j)) * dnorm(log_obs(i,j), log_pred(i,j), log_sd(i,j), 1);
        // jnll(1) -= keep.cdf_lower(keep_obs(i,j)) * log(squeeze(pnorm(log_obs(i,j), log_pred(i,j), log_obs_sd(i,j))));
        // jnll(1) -= keep.cdf_upper(keep_obs(i,j)) * log(1.0 - squeeze(pnorm(log_obs(i,j), log_pred(i,j), log_obs_sd(i,j))));
      }
      
    }
  }
   
   for(int i = 0; i < nyrs; i++) {
     for(int j = 0; j < n_strata; j++) {
       pred(i,j) = exp(log_pred(i,j));
     }
   }
   
  
  // report section
  ADREPORT(log_pred);
  
  // if(n_strata_biomass > 1) {
  //   vector<Type> tot_biomass_pred;
  //   tot_biomass_pred = biomass_pred.rowwise().sum();
  //   vector<Type> log_tot_biomass_pred;
  //   log_tot_biomass_pred = log(tot_biomass_pred);
  //   ADREPORT(log_tot_biomass_pred);
  // }
  // 

  REPORT(log_pred);
  REPORT(pred);
  REPORT(jnll);
  
  // jnll = dummy * dummy;        // Uncomment when debugging code
  nll = jnll.sum();
  return nll;
  
}
