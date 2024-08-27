write("data {
  int<lower=1> J; //total individuals
  int<lower=1> N_z; //total number of phoontype meausres (Z)
  array[N_z] int<lower=1> ind_z; //index of individual measurements (z)
  vector[N_z] x; //environmental covariate
  vector[N_z] z; //phenotype
}

parameters {
  //fixed population effects
  real mu_0; //z population intercept
  real beta_x; //z population slope
  real sigma_0; //z population dispersion
  
  //random effects for z
  vector<lower=0>[3] sd_RN; //RN parameter sds
  matrix[J,3] std_dev_RN; //individual-level RN deviations
  cholesky_factor_corr[3] R_chol; //RN parameter correlations
}

transformed parameters {
  matrix[J,3] RNj = std_dev_RN * diag_pre_multiply(sd_RN, R_chol)' ;
}

model{
  //separate RN parameters
  vector[J] mu_0j = col(RNj,1); //intercepts
  vector[J] beta_xj = col(RNj,2); //slopes
  vector[J] sigma_0j = col(RNj,3); //residuals
  
  //initialize vectors for response models
  vector[N_z] mu; //linear predictor of phenotype expectation
  vector[N_z] sigma; //linear predictor of phenotype dispersion
  
  //RN model
  mu = mu_0 + mu_0j[ind_z] + (beta_x + beta_xj[ind_z]) .* x;
  sigma = sqrt(exp(sigma_0 + sigma_0j[ind_z]));
  z ~ normal(mu, sigma);
  
  //model priors
  
  //fixed effects
  mu_0 ~ normal(0,1);
  beta_x ~ normal(0,1);
  sigma_0 ~ normal(0,1);
  
  //random effects
  sd_RN ~ exponential(2);
  R_chol ~ lkj_corr_cholesky(2);
  to_vector(std_dev_RN) ~ std_normal();
}

generated quantities{


matrix[3,3] R = R_chol * R_chol';        //RN correlation matrix
matrix[3,3] S = diag_matrix(sd_RN);      //RN SD matrix
matrix[3,3] P = S*R*S;                   //RN covariance matrix


vector<lower=0>[3] V_P = sd_RN .* sd_RN; //RN variances


}"
, file="Mod4.stan")

library(rstan)
  Mod4 <- stan("Mod4.stan", data = data, 
               chains = 1, iter = 5000,  warmup = 1000, 
               thin = 1, cores = 1)
  dat <- as.data.frame(round(summary(Mod4)$summary[,c(1,4,6,8,9,10)],3))
  taildat <- tail(dat)
}
