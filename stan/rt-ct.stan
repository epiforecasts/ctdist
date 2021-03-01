functions {
#include functions/gaussian_process.stan
#include functions/ct.stan
#include functions/rt.stan
}

data {
  int N; // number of ct samples
  int t; // time considered
  int ut; // time considered + ctmax
  int tt[N]; // time each sample taken
  real ct[N]; // count with ct value 
  real dt; // detection threshold
  real init_inf_prob; // initial probability of infection
  int ctmax; // maximum number of days post infection considered for ct values
  vector[ctmax] ct_inf_mean; // mean CT by day since infection
  vector[ctmax] ct_inf_sd; // standard deviation of CT by day of infection
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M;
  real L;
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
}

transformed data {
  // set up gaussian process
  matrix[ut, M] PHI = setup_gp(M, L, ut);  
  real intercept = logit(init_inf_prob);
  // calculate log density for each observed ct and day since infection
  vector[ctmax] ctlgd[N] = ct_log_dens(ct, ct_inf_mean, ct_inf_sd);
  // calculate log of probability CT below threshold
  vector[ctmax] ldtp = ct_threshold_prob(dt, ct_inf_mean, ct_inf_sd);
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
}

transformed parameters {
  vector[ut] growth;
  vector[ut] prob_inf;
  // relative probability of infection from growth
  growth = update_gp(PHI, M, L, alpha, rho, eta, 0);
  prob_inf = inv_logit(intercept + cumulative_sum(growth));
  prob_inf = prob_inf / sum(prob_inf);
}

model {
  vector[ctmax] lrit[t];
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  // calculate relative logged probability of infection for each t
  lrit = rel_inf_prob(prob_inf, ctmax, ut);
  // update likelihood (in parallel)
  target += reduce_sum(ct_mixture, ct, 1, tt, lrit, ctlgd, ldtp, ctmax);
}

generated quantities {
  vector[ut - 7] R;
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(prob_inf, 7, gtm_sample, gtsd_sample, gtmax, 1);
}
