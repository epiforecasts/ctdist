functions {
#include functions/gaussian_process.stan
#include functions/ct.stan
#include functions/rt.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // number of ct samples
  int t; // time considered
  int tt[N]; // time each sample taken
  real ct[N]; // count with ct value 
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
  int ut = t + ctmax;
  matrix[ut - 1, M] PHI = setup_gp(M, L, ut - 1);  
  real intercept = logit(init_inf_prob);
  vector[ctmax] ctlgd[N] = ct_log_dens(ct, ct_inf_mean, ct_inf_sd);
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
}

transformed parameters {
  vector[ut] growth;
  vector[ut] prob_inf;
  
  // Infections from growth
  growth[1] = 0;
  growth[2:ut] = update_gp(PHI, M, L, alpha, rho, eta, 0);
  prob_inf = inv_logit(intercept + cumulative_sum(growth));
}

model {
  vector[ctmax] lrit[t];

  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  
  lrit = rel_inf_prob(prob_inf, ctmax, ut);
  target += reduce_sum(ct_mixture, ct, 1, tt, lrit, ctlgd, ctmax);
}

generated quantities {
  vector[t-7] R;
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(prob_inf, 7, gtm_sample, gtsd_sample, gtmax, 1);
}
