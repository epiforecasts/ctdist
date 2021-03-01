functions {
#include functions/gaussian_process.stan
#include functions/ct.stan
#include functions/rt.stan
#include functions/generated_quantities.stan
}

data {
  int N; // number of ct samples
  int t; // time considered
  int ut; // time considered + ctmax
  int tt[N]; // time each sample taken
  real ct[N]; // count with ct value 
  real dt; // detection threshold
  real overall_prob; // overall probability of infection
  int ctmax; // maximum number of days post infection considered for ct values
  vector[ctmax] ct_inf_mean; // mean CT by day since infection
  vector[ctmax] ct_inf_sd; // standard deviation of CT by day of infection
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real gtm[2]; // mean and standard deviation (sd) of the mean generation time
  real gtsd[2]; // mean and sd of the sd of the generation time
  int gtmax; // maximum number of days to consider for the generation time
}

transformed data {
  // set up approximate gaussian process
  matrix[ut, M] PHI = setup_gp(M, L, ut);  
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
  vector[ut] gp;
  vector[ut] prob_inf;
  // update gaussian process
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // relative probability of infection
  prob_inf = inv_logit(gp);
  prob_inf = prob_inf / sum(prob_inf);
  prob_inf = overall_prob * prob_inf;
}

model {
  vector[ctmax] lrit[t];
  vector[t] ldtpt;
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  // relative log probability of infection for each t
  lrit = rel_inf_prob(prob_inf, ctmax, ut);
  // log prob of detection for each t
  ldtpt = rel_threshold_prob(ldtp, lrit, t, ctmax);
  // update likelihood (in parallel)
  target += reduce_sum(ct_loglik, ct, 1, tt, lrit, ctlgd, ldtpt, ctmax);
  // if using rstan/no reduce_sum comment out L66 and use L68 instead
  //target += ct_loglik(ct, 1, N, tt, lrit, ctlgd, ldtpt, ctmax);
}

generated quantities {
  vector[ut - 7] R;
  vector[ut - 1] r;
  // sample generation time
  real gtm_sample = normal_rng(gtm[1], gtm[2]);
  real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // calculate Rt using infections and generation time
  R = calculate_Rt(prob_inf, 7, gtm_sample, gtsd_sample, gtmax, 1);
  // calculate growth
  r = calculate_growth(prob_inf, 1);
}
