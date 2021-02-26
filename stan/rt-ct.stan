functions {
#include functions/gaussian_process.stan
#include functions/ct.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // number of ct samples
  int t; // time of considered
  int tt[N]; // time each sample taken
  vector[N] ct; // count with ct value 
  real init_inf_prob; // initial probability of infection
  int ctmax; // maximum number of days post infection considered for ct values
  vector[ctmax] ct_inf_mean; // mean CT by day since infection
  vector[ctmax] ct_inf_sd; // standard deviation of CT by day of infection
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M;
  real L;
}

transformed data {
  matrix[t - 1, M] PHI = setup_gp(M, L, t - 1);  
  real intercept = logit(init_inf_prob);
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
}

transformed parameters {
  vector[t] growth;
  vector[t] prob_inf;
  
  // Infections from growth
  growth[1] = 0;
  growth[2:t] = update_gp(PHI, M, L, alpha, rho, eta, 0);
  prob_inf = inv_logit(intercept + cumulative_sum(growth));
}

model {
  vector[ctmax] lrit[t - 1];

  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  
  lrit = rel_inf_prob(prob_inf, ctmax, t);
  
  for (n in 1:N) {
    vector[ctmax] lps = lrit[tt[n] - 1];
    for (k in 1:ctmax) {
      lps[k] += normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[k]);
    }
    target += log_sum_exp(lps);
  }
}

