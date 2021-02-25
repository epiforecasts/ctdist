functions {
#include functions/gaussian_process.stan
#include functions/ct.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // Number of ct samples
  int t;
  int <lower = 1> M;
  real L;
  int tt[N]; // time of each sample
  vector[N] ct; // count with ct value 
  int ctmax;
  vector[ctmax] ct_inf_mean;
  vector[ctmax] ct_inf_sd;
  real lengthscale_alpha;            // alpha for gp lengthscale prior
  real lengthscale_beta;             // beta for gp lengthscale prior
}

transformed data {
  matrix[t, M] PHI = setup_gp(M, L, t);  
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
  growth = update_gp(PHI, M, L, alpha, rho, eta, 0);
  prob_inf = exp(cumulative_sum(growth));
  //prob_inf = prob_inf / sum(prob_inf);
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

