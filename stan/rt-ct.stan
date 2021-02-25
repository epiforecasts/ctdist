functions {
#include functions/gaussian_process.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // Number of ct samples
  int t;
  int <lower = 1> M;
  real L;
  int tt[N]; // time of each sample
  vector ct[N]; // count with ct value 
  int i0;
  int ctmax;
  vector ct_inf_mean[ctmax];
  vector ct_inf_sd[ctmax];
  real lengthscale_alpha;            // alpha for gp lengthscale prior
  real lengthscale_beta;             // beta for gp lengthscale prior
}

transformed data {
  matrix[t - 1, M] PHI = setup_gp(M, L, t - 1);  
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  real <lower = 0> sigma;
}

transformed parameters {
  vector[t-1] growth;
  vector[t] infections;
  vector[ctmax] rel_inf_prob[t - 1];
  
  // Infections from growth
  growth = update_gp(PHI, M, L, alpha, rho, eta, 0);
  infections[1] = i0;
  infections[2:t] = i0 * exp(cumulative_sum(growth));
  
  for (i in 2:t) {
    int s = min(1, i - ct_max);
    rel_inf_prob[i] = rep_vector(0, ctmax);
      for ()
    rel_inf_prob[i] = infections[s:i] ./ sum(infections[s:i]);
  }
}

model {
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  
  for (n in 1:N) {
    vector[ctmax] lps = log(rel_inf_prob[[tt[n]]]);
    for (k in 1:ctmax) {
      lps[k] += normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[k]);
    }
    target += log_sum_exp(lps);
  }
}

