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
  vector[ct_max] rel_inf_prob[t - 1];
  
  // Infections from growth
  growth = update_gp(PHI, M, L, alpha, rho, eta, 0);
  infections[1] = i0;
  infections[2:t] = i0 * exp(cumulative_sum(growth));
  
  for (i in 2:t) {
    int s = min(1, i - ct_max);
    rel_inf_prob[i] = rep_vector(0, ct_max);
      for ()
    rel_inf_prob[i] = infections[s:i] ./ sum(infections[s:i]);
  }
}

model {
  i0 ~ exp(1/1000)
  sigma ~ normal(0, 1) T[0,];
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  
  for (n in 1:N) {
    vector[K] lps = log(rel_inf_prob[[tt[n]]]);
    for (k in 1:ct_max) {
      lps[k] += normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[i]);
    }
    target += log_sum_exp(lps);
  }
}

