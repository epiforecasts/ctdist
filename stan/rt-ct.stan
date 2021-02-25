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
  matrix[t, M] PHI = setup_gp(M, L, t);  
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  real <lower = 0> sigma;
}

transformed parameters {
  vector[t] growth;
  vector[t] infections;
  vector[ctmax] lrit[t - 1];
  
  // Infections from growth
  growth = update_gp(PHI, M, L, alpha, rho, eta, 0);
  prob_inf = i0 * inv_logit(cumulative_sum(growth));
  prob_inf = prob_inf / sum(prob_inf);
  
  for (i in 2:t) {
    int s = min(i, ctmax);
    lrit[i] = rep_vector(1e-8, ctmax);
    for (j in 1:s) {
      lrit[i][j] = prob_inf[i - s + 1];
    } 
    lrit[i] = log(lrit[i] / sum(lrit[i]));
  }
}

model {
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  
  for (n in 1:N) {
    vector[ctmax] lps = lrit[[tt[n]]];
    for (k in 1:ctmax) {
      lps[k] += normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[k]);
    }
    target += log_sum_exp(lps);
  }
}

