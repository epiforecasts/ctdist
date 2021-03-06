functions {
#include functions/gaussian_process.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // Number of ct samples
  int t;
  int AG;
  int <lower = 1> M;
  real L;
  int tt[N]; // time of each sample
  int agegrp[N]; // age group of each sample
  vector[N] ct; // ct value of each sample
  matrix[AG, t] vacc_cov; // vaccine coverage of each age group at each time point
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
  vector[AG] beta0;
  vector[AG] beta1;
  vector[AG] beta2;
}

transformed parameters {
  vector[t] gp;
  vector[N] mu;
  
  // Linear model using the spectral densities
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  
  for(i in 1:N){
    mu[i] = beta0[agegrp[i]] + beta1[agegrp[i]] * gp[tt[i]] + beta2[agegrp[i]] * vacc_cov[agegrp[i], tt[i]];
  }
}

model {
  ct ~ normal(mu, sigma);
  sigma ~ normal(0, 1) T[0,];
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, 1);
  eta ~ std_normal();
  beta0 ~ std_normal();
  beta1 ~ std_normal();
  beta2 ~ std_normal();
}

