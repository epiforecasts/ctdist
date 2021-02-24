functions {
  // exponential quadratic kernal
  real spd_SE(real alpha, real rho, real w) {
    real S;
    S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
    return S;
  }

  // basis function for approximate hilbert space gp
  // see here for details: https://arxiv.org/pdf/2004.11408.pdf
  vector phi_SE(real L, int m, vector x) {
    vector[rows(x)] fi;
    fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
    return fi;
  }

  // eigenvalues for approximate hilbert space gp
  // see here for details: https://arxiv.org/pdf/2004.11408.pdf
  real lambda(real L, int m) {
    real lam;
    lam = ((m*pi())/(2*L))^2;
    return lam;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int N; // Number of ct samples
  int t;
  int AG;
  int <lower = 1> M;
  real L;
  vector[t] time; // vector of 1:t
  int tt[N]; // time of each sample
  int agegrp[N]; // age group of each sample
  vector[N] ct; // ct value of each sample
  matrix[AG, t] vacc_cov; // vaccine coverage of each age group at each time point
  real lengthscale_alpha;            // alpha for gp lengthscale prior
  real lengthscale_beta;             // beta for gp lengthscale prior
}

transformed data {
  matrix[t, M] PHI;

  for(m in 1:M) {
    PHI[, m] = phi_SE(L, m, time);
  }
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
  vector[M] diagSPD;
  vector[M] SPD_beta;
  vector[N] mu;
  
  for(m in 1:M) {
    // Spectral density calculation
    diagSPD[m] =  sqrt(spd_SE(alpha, rho, sqrt(lambda(L, m))));
  }

  // Linear model using the spectral densities
  SPD_beta = diagSPD .* eta;
  
  
  for(i in 1:N){
    mu[i] = beta0[agegrp[i]] + beta1[agegrp[i]] * SPD_beta[tt[i]] + beta2[agegrp[i]] * vacc_cov[agegrp[i], tt[i]];
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

