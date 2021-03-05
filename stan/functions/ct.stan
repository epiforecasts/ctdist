// Calculate log normal cumulative density
// This is phi_a(X) which is the probability
// a ct value of X is detected a days after infection
vector ct_threshold_prob(real dt, vector ct_inf_mean, vector ct_inf_sd) {
  int ctmax = num_elements(ct_inf_mean);
  vector[ctmax] ldtp;
  for (k in 1:ctmax) {
    ldtp[k] = normal_lcdf(dt | ct_inf_mean[k], ct_inf_sd[k]);
  }
  return(ldtp);
}
vector rel_threshold_prob(vector[] lrit, int t, vector ldtp) {
  vector[t] ldtpt;
  for (i in 1:t) {
    // This log sums lrit (pi) and ldtp (phi) for
    // the denominator
    ldtpt[i] = log_sum_exp(lrit[i] + reverse(ldtp));
  }
  return(ldtpt);
}
// Calculate log normal density for each observation and day since infection
// This is p_a(X)  = P(ct = X when infected a days ago)
// p_a(X) = P(ct = X | detected) * phi(X)
vector[] ct_log_dens(real[] ct, vector ct_inf_mean, vector ct_inf_sd, vector ldtp) {
  int N = num_elements(ct);
  int ctmax = num_elements(ct_inf_mean);
  vector[ctmax] ctlgd[N];
  for (n in 1:N) {
    for (k in 1:ctmax) {
      // This now accounts for p_a(X) = P(ct = X | detected) * phi(X)
      ctlgd[n][k] = normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[k]) + ldtp[k];
    }
  }
  return(ctlgd);
}
// relative infection probability for each reference time look up
vector[] rel_inf_prob(vector prob_inf, vector ldtp, int ctmax, int t) {
  int p;
  vector[ctmax] lrit[t - ctmax];
  for (i in (ctmax + 1):t) {
    p = i - ctmax;
    lrit[p] = rep_vector(1e-8, ctmax);
    for (j in 1:ctmax) {
      lrit[p][j] += prob_inf[i - j];
    } 
   lrit[p] = log(lrit[p]);
  }
  return(lrit);
}
// log likelihood across CT observations
real ct_loglik(real[] ct, int start, int end, int[] tt, vector[] lrit,
               vector[] ctlgd, vector ldtpt) {
  real tar = 0;
  int t;
  for (n in start:end) {
    t = tt[n];
    // log_exp_sum on top - log denominator 
    tar += log_sum_exp(lrit[t] + ctlgd[n]) - ldtpt[t];
  }
  return(tar);
}
