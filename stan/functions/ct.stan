// Calculate log normal cumulative density
vector ct_threshold_prob(real dt, vector ct_inf_mean, vector ct_inf_sd) {
  int ctmax = num_elements(ct_inf_mean);
  vector[ctmax] ldtp;
  for (k in 1:ctmax) {
    ldtp[k] = normal_lcdf(dt | ct_inf_mean[k], ct_inf_sd[k]);
  }
  return(ldtp);
}
vector rel_threshold_prob(vector ldtp, vector[] lrit, int t, int ctmax) {
  vector[t] ldtpt;
  for (i in 1:t) {
    ldtpt[i] = log_sum_exp(lrit[i] + ldtp);
  }
  return(ldtpt);
}
// Calculate log normal density for each observation and day since infection
vector[] ct_log_dens(real[] ct, vector ct_inf_mean, vector ct_inf_sd) {
  int N = num_elements(ct);
  int ctmax = num_elements(ct_inf_mean);
  vector[ctmax] ctlgd[N];
  for (n in 1:N) {
    for (k in 1:ctmax) {
      ctlgd[n][k] = normal_lpdf(ct[n] | ct_inf_mean[k], ct_inf_sd[k]);
    }
  }
  return(ctlgd);
}
// relative infection probability for each reference time look up
vector[] rel_inf_prob(vector inf, int ctmax, int t) {
  int p;
  vector[ctmax] lrit[t - ctmax];
  for (i in (ctmax + 1):t) {
    p = i - ctmax;
    lrit[p] = rep_vector(1e-8, ctmax);
    for (j in 1:ctmax) {
      lrit[p][j] += inf[i - j];
    } 
   lrit[p] = log(lrit[p]);
  }
  return(lrit);
}
// create mixture of days since infection
real ct_loglik(real[] ct, int start, int end, int[] tt, vector[] lrit,
vector[] ctlgd, vector ldtpt, int ctmax) {
  real tar = 0;
  int t;
  for (n in start:end) {
    t = tt[n];
    tar += log_sum_exp(lrit[t] + ctlgd[n]) - ldtpt[t];
  }
  return(tar);
}
