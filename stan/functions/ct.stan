// Calculate relative infection probability for each reference time
vector[] rel_inf_prob(vector inf, int ctmax, int t) {
  int p;
  vector[ctmax] lrit[t - ctmax];
  for (i in (ctmax + 1):t) {
    p = i - ctmax;
    lrit[p] = rep_vector(1e-8, ctmax);
    for (j in 1:ctmax) {
      lrit[p][j] += inf[i - j];
    } 
    lrit[p] = log(lrit[p] / sum(lrit[p]));
  }
  return(lrit);
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
// create mixture of days since infection
real ct_mixture(real[] ct, int start, int end, int[] tt, vector[] lrit,
                vector[] ctlgd, int ctmax) {
  real tar = 0;
    for (n in start:end) {
    vector[ctmax] lps = lrit[tt[n] - 1];
    for (k in 1:ctmax) {
      lps[k] += ctlgd[n][k];
    }
    tar += log_sum_exp(lps);
  }
  return(tar);
}
