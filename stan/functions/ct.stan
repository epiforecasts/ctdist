// Calculate relative infection probability for each reference time
vector[] rel_inf_prob(vector inf, int ctmax, int t) {
  int p;
  int s;
  vector[ctmax] lrit[t - 1];
  for (i in 2:t) {
    p = i - 1;
    s = min(p, ctmax);
    lrit[p] = rep_vector(1e-8, ctmax);
    for (j in 1:s) {
      lrit[p][j] += inf[i - j];
    } 
    lrit[p] = log(lrit[p] / sum(lrit[p]));
  }
  return(lrit);
}

real ct_mixture(real[] ct, int start, int end, int[] tt, vector[] lrit,
                vector ct_inf_mean, vector ct_inf_sd, int ctmax) {
  real tar = 0;
    for (n in start:end) {
    vector[ctmax] lps = lrit[tt[n] - 1];
    for (k in 1:ctmax) {
      lps[k] += normal_lpdf(ct[n - start + 1] | ct_inf_mean[k], ct_inf_sd[k]);
    }
    tar += log_sum_exp(lps);
  }
  return(tar);
}
