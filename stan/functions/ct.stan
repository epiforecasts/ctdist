// Calculate relative infection probability for each reference time
vector[] rel_inf_prob(vector prob_inf, int ctmax, int t) {
  int p;
  int s;
  vector[ctmax] lrit[t - 1];
  for (i in 2:t) {
    p = i - 1;
    s = min(p, ctmax);
    lrit[p] = rep_vector(1e-8, ctmax);
    for (j in 1:s) {
      lrit[p][j] += prob_inf[i - j];
    } 
    lrit[p] = log(lrit[p] / sum(lrit[p]));
  }
  return(lrit);
}
