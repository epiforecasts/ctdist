# ctdist

> **Archived.** Exploratory work from February–March 2021. It is no longer maintained, and the code will not run as it stands because the data it reads are not in the repository.

Exploratory models for inferring changes in the distribution of SARS-CoV-2 cycle threshold (Ct) values over time, and for linking Ct values to the time-varying reproduction number.

- `R/direct-*.R` model Ct values directly over time and by covariates such as age group and vaccination, with Stan (`stan/direct-ct.stan`), brms and mgcv.
- `R/ct-rt.R` and `stan/rt-ct.stan` combine a Gaussian process for Rt with Ct value trajectories after infection, building on EpiNow2.
- `R/process-data.R` prepares the input data, which were not shared publicly.

Later work on estimating infections from prevalence and Ct data continued in [epiforecasts/inc2prev](https://github.com/epiforecasts/inc2prev).
