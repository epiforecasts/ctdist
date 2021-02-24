# Packages ----------------------------------------------------------------
library(brms)

# Parallel ----------------------------------------------------------------
options(mc.cores = 4)

# Load data ---------------------------------------------------------------
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))

# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 5000)
ep_raw_vacc <- ep_raw_vacc[samples, ]

# Fit ---------------------------------------------------------------------
# linear fit
fit <- brm(p2ch1cq ~
             age_group + sgene_result + s(time, k = 10) +
             First * sgene_result + Second * sgene_result,
           data = ep_raw_vacc,
           family = gaussian(link = "identity"))

saveRDS(fit, "output/ct_fit_linear.rds")

# log fit
fit <- brm(p2ch1cq ~
             age_group + sgene_result + s(time, k = 10) +
             First * sgene_result + Second * sgene_result,
           data = ep_raw_vacc,
           family = gaussian(link = "log"))

saveRDS(fit, "output/ct_fit_logrds")