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
linear_fit <- brm(p2ch1cq ~
             age_group + sgene_result + s(time, k = 5) +
             First * sgene_result + Second * sgene_result,
           data = ep_raw_vacc,
           family = gaussian(link = "identity"))

saveRDS(linear_fit, "output/ct_fit_linear_brms.rds")

# log fit
log_fit <- brm(p2ch1cq ~
             age_group + sgene_result + s(time, k = 5) +
             First * sgene_result + Second * sgene_result,
           data = ep_raw_vacc,
           family = gaussian(link = "log"))

saveRDS(log_fit, "output/ct_fit_log_brms.rds")