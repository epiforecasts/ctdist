# Packages -----------------------------------------------------------------
library("rstan")
library("EpiNow2")
library("here")
library("dplyr")

# Set up parallel ---------------------------------------------------------
options(mc.cores = 4)

# Load data ---------------------------------------------------------------
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))
ep_vacc <- readRDS(here("data", "ct_summarised.rds"))
vacc_grp <- readRDS(here("data", "national_vacc_coverage.rds"))

# Data for stan -----------------------------------------------------------
# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 5000)
ep_raw_vacc <- ep_raw_vacc[samples, ]

# define stan list
dat <- list()
dat$N <- nrow(ep_raw_vacc)
dat$t <- max(ep_raw_vacc$time) + 1
dat$tt <- ep_raw_vacc$time + 1
dat$ct <- ep_raw_vacc$p2ch1cq


# gaussian process data
dat$M <- ceiling(dat$t / 3)
dat$L <- 2
lsp <- tune_inv_gamma(floor(dat$t/3), dat$t)
dat$lengthscale_alpha <- lsp$alpha
dat$lengthscale_beta <- lsp$beta

# Load model --------------------------------------------------------------
mod <- stan_model(here::here("stan", "rt-ct.stan"))

# Fit model ---------------------------------------------------------------
res <- sampling(mod, data = dat)
