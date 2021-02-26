# Packages -----------------------------------------------------------------
library("rstan")
library("EpiNow2")
library("here")
library("dplyr")

# Set up parallel ---------------------------------------------------------
options(mc.cores = 4)

# Load data ---------------------------------------------------------------
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))

# Data for stan -----------------------------------------------------------
# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 500)
ep_raw_vacc <- ep_raw_vacc[samples, ]

# define observations
dat <- list()
dat$N <- nrow(ep_raw_vacc)
dat$t <- max(ep_raw_vacc$time) + 2
dat$tt <- ep_raw_vacc$time + 2
dat$ct <- ep_raw_vacc$p2ch1cq

# define initial probability of infection (1%)
dat$init_inf_prob <- 0.001
# define ct parameters
# assume ct lower than 30 threshold for 16 days 
# initial linear decrease followed by linear increase
dat$ctmax <- 16
dat$ct_inf_mean <- c(40 - 3*0:5, 20 + 1:10)
dat$ct_inf_sd <- rep(1, 16)

# gaussian process parameters
dat$M <- ceiling(dat$t / 4)
dat$L <- 2
lsp <- tune_inv_gamma(7, dat$t)
dat$lengthscale_alpha <- lsp$alpha
dat$lengthscale_beta <- lsp$beta

# Load model --------------------------------------------------------------
mod <- stan_model(here::here("stan", "rt-ct.stan"))

# Fit model ---------------------------------------------------------------
res <- sampling(mod, data = dat)
