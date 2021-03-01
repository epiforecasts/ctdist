# Packages -----------------------------------------------------------------
library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(EpiNow2)
library(dplyr)
library(tidyr)
library(ggplot2)
color_scheme_set("brightblue")


# Source utilities --------------------------------------------------------
source(here("R", "ct-rt-utils.R"))

# Set up parallel ---------------------------------------------------------
cores <- 4
threads <- 4

# Set up CmdStan if required ----------------------------------------------
# install_cmdstan(cores = cores)

# Load data ---------------------------------------------------------------
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))

# Data for stan -----------------------------------------------------------
# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 10)
ep_raw_vacc <- ep_raw_vacc[samples, ]

# define observations
dat <- stan_data(ep_raw_vacc, 
                 init_prob = 0.1, 
                 ct_mean =  c(40 - 4*0:6, 18 + (0:4)*3), 
                 ct_sd =  rep(1, 12),
                 dt = 30,
                 gt = get_generation_time(
                   disease = "SARS-CoV-2", source = "ganyani", max = 15
                   ), gp_m = 0.3, gp_ls = c(7, NA)
                 )

# Load model --------------------------------------------------------------
mod <- cmdstan_model(here("stan", "rt-ct.stan"), include_paths = "stan",
                     cpp_options = list(stan_threads = TRUE))

# Fit model ---------------------------------------------------------------
fit <- mod$sample(data = dat, parallel_chains = cores, 
                  threads_per_chain = threads)

# summarise fit
fit$summary()

# check
fit$cmdstan_diagnose()

# Plot variables over time ------------------------------------------------
plot_trend(fit, "prob_inf") +
  labs(y = "Relative probability of infection", x = "Date")

ggsave(here("figures", "ct-relative-prob-inf.pdf"), width = 7, height = 5)

plot_trend(fit, "growth") +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)

ggsave(here("figures", "ct-growth.pdf"), width = 7, height = 5)

plot_trend(fit, "R", min(ep_raw_vacc$date_specimen) - dat$ctmax + 7) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)

ggsave(here("figures", "ct-Rt.pdf"), width = 7, height = 5)
