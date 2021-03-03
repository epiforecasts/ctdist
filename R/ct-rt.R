# Packages -----------------------------------------------------------------
library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(EpiNow2)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
color_scheme_set("brightblue")


# Source utilities --------------------------------------------------------
source(here("R", "ct-rt-utils.R"))

# Set up parallel ---------------------------------------------------------
# number of cores (only uses up to 4 chains)
cores <- 4
# number of threads per core (so available cores / cores)
threads <- 4

# Set up CmdStan if required ----------------------------------------------
# install_cmdstan(cores = cores)

# Load data ---------------------------------------------------------------
# must contain, viral load or proxy, time, date
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))

# Data for stan -----------------------------------------------------------
# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 1000)
ep_raw_vacc <- ep_raw_vacc[samples, ]
min_date <- min(ep_raw_vacc$date_specimen, na.rm = TRUE)

# define CT post infection (loosely inspired Hay et al)
ct <- tibble(mean = c(40 - 0:4*5, 18 + 0:10*2),
             sd = c(rep(1, 5), rep(2, 11)))

# define observations
dat <- stan_data(ep_raw_vacc, 
                 load_vec = "p2ch1cq",
                 overall_prob = 1,
                 ct =  ct,
                 dt = 30,
                 gt = get_generation_time(
                   disease = "SARS-CoV-2", source = "ganyani", max = 15
                   ), gp_m = 0.1, gp_ls = c(7, NA)
                 )

# Load model --------------------------------------------------------------
mod <- cmdstan_model(here("stan", "rt-ct.stan"), include_paths = "stan",
                     cpp_options = list(stan_threads = TRUE))

# Fit model ---------------------------------------------------------------
fit <- mod$sample(data = dat, parallel_chains = cores, 
                  threads_per_chain = threads)
# check
fit$cmdstan_diagnose()

# summarise fit
fit$summary()

# Plot variables over time ------------------------------------------------
plot_trend(fit, "prob_inf", date_start = min_date - dat$ctmax) +
  labs(y = "Relative probability of infection", x = "Date")

ggsave(here("figures", "ct-relative-prob-inf.pdf"), width = 7, height = 5)

plot_trend(fit, "r", date_start = min_date - dat$ctmax - 1) +
  labs(y = "Daily growth rate", x = "Date") +
  geom_hline(yintercept = 0, linetype = 2)

ggsave(here("figures", "ct-growth.pdf"), width = 7, height = 5)

plot_trend(fit, "R",  date_start = min_date - dat$ctmax + 7) +
  labs(y = "Effective reproduction number", x = "Date") +
  geom_hline(yintercept = 1, linetype = 2)

ggsave(here("figures", "ct-Rt.pdf"), width = 7, height = 5)
