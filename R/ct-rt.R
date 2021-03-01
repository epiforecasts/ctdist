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

# Set up parallel ---------------------------------------------------------
cores <- 4
threads <- 4

# Set up CmdStan if required ----------------------------------------------
# install_cmdstan(cores = cores)

# Load data ---------------------------------------------------------------
ep_raw_vacc <- readRDS(here("data", "ct_covariates.rds"))

# Data for stan -----------------------------------------------------------
# subsample available data
samples <- sample(1:nrow(ep_raw_vacc), 1000)
ep_raw_vacc <- ep_raw_vacc[samples, ]

# define observations
dat <- list()
dat$N <- nrow(ep_raw_vacc)
dat$t <- max(ep_raw_vacc$time) + 1
dat$tt <- ep_raw_vacc$time + 1
dat$ct <- ep_raw_vacc$p2ch1cq

# define initial probability of infection (10%)
# I am not sure the absolute number is meaningful
dat$init_inf_prob <- 0.01

# define ct parameters
# assume ct lower than 30 threshold for 16 days 
# initial linear decrease followed by linear increase
dat$ctmax <- 12
dat$ct_inf_mean <- c(50 - 4*0:5, 20 + (0:5)*2)
dat$ct_inf_sd <- rep(0.5, 12)

# gaussian process parameters
dat$M <- ceiling(dat$t / 3)
dat$L <- 2
lsp <- tune_inv_gamma(7, dat$t)
dat$lengthscale_alpha <- lsp$alpha
dat$lengthscale_beta <- lsp$beta

# add generation time assumption
gt <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
dat$gtm <- unlist(gt[c("mean", "mean_sd")])
dat$gtsd <- unlist(gt[c("sd", "sd_sd")])
dat$gtmax <- 15

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
plot_trend <- function(fit, var, max_date = max(ep_raw_vacc$date_specimen)) {
  fit$summary(variables = var, 
              ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))) %>% 
    mutate(time = 1:n()) %>% 
    ggplot() +
    aes(x = time, y = `50%`, ymin = `5%`, ymax = `95%`) + 
    geom_line(col = "lightblue", size = 1.4) +
    geom_ribbon(fill = "lightblue", alpha = 0.4,
                col = "lightblue", size = 0.8) +
    geom_ribbon(fill = "lightblue", alpha = 0.4,
                col = "lightblue", size = 0.8,
                aes(ymin = `20%`, ymax = `80%`)) +
    theme_minimal()
}

plot_trend(fit, "prob_inf") +
  labs(y = "Probability of infection", x = "Time")

plot_trend(fit, "R") +
  labs(y = "Effective reproduction number", x = "Time") +
  geom_hline(yintercept = 1, linetype = 2)
