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
dat$AG <- length(unique(ep_raw_vacc$age_group))
dat$time <- 0:(dat$t - 1)
dat$tt <- ep_raw_vacc$time
dat$agegrp <- as.numeric(ep_raw_vacc$age_group)
dat$ct <- ep_raw_vacc$p2ch1cq
dat$M <- ceiling(dat$t / 3)
dat$L <- dat$t * 2

# tune lengthscale of the gamma process (Joel to replace)
lsp <- tune_inv_gamma(floor(dat$t/3), dat$t)
dat$lengthscale_alpha <- lsp$alpha
dat$lengthscale_beta <- lsp$beta


# add in vaccine coverage by age group
vacc_grp <- vacc_grp %>%
  select(vaccination_date, age_group, cum_prop) %>%
  group_by(age_group) %>%
  mutate(time = 0:(n() - 1)) %>%
  filter(time %in% dat$time)

dat$vacc_cov <- matrix(data = vacc_grp$cum_prop, nrow = dat$AG, ncol = dat$t)


# Load model --------------------------------------------------------------
mod <- stan_model(here::here("stan", "direct-ct.stan"))


# Fit model ---------------------------------------------------------------
res <- sampling(mod, data = dat)
