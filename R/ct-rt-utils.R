library(EpiNow2)
library(data.table)

# define required stan data
stan_data <- function(obs, load_vec = "p2ch1cq", overall_prob = 1, 
                      ct, dt = 30,
                      gt = get_generation_time(
                        disease = "SARS-CoV-2", source = "ganyani", max = 15
                        ), gp_m = 0.3, gp_ls = c(7, NA)
                      ) {
  
  obs <- as.data.table(obs)

  # define observations
  dat <- list()
  dat$n <- nrow(obs)
  dat$t <- max(obs$time) + 1
  dat$tt <- obs$time + 1
  dat$ct <- obs[["p2ch1cq"]]
  dat$dt <- dt
  
  # observations by time point
  dat$nt <- copy(obs)[, .(nobs = .N), by = time][order(time)]
  dat$nt <- dat$nt[, time := time + 1]
  dat$nt <- merge(data.table(time = 1:dat$t), dat$nt, by = "time", all = TRUE)
  dat$nt <- dat$nt[is.na(nobs), nobs := 0]
  dat$nt <- dat$nt$nobs

  # overall probability scaling factor
  dat$overall_prob <- overall_prob
  
  # define ct parameters + unobserved time
  dat$ctmax <- length(ct$mean)
  dat$ct_inf_mean <- ct$mean
  dat$ct_inf_sd <- ct$sd
  dat$ut <- dat$t + dat$ctmax
  
  # gaussian process parameters
  dat$M <- ceiling(dat$ut * gp_m)
  dat$L <- 2
  if (is.na(gp_ls[2])) {
    gp_ls[2] <- dat$ut
  }
  lsp <- tune_inv_gamma(gp_ls[1], gp_ls[2])
  dat$lengthscale_alpha <- lsp$alpha
  dat$lengthscale_beta <- lsp$beta
  
  #define generation time
  dat$gtm <- unlist(gt[c("mean", "mean_sd")])
  dat$gtsd <- unlist(gt[c("sd", "sd_sd")])
  dat$gtmax <- unlist(gt[c("max")])
  return(dat)
}
# plot packages
library(cmdstanr)
library(dplyr)
library(ggplot2)

# plot trend with date over time
plot_trend <- function(fit, var, date_start) {
  fit$summary(variables = var, 
              ~ quantile(.x, probs = c(0.05, 0.2, 0.5, 0.8, 0.95))) %>% 
    mutate(time = 1:n(),
           date = date_start + time - 1) %>% 
    ggplot() +
    aes(x = date, y = `50%`, ymin = `5%`, ymax = `95%`) + 
    geom_line(col = "lightblue", size = 1.4) +
    geom_ribbon(fill = "lightblue", alpha = 0.4,
                col = "lightblue", size = 0.6) +
    geom_ribbon(fill = "lightblue", alpha = 0.4,
                col = NA, aes(ymin = `20%`, ymax = `80%`)) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme_minimal()
}

