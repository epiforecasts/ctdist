# Packages -----------------------------------------------------------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("socialmixr")
library("readxl")
library("janitor")
library("here")
library("mgcv")
library("lubridate")
library("patchwork")


# Load data ---------------------------------------------------------------
# load and clean vaccination data
vacc_raw <- readRDS(here("data-raw", "english_vaccinations_raw.rds"))
en_pop <- read_excel(here::here("repos", "ctdist","uk_pop.xls"),
                     sheet = "MYE1", skip = 5) %>%
  janitor::clean_names() %>%
  rename(age_group = x1) %>%
  filter(grepl("^[0-9]", age_group)) %>%
  mutate(lower_age_limit = as.integer(sub("[ -].*$", "", age_group))) %>%
  select(-age_group) %>%
  select(lower_age_limit, pop = england)

# load and clean population data
en_pop_l <- read_excel(here("data-raw","uk_pop.xls"),
                       sheet = "MYE2 - Persons", skip = 4) %>%
  janitor::clean_names() %>%
  filter(!is.na(name)) %>%
  select(-all_ages) %>%
  pivot_longer(starts_with("x"), names_to = "age", values_to = "pop") %>%
  mutate(age = as.integer(sub("^x", "", age)))

# categorise age (consider using as continuous?)
max_age <- max(ep_raw$age, na.rm = TRUE)
age_limits <- c(seq(0, 59, by = 20), seq(60, max_age, by = 10))
if (max_age < max(age_limits)) {
  age_limits <- age_limits[-length(age_limits)]
}
age_limits <- age_limits[age_limits <= 80]

# apply age groups
iep <- ep_raw %>%
  filter(!is.na(age), date_specimen >= "2020-12-01") %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(date_specimen, age_group) %>%
  summarise(ct = mean(p2ch1cq), n = n(), .groups = "drop") %>%
  filter(date_specimen < max(date_specimen) - 3)

p_ct <- ggplot(iep, aes(x = date_specimen, y = ct,
                     colour = age_group, fill = age_group)) +
  geom_smooth() +
  geom_point(alpha = 0.25, shape = 19) +
  xlab("") +
  ylab("Mean CT value") +
  theme_minimal() +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  coord_cartesian(xlim = c(lubridate::dmy("01-12-2020"), lubridate::dmy("21-02-2021")))

# ggsave(here::here("figure", "ct_age.pdf"), p, width = 7, height = 5)

en_pop_g <- en_pop %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(lower_age_limit, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(age_group) %>%
  summarise(pop = sum(pop), .groups = "drop")

ivacc <- vacc_raw %>%
  filter(!is.na(age), string_dose_number == "First") %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(vaccination_date, age_group) %>%
  tally() %>%
  ungroup() %>%
  arrange(vaccination_date, age_group) %>%
  group_by(age_group) %>%
  mutate(cum = cumsum(n)) %>%
  ungroup() %>%
  inner_join(en_pop_g, by = "age_group") %>%
  mutate(cum_prop = cum / pop)

p_vacc <- ggplot(ivacc, aes(x = vaccination_date + 14, y = cum_prop,
                       colour = age_group, fill = age_group)) +
  geom_line() +
  geom_point(alpha = 0.25, shape = 19) +
  xlab("") +
  ylab("Proportion vaccinated") +
  theme_minimal() +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  ylim(c(0, 1)) +
  coord_cartesian(xlim = c(lubridate::dmy("01-12-2020"), lubridate::dmy("21-02-2021")))

# ggsave(here::here("figure", "vacc.pdf"), p, width = 7, height = 5)

en_pop_lg <- en_pop_l %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(code, name, geography1, lower_age_limit, age_group) %>%
  summarise(pop = sum(pop), .groups = "drop")

lvacc <- vacc_raw %>%
  filter(!is.na(age)) %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(vaccination_date, age_group, ltla_code, string_dose_number) %>%
  tally() %>%
  ungroup() %>%
  arrange(vaccination_date, age_group, ltla_code, string_dose_number) %>%
  group_by(age_group, ltla_code, string_dose_number) %>%
  mutate(cum = cumsum(n)) %>%
  ungroup() %>%
  inner_join(en_pop_lg %>%
               rename(ltla_code = code), by = c("age_group", "ltla_code")) %>%
  mutate(cum_prop = cum / pop) %>%
  select(vaccination_date, age_group, ltla_code, string_dose_number,
         cum_prop) %>%
  pivot_wider(names_from = string_dose_number, values_from = cum_prop) %>%
  replace_na(list(First = 0, Second = 0))

ep_raw_vacc <- ep_raw %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  inner_join(lvacc %>%
               mutate(date_specimen = vaccination_date + 14 + 7),
             by = c("date_specimen", "ltla_code", "age_group")) %>%
  mutate(time = as.integer(date_specimen - min(date_specimen)))

ep_raw_vacc_mean <- ep_raw_vacc %>%
  mutate(week_specimen = floor_date(date_specimen, "week", 1)) %>%
  group_by(week_specimen, ltla_code, age_group) %>%
  summarise(ct = mean(p2ch1cq), vacc = min(Second), .groups = "drop") %>%
  ungroup() %>%
  group_by(week_specimen, ltla_code, age_group)

ggplot(ep_raw_vacc_mean, aes(x = vacc, y = ct, colour = age_group)) +
  geom_jitter() +
  scale_colour_brewer(palette = "Dark2")

fit <- bam(p2ch1cq ~
             age_group + sgene_result + s(time) +
             First * sgene_result + Second * sgene_result,
           data = ep_raw_vacc,
           family = gaussian(link = "identity"))

ep_raw_vacc %>%
  select(p2ch1cq, age_group, sgene_result, time, First, Second) %>%
  group_by(sgene_result, age_group, time) %>%
  summarise(n()) %>%
  filter(age_group %in% c("80+", "70-79"), time > 40) %>%
  ggplot(aes(x = time, y = `n()`)) +
  geom_point() +
  facet_grid(sgene_result ~ age_group, scales = "free")

ep_raw_vacc %>%
  select(p2ch1cq, age_group, sgene_result, time, First, Second) %>%
  group_by(sgene_result, age_group, time) %>%
  summarise(med = median(p2ch1cq),
            uq = quantile(p2ch1cq, 0.975),
            lq = quantile(p2ch1cq, 0.025)) %>%
  ggplot(aes(x = time, y = med)) +
  geom_errorbar(aes(ymin = lq, ymax = uq)) + 
  geom_point() + 
  geom_smooth() +
  facet_grid(sgene_result ~ age_group, scales = "free_y") +
  labs(y = "Median Ct value")


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
source("~/repos/EpiEpi/R/lengthscale_prior.R")
lsp <- get_lengthscale_prior(floor(dat$t/3), dat$t)
dat$lengthscale_alpha <- lsp$alpha
dat$lengthscale_beta <- lsp$beta

library(rstan)
mod <- stan_model(here::here("ctdist", "attempt.stan"))

vacc_grp <- ivacc %>%
  select(vaccination_date, age_group, cum_prop) %>%
  group_by(age_group) %>%
  mutate(time = 0:(n() - 1)) %>%
  filter(time %in% dat$time)

dat$vacc_cov <- matrix(data = vacc_grp$cum_prop, nrow = dat$AG, ncol = dat$t)


res <- rstan::sampling(mod, data = dat)
