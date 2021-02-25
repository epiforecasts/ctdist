# Packages ----------------------------------------------------------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("socialmixr")
library("readxl")
library("janitor")
library("here")
library("lubridate")

# Load data ---------------------------------------------------------------
# load and clean pillars data
ep_raw <- readRDS(here("data-raw", "english_pillars_raw.rds")) %>%
  filter(!is.na(p2ch1cq), p2ch1cq > 0, p2ch1cq < 30,
         !is.na(p2ch2cq), p2ch2cq > 0, p2ch2cq < 30,
         !is.na(p2ch3cq), p2ch3cq > 0 | p2ch3cq < 30) %>%
  mutate(sgene_result = if_else(sgtf == 0, "positive", "negative"))

# load vaccination data
vacc_raw <- readRDS(here("data-raw", "english_vaccinations_raw.rds"))

# load population data
en_pop <- read_excel(here("data-raw", "uk_pop.xls"),
                     sheet = "MYE1", skip = 5) %>%
  janitor::clean_names() %>%
  rename(age_group = x1) %>%
  filter(grepl("^[0-9]", age_group)) %>%
  mutate(lower_age_limit = as.integer(sub("[ -].*$", "", age_group))) %>%
  select(-age_group) %>%
  select(lower_age_limit, pop = england)

# load next level population data
en_pop_l <- read_excel(here("data-raw", "uk_pop.xls"),
                       sheet = "MYE2 - Persons", skip = 4) %>%
  janitor::clean_names() %>%
  filter(!is.na(name)) %>%
  select(-all_ages) %>%
  pivot_longer(starts_with("x"), names_to = "age", values_to = "pop") %>%
  mutate(age = as.integer(sub("^x", "", age)))

# categorise age 
max_age <- max(ep_raw$age, na.rm = TRUE)
age_limits <- c(seq(0, 59, by = 20), seq(60, max_age, by = 10))
if (max_age < max(age_limits)) {
  age_limits <- age_limits[-length(age_limits)]
}
age_limits <- age_limits[age_limits <= 80]


# Summarise in plots ------------------------------------------------------
# apply age categorisation and summarise ct values
iep <- ep_raw %>%
  filter(!is.na(age), date_specimen >= "2020-12-01") %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  group_by(date_specimen, age_group) %>%
  summarise(ct = mean(p2ch1cq), n = n(), .groups = "drop") %>%
  filter(date_specimen < max(date_specimen) - 3)

# plot ct by age
p <- ggplot(iep, aes(x = date_specimen, y = ct,
                     colour = age_group, fill = age_group)) +
  geom_smooth() +
  geom_point(alpha = 0.25, shape = 19) +
  xlab("") +
  ylab("Mean CT value") +
  theme_minimal() +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2")

ggsave(here("figures", "ct_age.pdf"), p, width = 7, height = 5)

# plot proportion vaccinated by age group
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

p <- ggplot(ivacc, aes(x = vaccination_date + 14, y = cum_prop,
                       colour = age_group, fill = age_group)) +
  geom_line() +
  geom_point(alpha = 0.25, shape = 19) +
  xlab("") +
  ylab("Proportion vaccinated") +
  theme_minimal() +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  ylim(c(0, 1))

ggsave(here("figures", "vacc.pdf"), p, width = 7, height = 5)


# Plot summary for lower level geography ----------------------------------
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


# Save data ---------------------------------------------------------------
# raw covariate data + cts
ep_raw_vacc <- ep_raw %>%
  mutate(lower_age_limit =
           socialmixr::reduce_agegroups(age, age_limits),
         age_group =
           socialmixr::limits_to_agegroups(lower_age_limit)) %>%
  inner_join(lvacc %>%
               mutate(date_specimen = vaccination_date + 14 + 7),
             by = c("date_specimen", "ltla_code", "age_group")) %>%
  mutate(time = as.integer(date_specimen - min(date_specimen)))

saveRDS(ep_raw_vacc, here("data", "ct_covariates.rds"))

# summarised ct and vaccination data by LTLA
ep_raw_vacc_mean <- ep_raw_vacc %>%
  mutate(week_specimen = floor_date(date_specimen, "week", 1)) %>%
  group_by(week_specimen, ltla_code, age_group) %>%
  summarise(ct = mean(p2ch1cq), vacc = min(Second), .groups = "drop") %>%
  ungroup() %>%
  group_by(week_specimen, ltla_code, age_group)

saveRDS(ep_raw_vacc_mean, here("data", "ct_summarised.rds"))

# overall vaccine coverage by age
saveRDS(ivacc, here("data", "national_vacc_coverage.rds"))



