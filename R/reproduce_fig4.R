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

# missing packages
library(readr)

# Source utilities --------------------------------------------------------
source(here("R", "ct-rt-utils.R"))

# Set up parallel ---------------------------------------------------------
# number of cores (only uses up to 4 chains)
cores <- 4
# number of threads per core (so available cores / cores)
threads <- 1

## CODE PILLAGED FROM PAPER REPO ##

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))


obs_dat_all <- read_csv(here("data","panther_Ct_20200403-20201110.csv")) %>% rename(panther_Ct=ORF1ab_Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n())
#obs_dat_all <- read_csv("data/panther_Ct_20200403-20201110.csv") %>% rename(panther_Ct=ORF1ab_Ct) %>%
#  mutate(platform="Panther",first_pos=1) %>%
#  mutate(id=1:n())
obs_dat1 <- obs_dat_all

obs_dat1 <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct)

obs_dat_all <- obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  arrange(first_day) %>%
  rename(date = first_day, ct=panther_Ct)

comb_dat <- left_join(obs_dat1, obs_dat_all)
date_key <- distinct(comb_dat %>% dplyr::select(t, date))

date_min_date <- min(date_key$date)
date_min_t <- min(date_key$t)

date_max_date <- max(date_key$date)
date_max_t <- max(date_key$t)

integer_seq_times <- seq(0, date_max_t)
date_seq_times <- seq(date_min_date-date_min_t, date_max_date,by="1 day")
date_key <- tibble(t=integer_seq_times,date=date_seq_times)

# This plot shows data is the same as in Figure 4C

p_dat <- ggplot(obs_dat_all) + 
  geom_violin(aes(x=date,group=date,y=ct),scale="width",fill="grey1",
              alpha=0.5,color="black",size=0.1,
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=date, y=ct,group=date),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data=obs_dat_all %>% group_by(date) %>% summarize(median_ct=median(ct)),
              aes(x=date,y=median_ct),col="blue",se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(42, 10),expand=c(0,0)) +
  geom_hline(yintercept=40,linetype="dashed") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_date(limits=as.Date(c("2020-04-01","2020-11-15")),breaks="1 month",expand=c(0,0)) +
  xlab("Date of sample") +
  ylab("Ct value") +
  labs(tag="C")

## END OF PILLAGED CODE ##

obs_dat_all <- as.data.table(obs_dat_all)
obs_dat_all <- obs_dat_all[, time := date - min(date) + 1]

min_date <- min(obs_dat_all$date, na.rm = TRUE)

# define CT post infection (loosely inspired by Hay et al)
ct <- tibble(mean = c(40 - 0:4*5, 18 + 0:10*2),
             sd = c(rep(1, 5), rep(2, 11)))

# define observations
dat <- stan_data(obs_dat_all, 
                 load_vec = "ct",
                 overall_prob = 1,
                 ct =  ct,
                 dt = 43,
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
fit$cmdstan_summary()

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
