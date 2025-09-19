# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# Code for Simulating Weekly Flooding Epi Analysis
#
# We will create conditional Quasi-Poisson selecting on county and week. 
# We use two-stage selection with sliding strata on years before and after a 
# flood event, following the proces defined as in the Aggarwal et al. 
# pre-print: https://arxiv.org/abs/2309.13142
#
# P1:
#       Setting up data for simulated epi analysis.
#             1) Create Dummy Exposure Data
#                   i) Bring in 
# 
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

library(dlnm)
library(tidyverse)
library(lubridate)
library(gnm)
library(splines)
library(future)
library(future.apply)
library(patchwork)
library(zoo)
library(purrr)
library(tidyr)
library(data.table)

set.seed(123)
plan(multisession)

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# CREATE DUMMY EXPOSURE DATA
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# Investigate dummy data from dlnm package
head(chicagoNMMAPS)
summary(chicagoNMMAPS$date)

# reset cases based on a year trend
# This function is set up to take in information about the size of a county, 
# and relative risks of flooding/flood lags
#
set_cases <- function(df, 
                      baseline = 10000, 
                      variance = 500,
                      baseline_yr = 1987,
                      year_beta = 1.2,
                      RR = 2,
                      RR_lag = 1.0,
                      RR_nvis = NULL,
                      RR_lag_nvis = NULL) {
  
  # ******
  # df <- chicagoNMMAPS  
  # baseline = 100 
  # variance = 5
  # baseline_yr = 1987
  # year_beta = 2.0
  # RR = 1
  # ******
  
  df <- df %>% arrange(date)
  df_l <- split(df, f = df$year)
  n_years <- length(df_l)
  
  for(yr_i in 1:n_years) {
    
    # *****
    # yr_i = 1
    # *****
    
    n_rows <- nrow(df_l[[yr_i]])
    
    # need to work in log-space so the coefficients work out
    df_l[[yr_i]]$death_true = 
      # baseline + how much to increase by year
      baseline * (year_beta)^(df_l[[yr_i]]$year - baseline_yr)
    
    # add some random noise
    v1 <- sample(c(-1, 0, 1), size = n_rows, replace = T) 
    v2 <- rpois(n_rows, variance)
    df_l[[yr_i]]$v1 <- v1
    df_l[[yr_i]]$v2 <- v2
    df_l[[yr_i]]$death = df_l[[yr_i]]$death_true + v1 * v2
    
    # make sure its an integer
    df_l[[yr_i]]$death <- round(df_l[[yr_i]]$death)
    
    
  }
  
  df <- do.call(rbind, df_l)
  
  df$RR <- RR
  df$RR_lag <- RR_lag
  df$RR_nvis <- RR_nvis
  df$RR_lag_nvis <- RR_lag_nvis
  
  # ggplot(df) +
  #   geom_point(aes(x = date, y = death))
  
  return(df)
  
}

# create fake data for two counties
# first, lets try with no year trend (so year_beta = 1)
x1 <- set_cases(chicagoNMMAPS, 
                baseline = bl_1,
                variance = var_1, 
                year_beta = ybeta_1, 
                RR = RR_1, 
                RR_lag = RR_1_lag,
                RR_nvis = RR_1_nvis,
                RR_lag_nvis = RR_1_lag_nvis)

x2 <- set_cases(chicagoNMMAPS, 
                baseline = bl_2,
                variance = var_2, 
                year_beta = ybeta_2, 
                RR = RR_2, 
                RR_lag = RR_2_lag,
                RR_nvis = RR_2_nvis,
                RR_lag_nvis = RR_2_lag_nvis)

x1$county <- 'CountyA'
x2$county <- 'CountyB'

df <- rbind(x1, x2)

df_struct <- unique(df[, c('year', 'county')])
df_struct

# create a fake exposure that represents floods
# within each county and year there are 2 floods
# potentially they crowd, meaning they occur in back-to-back weeks

flood_days_l <- vector('list', nrow(df_struct))

for(i in 1:nrow(df_struct)) {
  
  # ******
  # i = 1
  # ******
  
  # pick N flood weeks in year
  is_flood_week_per_year <- sample(samp_set, 1)
  
  # If there are any floods:
  # when does the first one occur
  # IMPORTANT - make it so none happen around the year breakpoint
  dt_range <- 40:300
  flood_day <- sample(dt_range, 1)
  flood_days <- flood_day
  trac_i <- 1
  
  # if there is a second one, is it the next week or not
  while(is_flood_week_per_year > 1) {
    
    # there's a 70% chance its the following week
    is_next <- sample(c(0,1), 1, prob = c(rep_prob, 1- rep_prob))
    
    # if its next, just add 7 to first flood day or -7 if its > 350
    if(is_next == 1) {
      if(flood_day < 350) {
        flood_days = c(flood_days, flood_days[trac_i] + 7)
      } else {
        flood_days = c(flood_days, flood_days[trac_i] - 7)
      }
    } else {
      flood_days <- c(flood_days, sample(dt_range, 1))
    }
    
    is_flood_week_per_year <- is_flood_week_per_year - 1
    trac_i <- trac_i + 1
    
  }
  
  # adding flood_days to this year
  out <- data.frame(flood_day = flood_days)
  out$year <- df_struct$year[i]
  out$county <- df_struct$county[i]
  
  # output to list
  flood_days_l[[i]] <- out
  
}

# to data frame
flood_days_df <- do.call(rbind, flood_days_l)
flood_days_df$is_flood_day <- T

# keep one of each county-date
#
flood_days_df <- flood_days_df %>%
  group_by(flood_day, year, county) %>%
  slice(1)

# join to our dataset by doy
df_w_exp <- left_join(df, flood_days_df,
                      by = join_by(doy == flood_day, county == county,
                                   year == year))
# confirm the same number of rows
stopifnot(nrow(df_w_exp) == nrow(df))
head(df_w_exp)

# ** TRICKY PART #1, how you handle going around the year
earliest_day <- min(df_w_exp$date)
df_w_exp$date_int <- as.vector(df_w_exp$date - earliest_day)
df_w_exp$week_num <- floor(df_w_exp$date_int / 7)
df_w_exp$year_num <- floor(df_w_exp$date_int / 365)
head(df_w_exp)

# create the weekly dataset
df_weekly <- df_w_exp %>%
  mutate(is_flood_day_int = 1*(!is.na(is_flood_day))) %>%
  group_by(county, week_num, RR, RR_lag, RR_nvis, RR_lag_nvis) %>%
  summarize(
    .groups = 'keep',
    # hmm not all weeks have every day
    n_days_in_this_week = n(),
    # 
    week_start = min(date),
    # mark flood events as one
    is_flood_week = (sum(is_flood_day_int) > 0) * 1,
    # cases
    case_baseline = sum(death),
    # example of other covariates
    avg_tmp = mean(temp, na.rm = T)
  )

# we want to also add a variable for the n-vis in recent weeks (up to 4)
df_weekly$n_recent_floods <- NA
for(i in 5:nrow(df_weekly)) {
  df_weekly$n_recent_floods[i] <- sum(df_weekly$is_flood_week[(i-4):i])
}

# only maintain n_recent_floods for flood weeks. We will carry this through
# in lag for showing recent floods in non-weeks. Now we have a variable that
# essentially represents the week of consecutive or near-consecutive flooding
df_weekly$n_recent_floods <- ifelse(df_weekly$is_flood_week == 0, 0, df_weekly$n_recent_floods)
df_weekly$n_recent_floods <- pmax(0, df_weekly$n_recent_floods -1 )

# remove any week that doesn't have 7 days
df_weekly <- df_weekly %>% filter(n_days_in_this_week == 7)

# arrange, because order matters in the crossbasis ...
df_weekly <- df_weekly %>% arrange(county, week_num)

# set a row ID
# this is important for the strata cross-basis later
df_weekly$row_id <- 1:nrow(df_weekly)

# and add year back
df_weekly$year <- year(df_weekly$week_start)

# plot
# ggplot(df_weekly) + 
#   geom_point(aes(x = week_start, y = case_baseline, col = n_recent_floods)) +
#   facet_wrap(~county)

# now add the case spikes for RR, RR_lag, and accounting for RR_nvis
# ASSUMPTIONS:
#         the greatest impacts of flooding will be during and soon after the
#                 flood week
#         more events in recent weeks will lead to greater impacts
#         health impacts will decay over the course of lag time
# 
df_weekly_l <- split(df_weekly, f = df_weekly$county)
N_COUNTIES <- length(df_weekly_l)
for(county_i in 1:N_COUNTIES) {
  
  # Subset to county and arrange by date
  this_df_w_l <- df_weekly_l[[county_i]] %>% arrange(week_start)
  
  # Set n cases to baseline
  this_df_w_l$n_cases <- this_df_w_l$case_baseline
  
  # add lags (for is_flood)
  this_df_w_l$lag1 <- lag(this_df_w_l$is_flood_week, n = 1)
  this_df_w_l$lag2 <- lag(this_df_w_l$is_flood_week, n = 2)
  this_df_w_l$lag3 <- lag(this_df_w_l$is_flood_week, n = 3)
  this_df_w_l$lag4 <- lag(this_df_w_l$is_flood_week, n = 4)
  
  # add lags (for n_vis)
  this_df_w_l$lag_nvis_1 <- lag(this_df_w_l$n_recent_floods, n = 1)
  this_df_w_l$lag_nvis_2 <- lag(this_df_w_l$n_recent_floods, n = 2)
  this_df_w_l$lag_nvis_3 <- lag(this_df_w_l$n_recent_floods, n = 3)
  this_df_w_l$lag_nvis_4 <- lag(this_df_w_l$n_recent_floods, n = 4)
 
  for(i in 1:nrow(this_df_w_l)) {
    
    # To more accurately consider, for county 2 randomly add noise to the RR
    #
    #if (county_i == 2) {
    #  mult <- sample(c(0.5, 0.7, 1.4, 2), 1)
    #  this_df_w_l$RR[i] <- this_df_w_l$RR[i] * mult
    #  this_df_w_l$RR_lag[i] <- this_df_w_l$RR_lag[i] * mult
    #}
    
    ## main effect
    if (this_df_w_l$is_flood_week[i] == 1) {
      this_df_w_l$n_cases[i] = this_df_w_l$RR[i] * this_df_w_l$case_baseline[i] 
    } else {
      this_df_w_l$n_cases[i] = this_df_w_l$case_baseline[i]
    }
     
    if(i > 4) {
        for (j in c(1:4)) {
          if(this_df_w_l[[paste0("lag", j)]][i] == 1) {
            this_df_w_l$n_cases[i] = this_df_w_l$RR_lag[i] * this_df_w_l$n_cases[i] 
          }
        }
      }
      
      if(i > 2*4) {
        
        if (this_df_w_l$n_recent_floods[i] >= 1) {
          this_df_w_l$n_cases[i] = this_df_w_l$n_cases[i] +
            ((this_df_w_l$RR_nvis[i] ^ this_df_w_l$n_recent_floods[i]) - 1) * this_df_w_l$n_cases[i]
        } else {
          this_df_w_l$n_cases[i] = this_df_w_l$n_cases[i]
        }
        
        for (j in c(1:4)) {
          if(this_df_w_l[[paste0("lag_nvis_", j)]][i] >= 1) {
            this_df_w_l$n_cases[i] = this_df_w_l$n_cases[i] +
              ((this_df_w_l$RR_lag_nvis[i] ^ this_df_w_l[[paste0("lag_nvis_", j)]][i]) - 1) * this_df_w_l$n_cases[i] 
          }
        }
      }
  }
    
  df_weekly_l[[county_i]] <- this_df_w_l
}

df_weekly <- do.call(rbind, df_weekly_l)

# Apply method that will characterize types of floods
# Add flag for length of consecutive flood weeks
#
df_weekly <- df_weekly  %>%
  arrange(county, year, week_num) %>%
  group_by(county, year) %>%
  mutate(
    flood_run = ifelse(is_flood_week==1, rleid(is_flood_week), NA)
  ) %>% ungroup()

## Assign week_within_flood
df_weekly <- df_weekly %>%
  group_by(county, year, flood_run) %>%
  mutate(week_within_flood = ifelse(is_flood_week==1, row_number(), NA)) %>%
  ungroup()

## Make categorical: NoFlood, Flood_Week1, Flood_Week2, Flood_Week3plus
df_weekly <- df_weekly %>%
  mutate(
    flood_cat = case_when(
      is.na(week_within_flood) ~ "NoFlood",
      week_within_flood == 1 ~ "Flood_Week1",
      week_within_flood == 2 ~ "Flood_Week2",
      week_within_flood >= 3 ~ "Flood_Week3plus"
    ),
    flood_cat = factor(flood_cat, 
                       levels = c("NoFlood","Flood_Week1","Flood_Week2","Flood_Week3plus"))
  )

df_weekly$exp1 <- (df_weekly$flood_cat == "Flood_Week1")*1
df_weekly$exp2 <- (df_weekly$flood_cat == "Flood_Week2")*1
df_weekly$exp3 <- (df_weekly$flood_cat == "Flood_Week3plus")*1

