# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# Code for Simulating Weekly Flooding Epi Analysis
#
# We will create conditional Quasi-Poisson selecting on county and week. 
# We use two-stage selection with sliding strata on years before and after a 
# flood event, following the proces defined as in the Aggarwal et al. 
# pre-print: https://arxiv.org/abs/2309.13142
#
# P4:
#       Setting up data for simulated epi analysis.
#             1) Create Dummy Exposure Data
#                   i) Bring in dlnm base data and restructure for flood use
#             2) Link Dummy Flood Ocurrence to Exposure Data
#                   i) Convert 1) to a function for easier change of inputs
#                      like case count/flood effect/etc.
#             3) Add lag effects to Dummy Data, and prepare sliding windows.
#                   i) Add RR_lag to dummy data and carry through to 
#                      flood estimation.
#                   ii) Prepare sliding windows for case crossover approach.
#             4) Add cross basis, prepare other scripts to run from source.
#                   i) This script will now be updated based on parameters provided
#                   in script 04_run.R. To ensure those inputs are not over-
#                   written, the lines of code setting the RR, samp_set,
#                   rep_prob, etc. are replaced to take in the variables 
#                   provided in script 4.
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

# Investigate dummy data from dlnm package. The data includes a daily time 
# series, with death counts (including cause-specific counts) and different
# environmental variables (temp, dewpt, humidity, pm2.5, o3). We will use this 
# structure but replace the input to be suitable to changes in parameters for 
# a baseline sample size, and for flood-based exposure inputs.
#
head(chicagoNMMAPS)
summary(chicagoNMMAPS)

# Set parameters for input data, baseline case number, variance in case count,
# baseline year, RR (year_beta) to be used for time trend and RR (RR) to be used
# for flood effect. We are setting as variables here so that we can convert
# the dataset creation to a function moving forward.
# 
# Shifting from p1, we will now create our dataset using a function so we can 
# more easily change our inputs and receive an updated dataframe. We first
# create the function and then apply our inputs
#
set_cases <- function(df, 
                      baseline = 10000, 
                      variance = 500,
                      baseline_yr = 1987,
                      year_beta = 1.2,
                      RR = 2,
                      RR_lag = 1.3) {
  
  # ******
  # df <- chicagoNMMAPS  
  # baseline = 100 
  # variance = 5
  # baseline_yr = 1987
  # year_beta = 2.0
  # RR = 1
  # RR_lag = 1.3
  # ******
  
  df <- df %>% arrange(date)
  df_l <- split(df, f = df$year)
  n_years <- length(df_l)
  
  # Loop through the different years to
  for(yr_i in 1:n_years) {
    
    # Subset to single year as applicable
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
  
  return(df)
  
}

# We will create data for 2 counties. They have the same population, variance,
# and time trend, but with county 2 having a slightly stronger effect of flooding.
# (This will just create the dataset, the flood effect is not applied here)
#
x1 <- set_cases(chicagoNMMAPS, 
                baseline = bl,
                variance = var, 
                year_beta = ybeta_1, 
                RR = RR_1,
                RR_lag = RR_1_lag)

x2 <- set_cases(chicagoNMMAPS, 
                baseline = bl,
                variance = var, 
                year_beta = ybeta_2, 
                RR = RR_2, 
                RR_lag = RR_2_lag)


# The resulting df has the necessary components to now apply a flood effect.
#       We have the date time series.
#       We have death_true, v1, and v2 which together represent the baseline,
#           and the variance (increase or decrease)
#       With the death_true, v1, and v2 we have the number of deaths.
#
# The next step is the application of flood weeks to some set of the data.
# create a fake exposure that represents floods
# within each county and year there are 2 floods
# potentially they crowd, meaning they occur in back-to-back weeks
#
# Add county names
#
x1$county <- 'CountyA'
x2$county <- 'CountyB'

df <- rbind(x1, x2)

df_struct <- unique(df[, c('year', 'county')])
df_struct

# First we will set a list so that we can loop throguh every row in the 
# dataframe and assign some flood days
#
flood_days_l <- vector('list', nrow(df))

# We will also set some parameters here that we can adjust to affect the timing
# of flood occurrence.
# 
# The below is a random pool from which the number of floods in a given year 
# will be selected. With the below, there is an equal likelihood of 0, 1, 2, 
# or 3 flood weeks in a given year
#
# This is now passed by script 04_run

# The rep_prob is the likelihood of a second flood occurring the week following
# the first. Multiple consecutive floods can be useful for testing to ensure
# the model is effectiveyl accounting for these events
#
# This is now passed by script 04_run

for(i in 1:nrow(df_struct)) {
  
  # Use the below to test your code without running the full loop
  # ******
  # i = 1
  # ******
  
  # pick N flood weeks in year. 
  #
  is_flood_week_per_year <- sample(samp_set, 1)
  
  # If there are any floods:
  # When does the first one occur?
  # For simplicity we will avoid going around the year breakpoint
  #
  dt_range <- 40:300
  flood_day <- sample(dt_range, 1)
  flood_days <- flood_day
  trac_i <- 1
  
  # if there are multiple floods, we will go through the loop below and continue
  # adding new days to the list of flood days in the year
  #
  while(is_flood_week_per_year > 1) {
    
    # Determine if the next week is a flood (1) or not a flood (0) based on
    # the probability supplied above
    is_next <- sample(c(1, 0), 1, prob = c(rep_prob, 1- rep_prob))
    
    # if its next, just add 7 to first flood day or -7 if its > 350. This is
    # why we want to avoid wrapping around the year
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

# Form dataframe with the flood days and years. Mark as flood days
#
flood_days_df <- do.call(rbind, flood_days_l)
flood_days_df$is_flood_day <- T

# Remove any duplicates
#
flood_days_df <- flood_days_df %>%
  group_by(flood_day, year, county) %>%
  slice(1)

# Join to the dataset with cases
#
df_w_exp <- left_join(df, flood_days_df,
                      by = join_by(doy == flood_day, county == county,
                                   year == year))

# Confirm the same number of rows with added exposure data and view structure
#
stopifnot(nrow(df_w_exp) == nrow(df))
head(df_w_exp)

# Add markers for date/week/year as time from baseline
#
earliest_day <- min(df_w_exp$date)
df_w_exp$date_int <- as.vector(df_w_exp$date - earliest_day)
df_w_exp$week_num <- floor(df_w_exp$date_int / 7)
df_w_exp$year_num <- floor(df_w_exp$date_int / 365)
head(df_w_exp)

# For this analysis we will want to work as a weekly analysis. To do this we 
# will use dplyr to group by county and week, add an integer for flood marker
# and label the week as flooded or not
#
df_weekly <- df_w_exp %>%
  mutate(is_flood_day_int = 1*(!is.na(is_flood_day))) %>%
  group_by(county, week_num, RR, RR_lag) %>%
  summarize(
    .groups = 'keep',
    
    # Track number of days in week (365 days in a year, some weeks will have less)
    n_days_in_this_week = n(),
    
    # Add marker to track date and link additional data
    week_start = min(date),
    
    # Mark flood events as one
    is_flood_week = (sum(is_flood_day_int) > 0) * 1,
    
    # Note number of cases from data
    case_baseline = sum(death),
    
    # Example of other covariates
    avg_tmp = mean(temp, na.rm = T)
  )


# Remove any week that doesn't have 7 days
#
df_weekly <- df_weekly %>% filter(n_days_in_this_week == 7)

# Arrange dataset by county and week number. This will be especially relevant
# when we set up our crossbasis to run a distributed lag model
#
df_weekly <- df_weekly %>% arrange(county, week_num)

# Set a row ID
# This is important for the strata cross-basis later
#
df_weekly$row_id <- 1:nrow(df_weekly)

# Reintegrate year
#
df_weekly$year <- year(df_weekly$week_start)

# Now add the case spikes for RR for flooding
# We will split up the datafram by county and loop counties
#
df_weekly_l <- split(df_weekly, f = df_weekly$county)
N_COUNTIES <- length(df_weekly_l)

# Loop through counties to apply flood effect
#
for(county_i in 1:N_COUNTIES) {
  
  # Subset to county and arrange by date
  this_df_w_l <- df_weekly_l[[county_i]] %>% arrange(week_start)
  
  # Set n cases to baseline
  this_df_w_l$n_cases <- this_df_w_l$case_baseline
  
  # Add lags (for is_flood)
  this_df_w_l$lag1 <- lag(this_df_w_l$is_flood_week, n = 1)
  this_df_w_l$lag2 <- lag(this_df_w_l$is_flood_week, n = 2)
  this_df_w_l$lag3 <- lag(this_df_w_l$is_flood_week, n = 3)
  this_df_w_l$lag4 <- lag(this_df_w_l$is_flood_week, n = 4)
  
  for(i in 1:nrow(this_df_w_l)) {
    
    # Apply main effect
    if (this_df_w_l$is_flood_week[i] == 1) {
      this_df_w_l$n_cases[i] = this_df_w_l$RR[i] * this_df_w_l$case_baseline[i] 
    } else {
      this_df_w_l$n_cases[i] = this_df_w_l$case_baseline[i]
    }
    
    # Apply lag effect
    if(i > 4) {
      for (j in c(1:4)) {
        if(this_df_w_l[[paste0("lag", j)]][i] == 1) {
          this_df_w_l$n_cases[i] = this_df_w_l$RR_lag[i] * this_df_w_l$n_cases[i] 
        }
      }
    }
    
  }
  
  df_weekly_l[[county_i]] <- this_df_w_l
}

# Rebind county inputs
#
df_weekly <- do.call(rbind, df_weekly_l)

# Now we can view our flood effect!
#
ggplot(df_weekly) +
  geom_point(aes(x = week_start, y = n_cases, col = factor(is_flood_week))) +
  facet_wrap(~county) +
  labs(title = "Time Series of Weekly Case Counts in Two Fictional Counties",
       x = "Week",
       y = "Number of Cases",
       col = "Flood Week (1 = Flooded)") +
  theme_minimal()

# We can see a clear flood effect here. Now the effect looks a lot noisier - and
# really for CountyA where there is a lag effect. Let's Zoom in and see what the
# flood effect looks like.
#
df_weekly$lagflag <- ifelse(df_weekly$lag1 == 1 |
                              df_weekly$lag2 == 1 |
                              df_weekly$lag3 == 1 |
                              df_weekly$lag4 == 1, 1, 0)

# Add plot
#
ggplot(df_weekly[df_weekly$county == "CountyA" & df_weekly$year %in% c(1991:1993),]) +
  geom_point(aes(x = week_start, y = n_cases, col = factor(is_flood_week), size = factor(lagflag))) +
  labs(title = "Time Series of Weekly Case Counts in Two Fictional Counties",
       x = "Week",
       y = "Number of Cases",
       col = "Flood Week (1 = Flooded)",
       size = "Lag Week (1= Lag)") +
  theme_minimal() +
  facet_wrap(~year)

# We can see here that the effect is a result of the lag.
#




