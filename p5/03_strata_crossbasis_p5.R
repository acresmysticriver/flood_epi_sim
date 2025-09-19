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
#             4) Add cross basis estimation.
#                   i) Script 1 will now be updated based on parameters provided
#                   in script 04_run.R. To ensure those inputs are not over-
#                   written, the lines of code setting the RR, samp_set,
#                   rep_prob, etc. are replaced to take in the variables 
#                   provided in script 4.
#                   ii) No changes to script 2
#                   iii) Create cross basis to estimate lag effect of flooding.
#             5) Update scripts to account for added impact of consecutive
#               flood weeks, add script for cross-county meta regression.
#                   i) We will now pass additional RRs for the impact of
#                   multi-week floods (RR_nflood), as well as for the lag after a 
#                   multi-week flood (RR_lag_nflood).
#                   ii) No change to script 2
#                   iii) Add cross bases for n_recent floods and for categorical
#                   flood definitions
# 
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# CREATE STRATA CROSSBASIS
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# For our analysis, we want to estimate the effect of flooding across a lag
# period. To do this, we will establish a cross basis This script
# is meant to run after script 01 and will require the final output from that
# script (df_weekly) to be run. You can run script 1 or use the below to bring
# in the data:
#
# Note, ensuring dates are in order is needed for the cross basis. Arrange
# the data so as to ensure the set up is correct.
#
df_weekly <- df_weekly %>% arrange(county, week_start)

# We will establish our cross basis with our flood week variable.
#       We provide our flooding variable (is_flood_week)
#       We provide the lag period (4 weeks)
#       We note that the variable should be handled as categorical with a break
#           at 0.99 so that we have flood (1) and noflood (0)
#       We provide that the lag should be accounted for with a 4-degree polynomial
#       Lastly, we note that the cross-basis should be assessed by county
cb.flood.strat <- crossbasis(df_weekly$is_flood_week,
                       lag = 4,
                       argvar=list(fun="strata", breaks = (0.99)),
                       arglag = list(fun = "poly", degree = 4),
                       group = df_weekly$county)

# Now that we have a second lagged variable for multiple consecutive flood weeks,
# we will add additional cross basis. We will include one that is including 
# the n_recent floods over time, and then three separate additional cross basis
# that are yes/no for the different flood categories (single week, two week, 
# three or more weeks)
#
# N recent flood cross basis
#
cb.flood.nflood <- crossbasis(df_weekly$n_recent_floods,
                              lag = 4,
                              argvar=list(fun="lin"),
                              arglag = list(fun = "poly", degree = 4),
                              group = df_weekly$county)

# Categorical cross bases
#
cb.flood.flood1 <- crossbasis(df_weekly$exp1 ,
                              lag = 4,
                              argvar=list(fun="strata", breaks = (0.99)),
                              arglag = list(fun = "poly", degree = 4),
                              group = df_weekly$county)
cb.flood.flood2 <- crossbasis(df_weekly$exp2 ,
                              lag = 4,
                              argvar=list(fun="strata", breaks = (0.99)),
                              arglag = list(fun = "poly", degree = 4),
                              group = df_weekly$county)
cb.flood.flood3 <- crossbasis(df_weekly$exp3 ,
                              lag = 4,
                              argvar=list(fun="strata", breaks = (0.99)),
                              arglag = list(fun = "poly", degree = 4),
                              group = df_weekly$county)

# Assess our cross basis
#
head(cb.flood.strat)
cb.flood.strat
tail(cb.flood.strat)

