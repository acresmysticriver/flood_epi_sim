
# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# RUN DLNM AND MANUAL LAG ANALYSES
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# For our analysis, we want to estimate the effect of flooding across a lag
# period. In script 1, we set up our dataset. In script 2, we defined our
# case and control samples. In script 3 we set up a cross basis for use in 
# the dlnm model. 
#
# In this script we will source scripts 1-3, where our dataset and analytic
# structure are defined. We will modify the inputs for the functions applied
# in script 1 to assess how changing the flood frequency and clustering (meaning
# floods in back to back weeks) influence our result.
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
#                   iv) Perform initial analysis to evaluate whether system
#                       is robust to identify flooding
#             5) Update scripts to account for added impact of consecutive
#               flood weeks, add script for cross-county meta regression.
#                   i) We will now pass additional RRs for the impact of
#                   multi-week floods (RR_nflood), as well as for the lag after a 
#                   multi-week flood (RR_lag_nflood).
#                   ii) No change to script 2
#                   iii) Add cross bases for n_recent floods and for categorical
#                   flood definitions
#                   iv) Perform updated analysis to evaluate if updates 
#                   effectively allow for assessing impact of consecutive events
#                   v) Run as mvmeta meta-regression across counties
# 
# First, clear the R environment as we will reset the inputs to the scripts
#
rm(list = ls()); gc()

library(mvmeta)

# Set sample from which N floods per year will be set
#
samp_set <- c(0, 1, 2)

# Set likelihood of back-to-back flood weeks
#
rep_prob <- 0.4

# Set baseline and variance
#
bl <- 100000
var <- 500

# COUNTY A
RR_1     <- 1.50  # the RR on lag 0
RR_1_lag <- 1.20  # the RR on lags 1:4
ybeta_1  <- 0.90  # the year trend (in log space, so 2 means doubling every year)
RR_1_nflood <- 1.035  # the RR associated with each additional flood in recent weeks
RR_1_lag_nflood <- 1.02  # the RR associated with each additional flood in lag weeks

# COUNTY B
RR_2     <- 1.60  # the RR on lag 0
RR_2_lag <- 1.08 # the RR on lags 1:4
ybeta_2  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
RR_2_nflood <- 1.035  # the RR associated with each additional flood in recent weeks
RR_2_lag_nflood <- 1.02  # the RR associated with each additional flood in lag weeks

source('flood_epi_sim/p5_final_multiweek_run/01_create_dummy_data_p5.R')
source('flood_epi_sim/p5_final_multiweek_run/02_create_sliding_windows_p5.R')
source('flood_epi_sim/p5_final_multiweek_run/03_strata_crossbasis_p5.R')

# uncomment this if you want to save it as an image
# scale = 1.5
# png("demo_v1.png", units = 'in',
#     width = 15.6/scale,
#     height = 8.5/scale, res = 600)

par(mfrow = c(2, 2))

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# STAGE 1: CREATE A MODEL IN EACH COUNTY
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

COUNTIES <- unique(expanded_df$county)
N_COUNTIES <- length(COUNTIES)

# then resplit based on county
county_l <- split(expanded_df, f = expanded_df$county)

# Check data
#
ggplot(expanded_df) + 
  geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
             position = position_dodge(width = 1)) +
  facet_grid(cols = vars(county))

# Set list to capture coefficients and covariance
#
coef_list <- list()
vcov_list <- list()
 
# loop throguh counties
#
for (county_i in c(1:2)) {
  
  # Subset to county 
  #
  this_df <- county_l[[county_i]]
  dim(this_df)
  
  # Track runs
  #
  run <- paste0(county_i)
  cat(run, "\n")
  
  # Add week iteration
  #
  START_YEAR = 1987
  this_df$week_iter = this_df$week_num + 52*(this_df$year - START_YEAR)
  
  # Build formula to run
  #
  formula_run <- as.formula("n_cases ~ is_flood_week + 
                              lag1 + lag2 + lag3 + lag4 +
                              n_recent_floods +
                               lag_nflood_1  + lag_nflood_2 + lag_nflood_3  + lag_nflood_4  + 
                              avg_tmp +
                              year")
  
  # Run function
  #
  modelrun <- gnm(formula_run, data=this_df, 
                  family=quasipoisson,
                  eliminate = factor(strata),
                  na.action = 'na.exclude')
  
  # Add to list
  #
  coef_list[[run]] <- coef(modelrun)
  vcov_list[[run]]  <- vcov(modelrun)
  
  # Get outputs for model
  #
  t1 <- exp(confint(modelrun))
  t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(modelrun)))
  colnames(t1)[1:2] <- c('lb', 'ub')

  # Add markers
  t1$county <- county_i
  
  # Stack em up
  #
  if (county_i == 1) {
    stack_t <- t1
  } else {
    stack_t <- rbind(stack_t, t1)
  }
}

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# STAGE 2: RUN META-MODEL
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# Extract outputs from strat
#
coefs <- do.call(rbind, coef_list)
vcovs <- vcov_list

# Run meta model
#
meta_model <- mvmeta(coefs ~ 1, S = vcovs, method = "reml")

# View results
#
meta_model$model
summary(meta_model)

# Print output
#
print(exp(coef(meta_model)), digits=4)
print(exp(confint(meta_model)), digits=4)

# Save output with confidence intervals to view estimates and plot result
#
t1 <- exp(confint(meta_model))
t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(meta_model)))
colnames(t1)[1:2] <- c('lb', 'ub')

# Add markers
t1

# Plot results
ggplot(t1) +
  geom_point(aes(y = var, x = est))+
  geom_errorbar(aes(y = var, xmin = lb, xmax = ub))
