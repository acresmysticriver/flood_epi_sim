
library(mvmeta)

# SET SWITCHES FOR TESTING
# -- a 4 week lag (so 0, 1, 2, 3, 4) is hard-coded

# clear
rm(list = ls()); gc()

# COUNTY A
RR_1     <- 2.00  # the RR on lag 0
RR_1_lag <- 1.50  # the RR on lags 1:4
RR_1_nvis <- 3  # the RR associated with flood on day and in lag
ybeta_1  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl <- 10000       # n case baseline 
var <- 500        # variance in case by week

# COUNTY B
RR_2     <- 3.00  # the RR on lag 0
RR_2_lag <- 1.20 # the RR on lags 1:4
RR_2_nvis <- 6  # the RR associated with flood on day and in lag
ybeta_2  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl <- 100       # n case baseline 
var <- 75        # variance in case by week

source('01_create_dummy_data.R')
source('02_create_sliding_windows.R')
source('03_strata_crossbasis.R')

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
  
  # loop through model types
  #
  for (mod_type in c("double_flood")) {
    
    # Track runs
    #
    run <- paste0(county_i, "_", mod_type)
    cat(run, "\n")
    
    # *****
    # tricky part #1
    # in the model also control for time
    # - ns(Date, with 2 knots per year since we are doing weekly)
    # - but you can't do by day since we are aggregating to week
    START_YEAR = 1987
    this_df$week_iter = this_df$week_num + 52*(this_df$year - START_YEAR)
    
    # Specify type oformula from loop
    #
    if (mod_type == "year") {
      formula_run <- as.formula("n_cases ~ is_flood_week + 
                              lag1 + lag2 + lag3 + lag4 +
                              avg_tmp + year")
    } else if (mod_type == "spline") {
      formula_run <- as.formula("n_cases ~ is_flood_week + 
                              lag1 + lag2 + lag3 + lag4 +
                              avg_tmp +
                              ns(week_iter, df = 4)")
    } else if (mod_type == "doubleflood") {
      formula_run <- as.formula("n_cases ~ is_flood_week + 
                              lag1 + lag2 + lag3 + lag4 +
                              avg_tmp + year + double_flood")
    }
    
    # Run function
    #
    modelrun <- gnm(formula_run, data=this_df, # include mean of tmax for week as sensitivity, check what was done
                    family=quasipoisson,
                    eliminate = factor(strata),
                    na.action = 'na.exclude')
    
    # Add to list
    #
    coef_list[[run]] <- coef(modelrun)
    vcov_list[[run]]  <- vcov(modelrun)
    
    mean_var <- mean(vcov(modelrun))
    
    # AWESOME, this works
    # you see that the is_case_period coefficent and year reflects
    # what you set it as 
    
    t1 <- exp(confint(modelrun))
    t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(modelrun)))
    colnames(t1)[1:2] <- c('lb', 'ub')
    t1
    
    # Add markers
    t1$county <- county_i
    t1$mod <- mod_type
    t1$mean_var <- mean_var
    
    # Stack em up
    #
    if (county_i == 1 & mod_type == "year") {
      stack_t <- t1
    } else {
      stack_t <- rbind(stack_t, t1)
    }
  }
}

# Get metamodel for each model type 
#
for (pattern_in in c("year", "spline")) {
  
  # Get coefficients for a given outcome - strata
  #
  coef_list_n_visits <- coef_list[grepl(pattern_in, names(coef_list))]
  vcov_list_n_visits <- vcov_list[grepl(pattern_in, names(vcov_list))]
  
  # Extract outputs from strat
  #
  coefs <- do.call(rbind, coef_list_n_visits)
  vcovs <- vcov_list_n_visits
  
  # Run meta model
  #
  meta_model <- mvmeta(coefs ~ 1, S = vcovs, method = "reml")
  
  # MODEL FRAME 
  meta_model$model
  summary(meta_model)
  
  # Print output
  #
  print(exp(coef(meta_model)), digits=4)
  print(exp(confint(meta_model)), digits=4)
  
  t1 <- exp(confint(meta_model))
  t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(meta_model)))
  colnames(t1)[1:2] <- c('lb', 'ub')
  
  # Add markers
  t1$type <- pattern_in
  
  # Stack
  #
  if (pattern_in == "year") {
    all_stack_cumul <- t1
  } else {
    all_stack_cumul <- rbind(all_stack_cumul, t1)
  }
  
}
