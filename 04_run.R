
# SET SWITCHES FOR TESTING
# -- a 4 week lag (so 0, 1, 2, 3, 4) is hard-coded

# clear
rm(list = ls()); gc()

# COUNTY A
RR_1     <- 2.00  # the RR on lag 0
RR_1_lag <- 1.50  # the RR on lags 1:4
RR_1_nvis <- 1.03  # the RR associated with each additional flood in recent weeks
RR_1_lag_nvis <- 1.01  # the RR associated with each additional flood in lag weeks
ybeta_1  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl <- 10000       # n case baseline 
var <- 500        # variance in case by week

# COUNTY B
RR_2     <- 2.00  # the RR on lag 0
RR_2_lag <- 1.00 # the RR on lags 1:4
RR_2_nvis <- 1.03  # the RR associated with each additional flood in recent weeks
RR_2_lag_nvis <- 1.02  # the RR associated with each additional flood in lag weeks
ybeta_2  <- 1.00  # the year trend (in log space, so 2 means doubling every year)

bl <- 10000       # n case baseline 
var <- 500        # variance in case by week

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

county_i = 1

this_df <- county_l[[county_i]]
dim(this_df)

# first remove any empty strata
# this actually is redundant because of the way we are defining strata
# but good to check nonetheless
this_df_agg <- this_df %>%
  group_by(strata) %>%
  summarize(
    .groups = 'keep',
    n_weeks_in_strata = sum(is_flood_week)
  ) %>%
  mutate(keep = 1)

this_df <- left_join(this_df, this_df_agg,
                     by = join_by(strata))
dim(this_df)


# Subset to single events and those with event and event across lag
#
this_df <- this_df %>%
  group_by(strata) %>%
  mutate(max_n_rec = max(n_recent_floods))
this_df_single <- this_df[this_df$max_n_rec == 1, ]
this_df_multi <- this_df[this_df$max_n_rec > 1, ]


# test plot
 ggplot(this_df# [this_df$strata %in% c("CountyA_Flood01", "CountyA_Flood02") & this_df$year == 1987,]
        ) + 
   geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
              position = position_dodge(width = 1)) +
   facet_grid(cols = vars(max_n_rec))

 
# re-estimate case to make sure it matches expectation
 this_df_multi$flood_effect <- pmax(this_df_multi$case_baseline * this_df_multi$RR * this_df_multi$is_flood_week - this_df_multi$case_baseline, 0)
 this_df_multi$lag1_effect <- pmax(this_df_multi$case_baseline * this_df_multi$RR_lag * this_df_multi$lag1 - this_df_multi$case_baseline , 0)
 this_df_multi$lag2_effect <- pmax(this_df_multi$case_baseline * this_df_multi$RR_lag * this_df_multi$lag2 - this_df_multi$case_baseline, 0) 
 this_df_multi$lag3_effect <- pmax(this_df_multi$case_baseline * this_df_multi$RR_lag * this_df_multi$lag3 - this_df_multi$case_baseline, 0) 
 this_df_multi$lag4_effect <- pmax(this_df_multi$case_baseline * this_df_multi$RR_lag * this_df_multi$lag4 - this_df_multi$case_baseline, 0) 
 
 # 
 this_df_multi$total_case_exp <- this_df_multi$case_baseline + this_df_multi$flood_effect +
   this_df_multi$lag1_effect +
   this_df_multi$lag2_effect +
   this_df_multi$lag3_effect +
   this_df_multi$lag4_effect
 
 this_df_multi$case_diff <-  this_df_multi$total_case_exp - this_df_multi$n_cases
 
# *****
# tricky part #1
# in the model also control for time
# - ns(Date, with 2 knots per year since we are doing weekly)
# - but you can't do by day since we are aggregating to week
START_YEAR = 1987
this_df$week_iter = this_df$week_num + 52*(this_df$year - START_YEAR)
# *****

# and covariates
# - avg_tmp

# *****
# ok here's the other tricky part, 
# get the cross_basis externally based on the row_ids
# you are doing this so you can still do 
# cross-reduce on the model afterwards, 
# i think this should work 
cb.flood.nflood_local <- cb.flood.nflood[this_df$row_id, ]
cb.flood.strata_local <- cb.flood.strat[this_df$row_id, ]

# you have to transfer the attributes so crosspred works later
attr_transfer <- c('df', 'range', 'lag', 'argvar', 
                   'arglag', 'group', 'class')
for(att in attr_transfer) {
  attr(cb.flood.nflood_local, att) <- attr(cb.flood.nflood, att)
  attr(cb.flood.strata_local, att) <- attr(cb.flood.strat, att)
}
# *****

# make the model
N_YEARS = length(unique(this_df$year))

# Add classifier for whether a flood is multiweek within strata
#
this_df <- this_df %>%
  group_by(strata) %>%
  mutate(n_floodweeks = sum(is_flood_week),
         n_lag1flood = sum(lag1),
         n_lag2flood = sum(lag2),
         n_lag3flood = sum(lag3),
         n_lag4flood = sum(lag4))

# Mod v0 - use no cross basis
# **************
mod.v0 <- gnm(n_cases ~ is_flood_week + 
               # manually do lags
               lag1 + lag2 + lag3 + lag4 + 
               # account for covariates
               avg_tmp +
               # ********
               # and yearly time - I think just the decade, so just 2 knots?
               # can iterate here, uncertain what is best
               year,
             # ********
             data = this_df_multi,
             family = quasipoisson,
             eliminate = factor(strata),
             subset = keep == 1,
             na.action = 'na.exclude')

summary(mod.v0)

# AWESOME, this works
# you see that the is_case_period coefficent and year reflects
# what you set it as 

t0 <- exp(confint(mod.v0))
t0 <- data.frame(t0, var = row.names(t0), est = exp(coef(mod.v0)))
colnames(t0)[1:2] <- c('lb', 'ub')
t0

# Reverse the order of y-axis categories (to mimic ggplot behavior)
# Set up blank plot with appropriate limits
par(mar = c(5, 9, 4, 2))
y_vals <-1:nrow(t0)
plot(NULL,
     xlim = range(c(t0$lb, t0$ub)),
     ylim = range(1:nrow(t0)),
     yaxt = "n",
     xlab = "Estimate",
     ylab = "",
     main = "Base model")

# Add y-axis labels
axis(2, at = y_vals, labels = t0$var, las = 1)

# Add segments (error bars)
segments(x0 = t0$lb, x1 = t0$ub, y0 = y_vals, y1 = y_vals)

# Add points (estimates)
points(t0$est, y_vals, pch = 5, cex = 0.8)




# Mod v1 - use flood alone as cross basis
# **************
mod.v1 <- gnm(n_cases ~ cb.flood.strata_local + n_recent_floods +
                # account for covariates
                avg_tmp +
                # ********
                # and yearly time - I think just the decade, so just 2 knots?
                # can iterate here, uncertain what is best
                year,
              # ********
              data = this_df,
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')

summary(mod.v1)

# Ok this seems to work now as well
cp.strata <- crosspred(cb.flood.strata_local, mod.v1)

t1 <- data.frame(var = 'crosspred',
                 lb = cp.strata$allRRlow,
                 ub = cp.strata$allRRhigh,
                 est = cp.strata$allRRfit)

# with lags
par(mfrow = c(1, 2))
plot(cp.strata, 'slices', var = c(1),  main = 'CB with Lag, 1')
plot(cp.strata, 'overall')

# dev.off()


# Mod v2 - use cumulative flood lag as cross basis
# **************
mod.v2 <- gnm(n_cases ~ cb.flood.nflood_local + 
                # account for covariates
                avg_tmp +
                # ********
                # and yearly time - I think just the decade, so just 2 knots?
                # can iterate here, uncertain what is best
                year,
              # ********
              data = this_df,
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')


summary(mod.v2)

# Ok this seems to work now as well
cp_nvis <- crosspred(cb.flood.nflood_local, mod.v2)

t1_nvis <- data.frame(var = 'crosspred',
                 lb = cp_nvis$allRRlow,
                 ub = cp_nvis$allRRhigh,
                 est = cp_nvis$allRRfit)

# with lags
par(mfrow = c(1, 2))
plot(cp_nvis, 'slices', var = c(1),  main = 'CB with Lag, 1')
plot(cp_nvis, 'overall')

t1 <- data.frame(var = 'crosspred',
                 lb = cp_nvis$allRRlow,
                 ub = cp_nvis$allRRhigh,
                 est = cp_nvis$allRRfit)

# with lags
par(mfrow = c(2,3))
plot(cp_nvis, 'slices', lag = c(0),  main = 'CB with Lag, 0')
plot(cp_nvis, 'slices', lag = c(1),  main = 'CB with Lag, 1')
plot(cp_nvis, 'slices', lag = c(2),  main = 'CB with Lag, 2')
plot(cp_nvis, 'slices', lag = c(3),  main = 'CB with Lag, 3')
plot(cp_nvis, 'slices', lag = c(4),  main = 'CB with Lag, 4')

par(mfrow = c(1,1))
plot(cp_nvis)

# Mod v3 - use single flood and adjust for N_vis and in lag
# **************
mod.v3 <- gnm(n_cases ~ cb.flood.nflood_local + is_flood_week + 
                # account for covariates
                avg_tmp +
                # ********
                # and yearly time - I think just the decade, so just 2 knots?
                # can iterate here, uncertain what is best
                year,
              # ********
              data = this_df,
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')

summary(mod.v3)

# Ok this seems to work now as well
cp_incflag <- crosspred(cb.flood.nflood_local, mod.v3)

t1_incflaf <- data.frame(var = 'crosspred',
                      lb = cp_incflag$allRRlow,
                      ub = cp_incflag$allRRhigh,
                      est = cp_incflag$allRRfit)

# with lags
par(mfrow = c(2,3))
plot(cp_incflag, 'slices', lag = c(0),  main = 'CB with Lag, 0')
plot(cp_incflag, 'slices', lag = c(1),  main = 'CB with Lag, 1')
plot(cp_incflag, 'slices', lag = c(2),  main = 'CB with Lag, 2')
plot(cp_incflag, 'slices', lag = c(3),  main = 'CB with Lag, 3')
plot(cp_incflag, 'slices', lag = c(4),  main = 'CB with Lag, 4')

par(mfrow = c(1,1))
plot(cp_incflag)
