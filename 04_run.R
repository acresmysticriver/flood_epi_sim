
# SET SWITCHES FOR TESTING
# -- a 4 week lag (so 0, 1, 2, 3, 4) is hard-coded

# clear
rm(list = ls()); gc()

# set whether controls for windows should have no flood in the window (noflood),
# no lag in the window (nolag), or no restricitons on flood/lag (norest) 
control_select <- "nolag"

# set whether the effect of a multi-week flood

# set series of probabilities for flood in a year
samp_set <- c(1, 2, 6, 7, 7, 10)

# set likelihood of back-to-back flood weeks
rep_prob <- 0.95

# set how many weeks to use in lag
lag_weeks <- 4

# COUNTY A
RR_1     <- 2.00  # the RR on lag 0
RR_1_lag <- 1.15  # the RR on lags 1:4
RR_1_nvis <- 1.035  # the RR associated with each additional flood in recent weeks
RR_1_lag_nvis <- 1.02  # the RR associated with each additional flood in lag weeks
ybeta_1  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl_1 <- 1000000      # n case baseline 
var_1 <- 500       # variance in case by week

# COUNTY B
RR_2     <- 1.60  # the RR on lag 0
RR_2_lag <- 1.08 # the RR on lags 1:4
RR_2_nvis <- 1.20  # the RR associated with each additional flood in recent weeks
RR_2_lag_nvis <- 1.03  # the RR associated with each additional flood in lag weeks
ybeta_2  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl_2 <- 10000      # n case baseline 
var_2 <- 500       # variance in case by week

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

# See what the minimum case is within each case strata
#
check_mincase <- this_df %>%
  group_by(strata) %>% 
  filter(period_label == "case") %>%
  summarise(min_cases = min(n_cases))

# # test plot
 ggplot(this_df# [this_df$strata %in% c("CountyA_Flood01", "CountyA_Flood02") & this_df$year == 1987,]
        ) +
   geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
              position = position_dodge(width = 1)) +
   facet_grid(cols = vars(max_n_rec))
 # # 
 # ggplot(this_df[this_df$max_n_rec == 3 & this_df$strata == "CountyA_Flood19", ]
 # ) +
 #   geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
 #              position = position_dodge(width = 1)) +
 #   facet_grid(cols = vars(max_n_rec))
 # # 
 # ggplot(this_df[this_df$strata == "CountyA_Flood17" ,]
 # ) +
 #   geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
 #              position = position_dodge(width = 1)) +
 #   facet_grid(cols = vars(year), scales = "free_x")

#

# re-estimate case to make sure it matches expectation
 # this_df$n_case_exp <- this_df$case_baseline + (this_df$case_baseline * (this_df$RR - 1) * this_df$is_flood_week)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag - 1) * this_df$lag1)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag - 1) * this_df$lag2)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag - 1) * this_df$lag3)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag - 1) * this_df$lag4)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_nvis - 1) * this_df$n_recent_floods)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag_nvis - 1) * this_df$lag_nvis_1)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag_nvis - 1) * this_df$lag_nvis_2)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag_nvis - 1) * this_df$lag_nvis_3)
 # this_df$n_case_exp <- this_df$n_case_exp + (this_df$n_case_exp * (this_df$RR_lag_nvis - 1) * this_df$lag_nvis_4)
 # 
 # this_df$case_diff <-  this_df$n_case_exp - this_df$n_cases
 # 
 # this_df$prop <- this_df$case_diff / this_df$n_cases

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
cb.flood.flood1_local <- cb.flood.flood1[this_df$row_id, ]
cb.flood.flood2_local <- cb.flood.flood2[this_df$row_id, ]
cb.flood.flood3_local <- cb.flood.flood3[this_df$row_id, ]

# you have to transfer the attributes so crosspred works later
attr_transfer <- c('df', 'range', 'lag', 'argvar', 
                   'arglag', 'group', 'class')
for(att in attr_transfer) {
  # This crossbasis is based on the cumulative flood indicator
  # (ie: second week of flood = 2, third = 3...)
  attr(cb.flood.nflood_local, att) <- attr(cb.flood.nflood, att)
  
  # Crossbasis for flood_week indicator (0, 1)
  attr(cb.flood.strata_local, att) <- attr(cb.flood.strat, att)
  
  # Cross basis for exp1 indicator (only first week of flooding)
  attr(cb.flood.flood1_local, att) <- attr(cb.flood.flood1, att)
  
  # Cross basis for exp2 indicator (only second week of flooding)
  attr(cb.flood.flood2_local, att) <- attr(cb.flood.flood2, att)
  
  # Cross basis for exp3 indicator (only third week of flooding)
  attr(cb.flood.flood3_local, att) <- attr(cb.flood.flood3, att)
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

############### MOD V0 - INCLUDE HARD CODE FLOOD AND N RECENT ##################
# Mod v0 - use no cross basis. Include N recent floods and lag term 
# for N recent floods
# **************
mod.v0 <- gnm(n_cases ~ is_flood_week + n_recent_floods + 
               # manually do lags
                lag1 + lag2 + lag3 + lag4 +
                lag_nvis_1  + lag_nvis_2 + lag_nvis_3  + lag_nvis_4  + 
                
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

summary(mod.v0)

# Look at table outputs. For county 2:
#     is_flood_week appears accurate.
#     n_recent_floods 1.09 v 1.10 (small underestimat)
#     lag and lag_nvis both underestimates also.
t0 <- exp(confint(mod.v0))
t0 <- data.frame(t0, var = row.names(t0), est = exp(coef(mod.v0)))
colnames(t0)[1:2] <- c('lb', 'ub')
t0

############### MOD V1 - INCLUDE CROSS BASIS OF FLOODING AND N_Recent ##########
# Mod v1 - use flood alone as cross basis
# **************
mod.v1 <- gnm(n_cases ~ cb.flood.strata_local + cb.flood.nflood_local +
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

exp(coef(mod.v1))

# Ok this seems to work now as well
cp.strata <- crosspred(cb.flood.strata_local, mod.v1)
cp.nvis.strata <- crosspred(cb.flood.nflood_local, mod.v1)

t1_isflood <- data.frame(var = 'crosspred',
                 lb = cp.strata$allRRlow,
                 ub = cp.strata$allRRhigh,
                 est = cp.strata$allRRfit)

t1_nvis <- data.frame(var = 'crosspred',
                 lb = cp.nvis.strata$allRRlow,
                 ub = cp.nvis.strata$allRRhigh,
                 est = cp.nvis.strata$allRRfit)


# Look at plot outputs. For county 2:
#     is_flood_week appears accurate.
#     n_recent_floods 1.09 v 1.10 (small underestimat)
#     lag and lag_nvis both underestimates also.

# with lags
par(mfrow = c(2, 2))
plot(cp.strata, 'slices', var = c(1),  main = 'CB for Is Flood Week == 1')
plot(cp.strata, 'overall')

plot(cp.nvis.strata, 'slices', var = c(1),  main = 'CB For N Recent Floods == 2')
plot(cp.nvis.strata, 'overall')
# dev.off()

############### MOD V2 - INCLUDE CATEGORICAL CROSS BASIS #######################
# Mod v2 - use categorical flood crossbases
# **************
mod.v1 <- gnm(n_cases ~ cb.flood.flood1_local + cb.flood.flood2_local + 
                cb.flood.flood3_local +
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
cp_1.strata <- crosspred(cb.flood.flood1_local, mod.v1)
cp_2.strata <- crosspred(cb.flood.flood2_local, mod.v1)
cp_3.strata <- crosspred(cb.flood.flood3_local, mod.v1)

t_flood1 <- data.frame(var = 'crosspred',
                 lb = cp_1.strata$allRRlow,
                 ub = cp_1.strata$allRRhigh,
                 est = cp_1.strata$allRRfit)


t_flood2 <- data.frame(var = 'crosspred',
                       lb = cp_2.strata$allRRlow,
                       ub = cp_2.strata$allRRhigh,
                       est = cp_2.strata$allRRfit)

t_flood3 <- data.frame(var = 'crosspred',
                       lb = cp_3.strata$allRRlow,
                       ub = cp_3.strata$allRRhigh,
                       est = cp_3.strata$allRRfit)

t_flood1
t_flood2
t_flood3

# with lags
par(mfrow = c(3, 2))
plot(cp_1.strata, 'slices', var = c(1),  main = 'CB with Lag, One Week Flood')
plot(cp_1.strata, 'overall')

plot(cp_2.strata, 'slices', var = c(1),  main = 'CB with Lag, Second Week Flood')
plot(cp_2.strata, 'overall')

plot(cp_3.strata, 'slices', var = c(1),  main = 'CB with Lag, Third Week Flood')
plot(cp_3.strata, 'overall')


############### MOD V3 - RUN ON STRATIFIED DATA ##################
# Mod v3 - use no cross basis. Include N recent floods and lag term 
# for N recent floods
# **************
mod.v3_sing <- gnm(n_cases ~ is_flood_week + 
                # manually do lags
                lag1 + lag2 + lag3 + lag4 +
                # account for covariates
                avg_tmp +
                # ********
                # and yearly time - I think just the decade, so just 2 knots?
                # can iterate here, uncertain what is best
                year,
              # ********
              data = this_df[this_df$max_n_rec == 0, ],
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')

mod.v3_two <- gnm(n_cases ~ is_flood_week + n_recent_floods + 
                     # manually do lags
                     lag1 + lag2 + lag3 + lag4 +
                     lag_nvis_1  + lag_nvis_2 + lag_nvis_3  + lag_nvis_4  + 
                     # account for covariates
                     avg_tmp +
                     # ********
                     # and yearly time - I think just the decade, so just 2 knots?
                     # can iterate here, uncertain what is best
                     year,
                   # ********
                   data = this_df[this_df$max_n_rec == 1, ],
                   family = quasipoisson,
                   eliminate = factor(strata),
                   subset = keep == 1,
                   na.action = 'na.exclude')

mod.v3_three <- gnm(n_cases ~ is_flood_week + n_recent_floods + 
                    # manually do lags
                    lag1 + lag2 + lag3 + lag4 +
                    lag_nvis_1  + lag_nvis_2 + lag_nvis_3  + lag_nvis_4  + 
                    # account for covariates
                    avg_tmp +
                    # ********
                    # and yearly time - I think just the decade, so just 2 knots?
                    # can iterate here, uncertain what is best
                    year,
                  # ********
                  data = this_df[this_df$max_n_rec == 2, ],
                  family = quasipoisson,
                  eliminate = factor(strata),
                  subset = keep == 1,
                  na.action = 'na.exclude')

# Look at table outputs. For county 2:
#     is_flood_week appears accurate.
#     n_recent_floods 1.09 v 1.10 (small underestimat)
#     lag and lag_nvis both underestimates also.
t0_1 <- exp(confint(mod.v3_sing))
t0_1 <- data.frame(t0_1, var = row.names(t0_1), est = exp(coef(mod.v3_sing)))
colnames(t0_1)[1:2] <- c('lb', 'ub')
t0_1

t0_2 <- exp(confint(mod.v3_two))
t0_2 <- data.frame(t0_2, var = row.names(t0_2), est = exp(coef(mod.v3_two)))
colnames(t0_2)[1:2] <- c('lb', 'ub')
t0_2

t0_3 <- exp(confint(mod.v3_three))
t0_3 <- data.frame(t0_3, var = row.names(t0_3), est = exp(coef(mod.v3_three)))
colnames(t0_3)[1:2] <- c('lb', 'ub')
t0_3

