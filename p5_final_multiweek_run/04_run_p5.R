
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
# 
# First, clear the R environment as we will reset the inputs to the scripts
#
rm(list = ls()); gc()

# Set sample from which N floods per year will be set
#
samp_set <- c(0, 1, 2)

# Set likelihood of back-to-back flood weeks
#
rep_prob <- 0.4

# Set baseline and variance
# NOTE: script 1_p4 has modified function calls to take in these inputs
#
bl <- 100000
var <- 500

# Set flood effects and lag effects for each county
# NOTE: script 1_p4 has modified function calls to take in these inputs
#
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

# Source scripts.
# Note: These will be updated based on the above!
# Note: Be sure that scripts 1 through 3 do not include the sourcing lines from
# p1 to p3! This will overwrite the inputs above.
#
source('flood_epi_sim/p5_final_multiweek_run/01_create_dummy_data_p5.R')
source('flood_epi_sim/p5_final_multiweek_run/02_create_sliding_windows_p5.R')
source('flood_epi_sim/p5_final_multiweek_run/03_strata_crossbasis_p5.R')

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# STAGE 1: CHECK OUTPUTS FROM SCRIPT 1-3, LOCALIZE CROSS BASIS
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# For the dlnm structure, effects are modeled separately by group, in this case
# by county. We will extract the number of counties and split our input data
# accordingly
#
COUNTIES <- unique(expanded_df$county)
N_COUNTIES <- length(COUNTIES)

# Split data base on county
#
county_l <- split(expanded_df, f = expanded_df$county)

# Specify which county to process. For this example, we will look at a single
# county. More on multi-county analyses to come!
#
county_i = 1

# Get subset based on county
#
this_df <- county_l[[county_i]]
dim(this_df)

# Confirm that we are only keeping strata that have a flood week. We will
# be analyzing as a conditional quasi-Poisson and therefore will not assess
# all weeks, only those that are cases and their corresponding controls
#
this_df_agg <- this_df %>%
  group_by(strata) %>%
  summarize(
    .groups = 'keep',
    n_weeks_in_strata = sum(is_flood_week)
  ) %>%
  mutate(keep = 1)

# Add these characteristics to the larger data
#
this_df <- left_join(this_df, this_df_agg,
                     by = join_by(strata))
dim(this_df)

# Assess flood effect using ggplot2
ggplot(this_df) +
   geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
              position = position_dodge(width = 1)) +
   facet_wrap(~strata)

# Let's see what the effect of consecutive floods looks like
# Subset to single events and those with event and event across lag
#
this_df <- this_df %>%
  group_by(strata) %>%
  mutate(max_n_rec = max(n_recent_floods))

# Assess flood effect using ggplot2, assessing impact o
ggplot(this_df) +
  geom_point(aes(x = week_start, y = n_cases, col = is_flood_week, group = strata),
             position = position_dodge(width = 1)) +
  facet_grid(cols = vars(max_n_rec))

# Note - part of what we will test here is the effect of a time trend. We 
# can adjust for time in multiple ways. In this repo we will control by including
# a natural spline for time
#
# To do this, we want a continuous indicator of the week from start. The below
# code can derive this variable
#
START_YEAR = 1987
this_df$week_iter = this_df$week_num + 52*(this_df$year - START_YEAR)

# We built our cross basis based on our full dataset. Because we run in a single
# county, we need to derive a subset. To do this we can use the rowID in 
# our now subset data. We want to do this for all of the cross bases we set
# up in script 3
#
cb.flood.strat_local <- cb.flood.strat[this_df$row_id, ]
cb.flood.nflood_local <- cb.flood.nflood[this_df$row_id, ]
cb.flood.flood1_local <- cb.flood.flood1[this_df$row_id, ]
cb.flood.flood2_local <- cb.flood.flood2[this_df$row_id, ]
cb.flood.flood3_local <- cb.flood.flood3[this_df$row_id, ]

# Within our local cb, we need to transfer the relevant attributes from the 
# larger crossbasis
#
attr_transfer <- c('df', 'range', 'lag', 'argvar', 
                   'arglag', 'group', 'class')
for(att in attr_transfer) {
  # Crossbasis for flood_week indicator (0, 1)
  attr(cb.flood.strat_local, att) <- attr(cb.flood.strat, att)
  
  # Cross basis for consecutive flood weeks
  attr(cb.flood.nflood_local, att) <- attr(cb.flood.nflood, att)
  
  # Cross basis for exp1 indicator (only first week of flooding)
  attr(cb.flood.flood1_local, att) <- attr(cb.flood.flood1, att)
  
  # Cross basis for exp2 indicator (only second week of flooding)
  attr(cb.flood.flood2_local, att) <- attr(cb.flood.flood2, att)
  
  # Cross basis for exp3 indicator (only third week of flooding)
  attr(cb.flood.flood3_local, att) <- attr(cb.flood.flood3, att)
}

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# STAGE 2: RUN MODELS!
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# To assess whether we have a model set up to effectively estimate our provided
# flood and lag effect. We are going to run multiple models. 
#
#       MOD V0 - First, we will run a model that is not using the crossbasis
#                or dlnm, but instead provide the lag and flood explicitly,
#                along with a year effect and a covariate with no associated
#                effect (temperature) in this analysis. We'll first check
#                what the result looks like without adjusting for the added
#                impact of consecutive flood weeks we have injected in the data

############### MOD V0 - INCLUDE HARD CODE FLOOD AND N RECENT ##################
# Mod v0 - use no cross basis. Include N recent floods and lag term 
# for N recent floods
# 
mod.v0_unadj <- gnm(n_cases ~ is_flood_week + 
                lag1 + lag2 + lag3 + lag4 +
               avg_tmp +
               year,
             data = this_df,
             family = quasipoisson,
             eliminate = factor(strata),
             subset = keep == 1,
             na.action = 'na.exclude')

# Assess whether model appears to have run effectively. Note that within the 
# gnm application, we are looking for exponentiated outputs, so this summary
# can be non-intuitive to interpret.
#
summary(mod.v0_unadj)

# To get a better sense of our effect size, let's exponentiate!
# Look at table outputs. We are running for county 1. We expect
# that we will get back exactly what we put in for flood and lag effects
# Check these inputs
#
RR_1
RR_1_lag

# Exponentiate the confidence intervals from effect
#
t0_unadj <- exp(confint(mod.v0_unadj))
t0_unadj <- data.frame(t0_unadj, var = row.names(t0_unadj), est = exp(coef(mod.v0_unadj)))
colnames(t0_unadj)[1:2] <- c('lb', 'ub')
t0_unadj

# Without adjusting, the RR and lag are all inflated! Now let's assess the 
# effect when accounting for n_recent floods
#
# Mod v0 - use no cross basis. Include N recent floods and lag term 
# for N recent floods
# 
mod.v0 <- gnm(n_cases ~ is_flood_week + 
                      lag1 + lag2 + lag3 + lag4 +
                      n_recent_floods +
                      lag_nflood_1 + lag_nflood_2 + lag_nflood_3 + lag_nflood_4 +
                      avg_tmp +
                      year,
                    data = this_df,
                    family = quasipoisson,
                    eliminate = factor(strata),
                    subset = keep == 1,
                    na.action = 'na.exclude')

# Assess whether model appears to have run effectively. Note that within the 
# gnm application, we are looking for exponentiated outputs, so this summary
# can be non-intuitive to interpret.
#
summary(mod.v0)

# To get a better sense of our effect size, let's exponentiate!
# Look at table outputs. We are running for county 1. We expect
# that we will get back exactly what we put in for flood and lag effects
# Check these inputs
#
RR_1
RR_1_lag
RR_1_nflood
RR_1_lag_nflood

# Exponentiate the confidence intervals from effect
#
t0 <- exp(confint(mod.v0))
t0 <- data.frame(t0, var = row.names(t0), est = exp(coef(mod.v0)))
colnames(t0)[1:2] <- c('lb', 'ub')
t0

# Do the output estimates match the inputs? If so, it seems that the RRs in 
# script 1 were applied effectively and this model can effectively assess
# the risk of flooding in this context.
# Now, let's assess how use of the dlnm compares.

############### MOD V1 - INCLUDE CROSS BASIS OF FLOODING AND N_Recent ##########
# Mod v1 - use flood crossbasis
#
mod.v1 <- gnm(n_cases ~ cb.flood.strat_local + 
                cb.flood.nflood_local +
                avg_tmp +
                year,
              data = this_df,
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')

summary(mod.v1)

# Exponentiate coefficients
#
exp(coef(mod.v1))

# To visualize the effect, we can use a crossprediction. Note, this has to be 
# done separately for each cross basis
#
#
cp.strata <- crosspred(cb.flood.strat_local, mod.v1)
cp.nvis.strata <- crosspred(cb.flood.nflood_local, mod.v1)

t1_isflood <- data.frame(var = 'crosspred',
                 lb = cp.strata$allRRlow,
                 ub = cp.strata$allRRhigh,
                 est = cp.strata$allRRfit)

t1_nvis <- data.frame(var = 'crosspred',
                      lb = cp.nvis.strata$allRRlow,
                      ub = cp.nvis.strata$allRRhigh,
                      est = cp.nvis.strata$allRRfit)

# Look at plot outputs. 
# The below will generate a left- and right-hand plot. 
# On the left hand, the effect is shown at weeks 0 through 4, with week 0 
# representing the effect the week of flooding and weeks 1 through 4 
# representing the lagged effect for each week. On the right side, we will 
# see a plot of the cumulative risk. Because we have set up the flood exposure
# as a strata with breaks at 0.99, this will be a flat line at 1, up until the 
# x axis reaches one, where the peak represents the cumulative effect across
# the lag period. For inputs of RR 2 and RR_lag 1.15, this results in a 
# cumulative effect around 3.5 (2 * 1.15 * 1.15 * 1.15 * 1.15)
#
par(mfrow = c(2, 2))
plot(cp.strata, 'slices', var = c(1),  main = 'CB for Is Flood Week == 1')
plot(cp.strata, 'overall')

plot(cp.nvis.strata, 'slices', var = c(1),  main = 'CB For N Recent Floods')
plot(cp.nvis.strata, 'overall')

# It looks like the dlnm was able to effectively address these inputs!
# To ensure our model is robust, we can modify the inputs at top and ensure 
# they carry throguh the model outputs.
#
#       How does the addition of a reducing or increasing time effect impact
#       our results?
#       What about increasing the likelihood of repeat flood events?
#       Increasing the number of flood events?
#       Removing the lag effect?
#
# While the cross basis appears to be effective, we might want to also test
# in our 'real' data in the future the effect of different flood effects. Let's
# test the impact of differing categorical floods
#
############### MOD V2 - INCLUDE CATEGORICAL CROSS BASIS #######################
# Mod v2 - use categorical flood crossbases
# **************
mod.v1 <- gnm(n_cases ~ cb.flood.flood1_local + cb.flood.flood2_local + 
                cb.flood.flood3_local +
                avg_tmp +
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

# with lags
par(mfrow = c(3, 2))
plot(cp_1.strata, 'slices', var = c(1),  main = 'CB with Lag, One Week Flood')
plot(cp_1.strata, 'overall')

plot(cp_2.strata, 'slices', var = c(1),  main = 'CB with Lag, Second Week Flood')
plot(cp_2.strata, 'overall')

plot(cp_3.strata, 'slices', var = c(1),  main = 'CB with Lag, Third Week Flood')
plot(cp_3.strata, 'overall')

# We cannot extract the impacts in the same way as with the n_recent_flood 
# cross basis, but this approach allows for quantifying the different impacts
# of floods of differing durations

