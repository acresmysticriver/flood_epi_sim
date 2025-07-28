
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
RR_2_lag_nvis <- 1.01  # the RR associated with each additional flood in lag weeks
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

# test plot
 ggplot(this_df) + 
   geom_point(aes(x = week_start, y = n_cases, col = n_recent_floods)) +
   facet_wrap(~strata)

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
cb.flood_local <- cb.flood[this_df$row_id, ]

# you have to transfer the attributes so crosspred works later
attr_transfer <- c('df', 'range', 'lag', 'argvar', 
                   'arglag', 'group', 'class')
for(att in attr_transfer) {
  attr(cb.flood_local, att) <- attr(cb.flood, att)
}
# *****

# make the model
N_YEARS = length(unique(this_df$year))


# Run the cross-pred
# **************
mod.v1 <- gnm(n_cases ~ cb.flood_local + is_flood_week +
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
cp <- crosspred(cb.flood_local, mod.v1)

t1 <- data.frame(var = 'crosspred',
                 lb = cp$allRRlow,
                 ub = cp$allRRhigh,
                 est = cp$allRRfit)

# with lags
par(mfrow = c(3, 2))
plot(cp, 'slices', var = c(1),  main = 'CB with Lag, 1')
plot(cp, 'slices', var = c(2),  main = 'CB with Lag, 2')
plot(cp, 'slices', var = c(3),  main = 'CB with Lag, 3')
plot(cp, 'slices', var = c(4),  main = 'CB with Lag, 4')
plot(cp, 'slices', var = c(5),  main = 'CB with Lag, 5')

plot(cp, 'overall')

# dev.off()




