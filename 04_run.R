
# SET SWITCHES FOR TESTING
# -- a 4 week lag (so 0, 1, 2, 3, 4) is hard-coded

# clear
rm(list = ls()); gc()

# COUNTY A
RR_1     <- 2.10  # the RR on lag 0
RR_1_lag <- 1.20  # the RR on lags 1:4
RR_1_double <- 2.10  # the RR associated with flood on day and in lag
ybeta_1  <- 1.00  # the year trend (in log space, so 2 means doubling every year)
bl <- 10000       # n case baseline 
var <- 500        # variance in case by week

# COUNTY B
RR_2     <- 2.00  # the RR on lag 0
RR_2_lag <- 1.00 # the RR on lags 1:4
RR_2_double <- 2.30  # the RR associated with flood on day and in lag
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
   geom_point(aes(x = week_start, y = n_cases, col = double_flood)) +
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

# **************
# OK SO FIRST, JUST DO A SIMPLE MODEL TO MAKE SURE YOU KNOW
# WHAT IS GOING ON
mod.v0 <- gnm(n_cases ~ is_flood_week + 
                # manually do lags
                lag1 + lag2 + lag3 + lag4 +
                # account for covariates
                avg_tmp +
                # ********
                # and yearly time - I think just the decade, so just 2 knots?
                # can iterate here, uncertain what is best
                year +
                ## add double flood marker
                double_flood,
              # ********
              data = this_df,
              family = quasipoisson,
              eliminate = factor(strata),
              subset = keep == 1,
              na.action = 'na.exclude')

summary(mod.v0)

# AWESOME, this works
# you see that the is_case_period coefficent and year reflects
# what you set it as 

t1 <- exp(confint(mod.v0))
t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(mod.v0)))
colnames(t1)[1:2] <- c('lb', 'ub')
t1
# Example data structure assumed for conversion
# t1 <- data.frame(var = c("A", "B", "C"), lb = c(1, 2, 3), est = c(2, 3, 4), ub = c(3, 4, 5))

# Reverse the order of y-axis categories (to mimic ggplot behavior)
# Set up blank plot with appropriate limits
par(mar = c(5, 9, 4, 2))
y_vals <-1:nrow(t1)
plot(NULL,
     xlim = range(c(t1$lb, t1$ub)),
     ylim = range(1:nrow(t1)),
     yaxt = "n",
     xlab = "Estimate",
     ylab = "",
     main = "Base model")

# Add y-axis labels
axis(2, at = y_vals, labels = t1$var, las = 1)

# Add segments (error bars)
segments(x0 = t1$lb, x1 = t1$ub, y0 = y_vals, y1 = y_vals)

# Add points (estimates)
points(t1$est, y_vals, pch = 5, cex = 0.8)


# **************
# Next, try and make the year variable using ns()
mod.v0b <- gnm(n_cases ~ is_flood_week +
                 # manually do lags
                 lag1 + lag2 + lag3 + lag4 +
                 # is it double
                double_flood +
                 # account for covariates
                 avg_tmp +
                 # ********
                 # and yearly time - I think just the decade, so just 2 knots?
                 # can iterate here, uncertain what is best
                 ns(week_iter, df = 4),
               # ********
               data = this_df,
               family = quasipoisson,
               eliminate = factor(strata),
               subset = keep == 1,
               na.action = 'na.exclude')

summary(mod.v0b)

# seems to work
# you see that the is_case_period coefficent reflects
# what you set it as
exp(coef(mod.v0b))
exp(confint(mod.v0b))
# *************

t1 <- exp(confint(mod.v0b))
t1 <- data.frame(t1, var = row.names(t1), est = exp(coef(mod.v0b)))
colnames(t1)[1:2] <- c('lb', 'ub')
t1
y_vals <-1:nrow(t1)
par(mar = c(5, 9, 4, 2))
plot(NULL,
     xlim = range(c(t1$lb, t1$ub)),
     ylim = range(1:nrow(t1)),
     yaxt = "n",
     xlab = "Estimate",
     ylab = "",
     main = "Base model with ns(year)")

# Add y-axis labels
axis(2, at = y_vals, labels = t1$var, las = 1)

# Add segments (error bars)
segments(x0 = t1$lb, x1 = t1$ub, y0 = y_vals, y1 = y_vals)

# Add points (estimates)
points(t1$est, y_vals, pch = 5, cex = 0.8)


# Ok now try the cross-pred
# **************
mod.v1 <- gnm(n_cases ~ cb.flood_local +
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
cp <- crosspred(cb.flood_local, mod.v1, cen = 0, at = 1)

t1 <- data.frame(var = 'crosspred',
                 lb = cp$allRRlow,
                 ub = cp$allRRhigh,
                 est = cp$allRRfit)

t1
y_vals <-1
par(mar = c(5, 9, 4, 2))
plot(NULL,
     xlim = range(c(t1$lb, t1$ub)),
     ylim = range(1:nrow(t1)),
     yaxt = "n",
     xlab = "Estimate",
     ylab = "",
     main = "Overall with CB")

# Add y-axis labels
axis(2, at = y_vals, labels = t1$var, las = 1)

# Add segments (error bars)
segments(x0 = t1$lb, x1 = t1$ub, y0 = y_vals, y1 = y_vals)

# Add points (estimates)
points(t1$est, y_vals, pch = 5, cex = 0.8)

# with lags
plot(cp, 'slices', var = c(1),  main = 'CB with Lag')


# dev.off()




