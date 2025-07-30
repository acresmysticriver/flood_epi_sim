# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# CREATE STRATA CROSSBASIS
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# Create crossbasis
# - important to add a GROUP variable for County
# - ORDER MATTERS !! see test code 

# *** TEST CODE ***
# cb.test <- crossbasis(df_weekly$avg_tmp[1:10], 
#                        lag=6, 
#                        argvar=list(fun = 'ns'),
#                        arglag=list(fun="poly",degree=4))
# 
# head(cb.test, 10)
# 
# cb.test <- crossbasis(df_weekly$avg_tmp[c(1:8, 10, 9)], 
#                        lag=6, 
#                        argvar=list(fun = 'ns'),
#                        arglag=list(fun="poly",degree=4))
#
# head(cb.test, 10)
# ******

# arrange
df_weekly <- df_weekly %>% arrange(county, week_start)

cb.flood.nflood <- crossbasis(df_weekly$n_recent_floods,
                       lag = 4,
                       argvar=list(fun="ns"),
                       arglag = list(fun = "poly", degree = 4),
                       group = df_weekly$county)


cb.flood.strat <- crossbasis(df_weekly$is_flood_week,
                       lag = 4,
                       argvar=list(fun="strata", breaks = (0.99)),
                       arglag = list(fun = "poly", degree = 4),
                       group = df_weekly$county)


head(cb.flood.nflood)
cb.flood.nflood
tail(cb.flood.nflood)

head(cb.flood.strat)
cb.flood.strat
tail(cb.flood.strat)

