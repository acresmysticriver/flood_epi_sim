# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# CREATE SLIDING WINDOWS
# - This is where you'd start with your data
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# the current year, + 4 weeks
# then both the previous and following year, the same weeks
# but only if there isn't a flood
# if there is, then they iterate outwards to get the next year

# 1) Give each flood a number
df_weekly$flood_name <- NA
flood_num = 1
for(i in 1:nrow(df_weekly)) {
  if(df_weekly$is_flood_week[i] == 1) {
    df_weekly$flood_name[i] <- paste0('Flood_', flood_num)
    flood_num = flood_num + 1
  }
}

N_STRATA = flood_num - 1
N_WEEKS_PER_PERIOD <- 5  # (main week + 4 weeks)

# 2) a function to make the strata
make_strata <- function(strata_i) {
  
  # ****
  # strata_i = 5
  # ****
  
  cat("Flood_", strata_i, "\t")
  
  # find the flood
  this_flood <- df_weekly %>%
    filter(flood_name == paste0("Flood_", strata_i))
  this_year <- this_flood$year
  this_month <- month(this_flood$week_start)
  this_day <- day(this_flood$week_start)
  
  # define the current year's exposure period
  end_of_this_period = this_flood$week_start + N_WEEKS_PER_PERIOD * 7
  
  case_data <- df_weekly %>%
    filter(week_start >= this_flood$week_start,
           week_start < end_of_this_period,
           county == this_flood$county) %>%
    mutate(is_case_period = 1,
           period_label = 'case')
  
  # ----------
  # now define the prior period
  # lots of edge cases to worry about here
  search_for_prior = T
  year_i = 1
  while(search_for_prior) {
    
    # make a new time period
    # its - year_i because you are searching backwards
    start_of_this_period <- this_flood$week_start - year_i * 365
    end_of_this_period <- start_of_this_period + N_WEEKS_PER_PERIOD * 7
    
    prior_data <- df_weekly %>%
      filter(week_start >= start_of_this_period,
             week_start < end_of_this_period,
             county == this_flood$county) %>%
      mutate(is_case_period = 0,
             period_label = 'prior')
    
    search_for_prior = F
    
    # first is it runs off the end
    if(nrow(prior_data) == 0) {
      search_for_prior = F
    }
    
    # next, it has a flood in it, so iterate year and try again
    if(any(prior_data$is_flood_week == 1)) {
      search_for_prior = T
      year_i = year_i + 1 # search back further
    } 
  }
  
  # ----------
  # now define the post period
  # lots of edge cases to worry about here
  search_for_post = T
  year_i = 1
  while(search_for_post) {
    
    # make a new time period
    # this time its + year_i
    start_of_this_period <- this_flood$week_start + year_i * 365
    end_of_this_period <- start_of_this_period + N_WEEKS_PER_PERIOD * 7
    
    post_data <- df_weekly %>%
      filter(week_start >= start_of_this_period,
             week_start < end_of_this_period,
             county == this_flood$county) %>%
      mutate(is_case_period = 0,
             period_label = 'post')
    
    search_for_post = F
    
    # first is it runs off the end
    if(nrow(post_data) == 0) {
      search_for_post = F
    }
    
    # next, it has a flood in it, so iterate year and try again
    if(any(post_data$is_flood_week == 1)) {
      search_for_post = T
      year_i = year_i + 1 # + 1 because searching forwards
    } 
  }
  
  # ----- 
  # combine
  out <- do.call(rbind, list(prior_data, case_data, post_data))
  out$strata <- paste0(unique(out$county), sprintf("_Flood%02i", strata_i))
  return(out)
  
}

# doing the strata
STRATA_L <- future_lapply(1:N_STRATA, make_strata)

# get expanded df
expanded_df <- do.call(rbind, STRATA_L)

# check that it worked
# great, so you're seeing the increase and also the year-trend here
# d <- expanded_df %>% filter(strata == 'CountyA_Flood14')
# d
# 

