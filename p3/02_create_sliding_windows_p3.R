# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# Code for Simulating Weekly Flooding Epi Analysis
#
# We will create conditional Quasi-Poisson selecting on county and week. 
# We use two-stage selection with sliding strata on years before and after a 
# flood event, following the proces defined as in the Aggarwal et al. 
# pre-print: https://arxiv.org/abs/2309.13142
#
# P3:
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
# 
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# /////////////////////////////////////////////////////////////////////////////
# CREATE SLIDING WINDOWS
# /////////////////////////////////////////////////////////////////////////////
# -----------------------------------------------------------------------------

# To run a condiitonal quasi-Poission analysis, we will need to select matched
# controls for our flood 'case' weeks as indicated in script 01. This script
# is meant to run after script 01 and will require the final output from that
# script (df_weekly) to be run. You can run script 1 or use the below to bring
# in the data:
#
source('flood_epi_sim/p3/01_create_dummy_data_p3.R')

# For our case sample, we have a flood and four lag weeks.
# For our control sample, because of the potential for longer-lasting flood 
# effects, we are going to choose weeks from preceding and following years.
# For the second week in 2003, we will have lag weeks 3, 4, 5, 6.
# If there was a flood in week 2, 2003, then we would want to select two 
# controls - one from 2002 (week 2) and one from 2004 (week 2). We will also
# select the weeks 3, 4, 5, 6 from these years to cross reference the lag effect.
# If there was a flood during this timeframe in 2002, then we will look further
# back to 2001, 2000, etc. until we find a 'clean' control option. We will do 
# the same when looking forward for a control week.


# 1) To start, we will give each flood a number, this unique identifier can be
# used to explore specific events, and to subset to a specific event in order
# to identify its control weeks
#
df_weekly$flood_name <- NA
flood_num = 1

# Loop through all rows. If there is a flood, add the flood name (iterating to
# increase the flood index number as we go)
#
for(i in 1:nrow(df_weekly)) {
  if(df_weekly$is_flood_week[i] == 1) {
    df_weekly$flood_name[i] <- paste0('Flood_', flood_num)
    flood_num = flood_num + 1
  }
}

# We can now see how many strata we have based on our flood_num indicator
#
N_STRATA = flood_num - 1

# State how many weeks we will be exploring in the lag period
#
N_WEEKS_PER_PERIOD <- 5  # (main week + 4 weeks)

# 2) Next we will make our strata. Because we are working with many strata, and
# will seek to find controls for each, we are going to do this as a function. 
# When you are first setting up your analysis, work within the function to 
# explore the process. (Use the stata_i = 5 test input to explore)
#
make_strata <- function(strata_i) {
  
  # ****
  # strata_i = 5
  # ****
  
  cat("Flood_", strata_i, "\t")
  
  # As we pass strata to this function, we will subset our data to the event
  # of interest. The result will be the 'this_flood' dataframe
  this_flood <- df_weekly %>%
    filter(flood_name == paste0("Flood_", strata_i))
  
  # Extract the temporal components of the event
  this_year <- this_flood$year
  this_month <- month(this_flood$week_start)
  this_day <- day(this_flood$week_start)
  
  # Define the current year's exposure period. We have set the N weeks above,
  # and to find the end of the period we will multiply the n weeks by 7 days
  # per week
  end_of_this_period = this_flood$week_start + N_WEEKS_PER_PERIOD * 7
  
  # Now, restricting on these dates, we can find our case period, including 
  # both the flood event and the resulting lag period
  case_data <- df_weekly %>%
    filter(week_start >= this_flood$week_start,
           week_start < end_of_this_period,
           county == this_flood$county) %>%
    mutate(is_case_period = 1,
           period_label = 'case')
  
  # Next, we will define our 'prior' control
  # Because we want to search until we find a period without a flood, we will
  # use the while loop to continue until we have found a 'clean' control
  search_for_prior = T
  year_i = 1
  while(search_for_prior) {
    
    # Make a new time period based on the previous
    # We subtract year_i because we are searching backwards
    start_of_this_period <- this_flood$week_start - year_i * 365
    end_of_this_period <- start_of_this_period + N_WEEKS_PER_PERIOD * 7
    
    # Get cut of prior data based on these parameters
    prior_data <- df_weekly %>%
      filter(week_start >= start_of_this_period,
             week_start < end_of_this_period,
             county == this_flood$county) %>%
      mutate(is_case_period = 0,
             period_label = 'prior')
    
    # Reset to F so our while loop advances
    search_for_prior = F
    
    # First, we will check if there are no previous years left to get a control
    # from. If this is the case, we have no prior period but end our search
    if(nrow(prior_data) == 0) {
      search_for_prior = F
    }
    
    # Next, if it has a flood in it, iterate year and try again
    if(any(prior_data$is_flood_week == 1) ) {
      
      # Set search indicator to true as we will want to look further back
      search_for_prior = T
      year_i = year_i - 1 # search back further
    } 
  }
  
  # Next, we can define our post-flood control. We will search forward through
  # the years to find this period
  search_for_post = T
  year_i = 1
  while(search_for_post) {
    
    # Make a new time period based on the previous
    # We add year_i because we are searching forwards
    start_of_this_period <- this_flood$week_start + year_i * 365
    end_of_this_period <- start_of_this_period + N_WEEKS_PER_PERIOD * 7
    
    # Get cut of post data based on these parameters
    post_data <- df_weekly %>%
      filter(week_start >= start_of_this_period,
             week_start < end_of_this_period,
             county == this_flood$county) %>%
      mutate(is_case_period = 0,
             period_label = 'post')
    
    # Reset to F so our while loop advances
    search_for_post = F
    
    # Check if there are no following years left to get a control
    # from. If this is the case, we have no post period but end our search
    if(nrow(post_data) == 0) {
      search_for_post = F
    }
    
    # Next, if it has a flood in it, iterate year and try again
    if(any(post_data$is_flood_week == 1) ) {
      
      # Set search indicator to true as we will want to look further back
      search_for_post = T
      year_i = year_i + 1 # search forward further
    } 
  }
  
  # For our output we want to include our case with its two control periods
  # if available. We will include with this data cut an indicator of the county,
  # and flood name
  out <- do.call(rbind, list(prior_data, case_data, post_data))
  out$strata <- paste0(unique(out$county), sprintf("_Flood%02i", strata_i))
  return(out)
  
}

# After running through the above with a single strata, we can apply the function
# to create case/control sets for each of our floods
#
STRATA_L <- future_lapply(1:N_STRATA, make_strata)

# We want to have our output include all of the case and controls for all floods,
# so we combine
#
expanded_df <- do.call(rbind, STRATA_L)

# Did this work? Let's look at a single flood
#
d <- expanded_df %>% filter(strata == 'CountyA_Flood14')
d

# Let's visualize our single event. We Can see how the timing aligns for 
# distinct years
#
ggplot(d) +
  geom_point(aes(x = week_start, y = n_cases, col = factor(is_flood_week), size = factor(lagflag))) +
  labs(title = "Time Series of Weekly Case Counts in Two Fictional Counties",
       x = "Week",
       y = "Number of Cases",
       col = "Flood Week (1 = Flooded)",
       size = "Lag Week (1= Lag)") +
  theme_minimal() +
  facet_wrap(~year, scale = "free_x")
