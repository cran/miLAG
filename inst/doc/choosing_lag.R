## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  # Load required libraries
#  library(dplyr)
#  
#  # Generate example growth curve data
#  set.seed(123)
#  time <- 1:10
#  LOG10N <- c(2.0, 2.8, 3.6, 5.0, 7.0, 9.0, 12.0, 15.8, 20.0, 25.5)
#  gr_curve <- data.frame(t = time, LOG10N = LOG10N)
#  
#  # Fit the Baranyi model using the best algorithm
#  best_fit <- choose_lag_fit_algorithm_baranyi(gr_curve, LOG10N0 = 2.0, init_lag = 0.5, init_mumax = 0.3, init_LOG10Nmax = 30, max_iter = 100, lower_bound = 0)
#  
#  # Print the results
#  best_fit

## ----eval = FALSE-------------------------------------------------------------
#  # Load required libraries
#  library(dplyr)
#  
#  # Generate example growth curve data
#  set.seed(123)
#  time <- 1:10
#  biomass <- c(0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0)
#  gr_curve <- data.frame(time = time, biomass = biomass)
#  
#  # Fit the Logistic model using the best algorithm
#  best_fit <- choose_lag_fit_algorithm_logistic(gr_curve, n0 = 0.1, init_gr_rate = 0.5, init_K = 30, init_lag = 0.5, max_iter = 100, lower_bound = c(0, 0, 0))
#  
#  # Print the results
#  best_fit

