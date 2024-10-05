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
#  biomass <- c(0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0)
#  gr_curve <- data.frame(time = time, biomass = biomass)
#  
#  # Get initial parameters for the Logistic model
#  initial_params <- get_init_pars_logistic(gr_curve, this_n0 = 0.1, init_K = 30, init_lag = 0.5, init_gr_rate = 0.5)
#  
#  # Print the initial parameter approximations
#  initial_params

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
#  # Get initial parameters for the Baranyi model
#  initial_params <- get_init_pars_baranyi(gr_curve, this_n0 = 0.1, init_lag = 0.5, init_gr_rate = 0.5)
#  
#  # Print the initial parameter approximations
#  initial_params

