## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Load required libraries
# library(dplyr)
# 
# # Generate example growth curve data
# set.seed(123)
# time <- 1:10
# biomass <- c(0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0)
# gr_curve <- data.frame(time = time, biomass = biomass)
# 
# # Calculate lag for the Logistic model
# lag_result <- calc_lag_fit_to_logistic_with_lag(gr_curve, n0 = 0.1, init_gr_rate = 0.5, init_K = 30, init_lag = 0.5)
# 
# # Print the calculated lag
# lag_result$lag_N

## ----eval= FALSE--------------------------------------------------------------
# # Load required libraries
# library(dplyr)
# 
# # Generate example growth curve data
# set.seed(123)
# time <- 1:10
# biomass <- c(0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0)
# gr_curve <- data.frame(time = time, biomass = biomass)
# 
# # Calculate lag for the Baranyi model
# lag_result <- calc_lag_fit_to_baranyi_with_lag(gr_curve, LOG10N0 = NULL, init_lag = NULL, init_mumax = NULL, init_LOG10Nmax = NULL, algorithm = "auto")
# 
# # Print the calculated lag
# lag_result$lag_N

