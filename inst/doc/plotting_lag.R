## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# # Load required libraries
# library(dplyr)
# 
# # Generate example data using Calculate.Lag function
# set.seed(123)
# time <- 1:10
# biomass <- c(0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0)
# tangent.point <- c(0.3, 0.5, 0.9, 2.0, 4.0, 6.0, 9.0, 12.0, 17.0, 24.0)
# predicted.data <- c(0.1, 0.4, 0.8, 1.6, 3.2, 6.0, 8.8, 12.5, 18.2, 25.1)
# threshold <- c(0.3, 0.8, 1.3, 2.3, 4.3, 7.0, 10.2, 15.0, 21.0, 28.0)
# N0 <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# second.deriv.b <- c(0.02, 0.04, 0.09, 0.2, 0.4, 0.6, 0.9, 1.2, 1.7, 2.4)
# line.intercept <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# line.slope <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
# 
# data_new <- data.frame(
#   time = time,
#   biomass = biomass,
#   tangent.point = tangent.point,
#   predicted.data = predicted.data,
#   threshold = threshold,
#   N0 = N0,
#   second.deriv.b = second.deriv.b,
#   line.intercept = line.intercept,
#   line.slope = line.slope
# )
# 
# # Plot the growth curve with lag information
# plot <- plot_lag_fit(data_new, print_lag_info = TRUE, log10_transform = TRUE)
# 
# # Print the plot
# print(plot)

