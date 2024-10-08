% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milags_functions.R
\name{plot_lag_fit}
\alias{plot_lag_fit}
\title{plot_lag_fit}
\usage{
plot_lag_fit(data_new, print_lag_info = TRUE, log10_transform = TRUE)
}
\arguments{
\item{data_new}{a data frame output by Calculate.Lag function: it needs to have the following columns: "time", "biomass", "tangent.point", "predicted.data", "threshold", "N0", "second.deriv.b", "line.intercept", "line.slope"}

\item{print_lag_info}{if set to "TRUE" prints the lag length on the graph}

\item{log10_transform}{if to plot y axis (biomass) on log10 scale}
}
\value{
ggplot object with a growth curve
}
\description{
Plots the provided growth curve (one single growth curve) together with the calculated lag and and the rationale for lag calculation
}
