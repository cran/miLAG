% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milags_functions.R
\name{calc_lag_fit_to_baranyi_with_lag}
\alias{calc_lag_fit_to_baranyi_with_lag}
\title{calc_lag_fit_to_baranyi_with_lag}
\usage{
calc_lag_fit_to_baranyi_with_lag(
  gr_curve,
  LOG10N0 = NULL,
  init_lag = NULL,
  init_mumax = NULL,
  init_LOG10Nmax = NULL,
  algorithm = "auto",
  max_iter = 100,
  lower_bound = c(0, 0, 0, 0)
)
}
\arguments{
\item{gr_curve}{data from one specific growth curve with these two columns: time and biomass}

\item{LOG10N0}{the decimal logarithm of initial biomass}

\item{init_lag}{initial value for the lag parameter}

\item{init_mumax}{initial value for the mumax parameter}

\item{init_LOG10Nmax}{initial value for the LOG10Nmax parameter}

\item{algorithm}{defaults to "auto" which chooses between bounded and unbounded Levenberg-Marquardt method and the bounded port method}

\item{max_iter}{max. number of itertaions; defaults to 100}

\item{lower_bound}{lower.bound for the bounded nls optimisation; defaults to 0}
}
\value{
lag and the nls fitting object with parameters fitted to logistic model
}
\description{
Runs nlsLM/nls algorithms with three different parameter setups to fit the best Logistic model parameters to our data and chooses the best model
}
