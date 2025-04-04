% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aj_estimate_helpers.R
\name{extract_aj_estimate_by_assumptions}
\alias{extract_aj_estimate_by_assumptions}
\title{Extract adjusted estimate by assumptions}
\usage{
extract_aj_estimate_by_assumptions(
  data_to_adjust,
  fixed_time_horizons,
  censoring_assumption = "excluded",
  competing_assumption = "excluded"
)
}
\arguments{
\item{data_to_adjust}{A data frame contains original estimated risks,
estimated risks stratas, the time to events, and the stratification technique}

\item{fixed_time_horizons}{Fixed time horizons.}

\item{censoring_assumption}{One of `"excluded"` or `"adjusted"`.}

\item{competing_assumption}{One of `"excluded"` or `"adjusted"`.}
}
\value{
An adjusted estimate.
}
\description{
A short description...
}
