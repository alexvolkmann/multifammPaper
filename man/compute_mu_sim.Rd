% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperfunctions_results.R
\name{compute_mu_sim}
\alias{compute_mu_sim}
\title{Compute the Fitted Mean Trajectories}
\usage{
compute_mu_sim(fitted_cu, ran_preds, I = 10, J = 16, reps = 5,
  nested = FALSE, mul = TRUE)
}
\arguments{
\item{fitted_cu}{Object of fitted curves saved from the simulation.}

\item{ran_preds}{Object of random effects saved from the simulation.}

\item{I}{Number of levels for component B. Defaults to 9.}

\item{J}{Number of levels for component C. Defaults to 16.}

\item{reps}{Number of repetitions of the levels of B and C. Defaults to 5.}

\item{nested}{TRUE if the model component C is nested. Defaults to FALSE.}

\item{mul}{FALSE if the univariate estimate is extracted from fitted_cu.
Defaults to TRUE which is the multivariate estimate.}
}
\description{
This function takes the fitted curves object of the
simulation and the random effects object and calculates the estimated mu
functions.
}
