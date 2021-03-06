% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initializemcmcST.R
\name{initialize_mcmc_ST}
\alias{initialize_mcmc_ST}
\title{initialize_mcmc_ST}
\usage{
initialize_mcmc_ST(ns, STu, MCMC, Nisolates)
}
\arguments{
\item{ns}{The number of sources}

\item{STu}{A vector of the unique ST types}

\item{MCMC}{The number of mcmc iterations}

\item{Nisolates}{The number of isolates}
}
\value{
A list to be passed to the simulator
}
\description{
A function to initialize the mcmc in R for the ST model
}
\author{
Jukka Ranta
}

