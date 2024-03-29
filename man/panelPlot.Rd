% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{panelPlot}
\alias{panelPlot}
\title{Plots implied and actual correlations}
\usage{
panelPlot(pcOutput)
}
\arguments{
\item{pcOutput}{An object created from running the \code{panelcoder()} command.}
}
\description{
\code{panelPlot()} plots implied and actual stability coefficients for
increasingly long lags for the variables analyzed by the panelcoder
function. The function takes a panelcoder output object as an argument and
then plots these stabilities. Currently, the model only works for observed-
variable models. Eventually it will be able to handle latent-variable models,
too.
}
