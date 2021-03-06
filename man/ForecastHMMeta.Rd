% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ForecastHMMeta.R
\name{ForecastHMMeta}
\alias{ForecastHMMeta}
\title{Estimated probabilities of the regimes given new observations}
\usage{
ForecastHMMeta(ynew, mu, sigma, Q, eta)
}
\arguments{
\item{ynew}{new observations (mx1);}

\item{mu}{vector of means for each regime (r x 1);}

\item{sigma}{vector of standard deviations for each regime (r x 1);}

\item{Q}{transition probality matrix (r x r);}

\item{eta}{vector of the estimated probability of each regime (r x 1) at time n;}
}
\value{
\item{etanew}{values of the estimated probabilities at times n+1 to n+m, using the new observations}

\item{w}{weights of the mixtures for periods n+1 to n+m}
}
\description{
This function computes the estimated probabilities of the regimes for a Gaussian HMM
given new observation after time n. it also computes the associated weight of the Gaussian mixtures
that can be used for forecasted density, cdf, or quantile function.
}
\examples{
mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); eta <- c(.1,.9);
x <- c(0.2,-0.1,0.73)
out <- ForecastHMMeta(x,mu,sigma,Q,eta)

}
\references{
Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
}
\author{
Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
}
