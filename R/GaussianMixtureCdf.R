#'@title Distribution function of a mixture of Gaussian univariate distributions
#'
#'@description This function computes the distribution function of a mixture of
#' Gaussian univariate distributions
#'
#'@param x       Points at which the distribution function is comptuted (nx1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param w      vector of the probability of each regime (r x r).
#
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{F}{values of the distribution function}
#'
#'@examples mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); w <-c(0.8, 0.2);
#'x <- seq(-1, 1, by = 0.01)
#'F <- GaussianMixtureCdf(x,mu,sigma,w)
#'plot(x,F,type="l")
#'
#'@importFrom stats pnorm
#'@export
#'
GaussianMixtureCdf=function(x,mu,sigma,w){

r = length(w);

n = length(x);

z = matrix(0,n,r);

for (k in 1:r){
z[,k] = (x-mu[k])/sigma[k]
}

F = as.vector(pnorm(z)%*%w);
return(F)

}
