#'@title Density function of a mixture of Gaussian univariate distributions
#'
#'@description This function computes the density function of a mixture of
#'Gaussian univariate distributions
#'
#'@param x       Points at which the density is comptuted (n x 1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param w       vector of the probability of each regime (r x 1).
#
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{f}{Values of the distribution function}
#'
#'@examples mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); w <-c(0.8, 0.2);
#'x <- seq(-1, 1, by = 0.01)
#'f <- GaussianMixturePdf(x,mu,sigma,w)
#'plot(x,f,type="l")
#'
#'@importFrom stats dnorm
#'@export
#
GaussianMixturePdf=function(x,mu,sigma,w){

r = length(w);

n = length(x);

D = diag(1./sigma);

z = matrix(0,n,r);

for (k in 1:r){
z[,k] = (x-mu[k])/sigma[k];
}

f = as.vector(dnorm(z)%*%D%*%w)

return(f)

}

