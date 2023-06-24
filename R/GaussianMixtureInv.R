#'@title Inverse distribution function of a mixture of Gaussian univariate distributions
#'
#'@description This function computes the inverse distribution function of a mixture of
#' Gaussian univariate distributions
#'
#'@param p       Points in (0,1) at which the distribution function is computed (nx1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param w       vector of the probability of each regime (r x 1).
#
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{q}{values of the quantile function}
#'
#'@examples mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); w <-c(0.8, 0.2);
#'p <- seq(0.01, 0.99, by = 0.01)
#'q <- GaussianMixtureInv(p,mu,sigma,w)
#'plot(p,q,type="l")
#'
#'@importFrom stats qnorm
#'@export
#'
GaussianMixtureInv=function(p,mu,sigma,w){


n = length(p)
q = rep(0,n);
for(i in 1:n)
{
u = p[i];

a = min(mu)+min(sigma*qnorm(u));

b = max(mu)+max(sigma*qnorm(u));

x0 = (a+b)/2;

for (k in 1:20){
u0 = GaussianMixtureCdf(x0,mu,sigma,w);
if(u0<u){a = x0}else {b=x0}
x0 = (a+b)/2;
}

q[i] = x0
}
return(q)
}

