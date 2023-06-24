#'@title Simulation of a univariate Gaussian Hidden Markov Model (HMM)
#'
#'@description Generates a univariate regime-switching random walk with Gaussian regimes starting from a given state eta0, using the inverse method from noise u.Can be useful when generating multiple time series.
#'
#'@param u       series of uniform i.i.d. series (n x 1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param Q       Transition probality matrix (r x r);
#'@param eta0    Initial value for the regime;
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'@references Nasri & Remillard (2019). Copula-based dynamic models for multivariate time series. JMVA, vol. 172, 107--121.
#'@return \item{x}{Simulated Data}
#'@return \item{eta}{Probability of regimes}
#'
#'@examples 
#'Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) 
#'set.seed(1)
#'u <-runif(250)
#'mu <- c(-0.3 ,0.7) 
#'sigma <- c(0.15,0.05);
#'eta0=1
#'x <- SimHMMGaussianInv(u,mu,sigma,Q,eta0)
#'
#'@importFrom stats runif dnorm
#'
#'@export
#'
SimHMMGaussianInv  = function(u,mu,sigma,Q,eta0)
{

  n=length(u)
  
  x = rep(0,n)
  r = length(mu)
  eta=matrix(0,nrow=n,ncol=r)
  f = rep(0,r)
  eps = qnorm(u)
  w = Q[eta0,]
  
  
  for(i in 1:n)
  {
   x[i] = GaussianMixtureInv(u[i],mu,sigma,w);
  
  for(k in 1:r){f[k] = w[k] * dnorm(x[i],mu[k],sigma[k])}
  
  eta[i,] = f/sum(f);
  w = t(eta[i,]%*% Q) 
  }
  out=list(x=x,eta=eta)
  out
}
 

  