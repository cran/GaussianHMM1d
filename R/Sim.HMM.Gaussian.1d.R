#'@title Simulation of a univariate Gaussian Hidden Markov Model (HMM)
#'
#'@description This function simulates observations from a univariate Gaussian HMM
#'
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param Q       Transition probality matrix (r x r);
#'@param eta0    Initial value for the regime;
#'@param n       number of simulated observations.
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{x}{Simulated Data}
#'@return \item{reg}{Markov chain regimes}
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ; mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05);
#'sim <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,n=100)
#'
#'@importFrom stats rnorm
#'
#'@export
Sim.HMM.Gaussian.1d  = function(mu,sigma,Q,eta0, n)
{


  r=length(mu)
  x=matrix(0,n,1)

  reg = Sim.Markov.Chain(Q,n,eta0);

  x0 = matrix(0,n,r);

  for (k in 1:r){
  x0[,k] = mu[k]+ sigma[k]*rnorm(n,1);
  }

  for (i in 1:n){
  x[i] = x0[i,reg[i] ];
  }

  out = list(x=x,reg=reg);
}
