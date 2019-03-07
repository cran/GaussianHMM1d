#'@title Estimated probabilities of the regimes given new observations
#'
#'@description This function computes the estimated probabilities of the regimes for a Gaussian HMM
#' given new observation after time n. it also computes the associated weight of the Gaussian mixtures
#' that can be used for forecasted density, cdf, or quantile function.
#'
#'@param ynew    new observations (mx1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param  Q      transition probality matrix (r x r);
#'@param eta     vector of the estimated probability of each regime (r x 1) at time n;
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{etanew}{values of the estimated probabilities at times n+1 to n+m, using the new observations}
#'@return \item{w}{weights of the mixtures for periods n+1 to n+m}
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'@examples mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); eta <- c(.1,.9);
#'x <- c(0.2,-0.1,0.73)
#'out <- ForecastHMMeta(x,mu,sigma,Q,eta)
#'
#'@importFrom stats dnorm
#'@export
#'
ForecastHMMeta=function(ynew,mu,sigma,Q,eta){

  m      = length(ynew);
  r      = length(mu);
  ww     = matrix(0,m+1,r);
  etanew = matrix(0,m,r);
  f      = matrix(0,m,r);

  for(i in 1:r){
    f[,i] = dnorm(ynew,mu[i],sigma[i]) + 1E-10;
  }

  ww[1,] =  eta %*% Q ;
  for(k in 1:m){
    num        = f[k,] * ww[k,];
    etanew[k,] = num/sum(num);
    ww[k+1,] =  etanew[k,] %*% Q ;
  }

  w = ww[-(m+1),]
  out = list(etanew=etanew,w=w)
  return(out)

}
