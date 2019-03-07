#'@title Density function of a Gaussian HMM at time n+k
#'
#'@description This function computes the density function of a Gaussian HMM
#' at time n+k, given observation up to time n.
#'
#'@param x       points at which the density function is comptuted (mx1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param  Q      transition probality matrix (r x r);
#'@param eta     vector of the estimated probability of each regime (r x 1) at time n;
#'@param k       time of prediction.
#
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{f}{values of the density function at time n+k}
#'@return \item{w}{weights of the mixture}
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'@examples mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05); Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ;
#'eta <- c(.9,.1);
#'x <- seq(-1, 1, by = 0.01)
#'out <- ForecastHMMPdf(x,mu,sigma,Q,eta,3)
#'plot(x,out$f,type="l")
#'
#'@export
#'
ForecastHMMPdf=function(x,mu,sigma,Q,eta,k){

  Q0=Q;
  if(k>1){
  for(i in 2:k)
  {Q0=Q0 %*% Q;}
  }
  w =  as.vector(eta %*% Q0);

  f = GaussianMixturePdf(x,mu,sigma,w);
  out = list(f=f,w=w)
  return(out)

}
