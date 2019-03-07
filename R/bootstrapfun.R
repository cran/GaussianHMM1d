#'@title Function to perform parametric bootstrap
#'
#'@description This function simulates the data under the  null hypothesis of a Gaussian HMM
#'and compute the Cramér-von Mises test statistic.
#'
#'@param mu         vector of means for each regime (r x 1);
#'@param sigma      vector of standard deviations for each regime (r x 1);
#'@param  Q         transition probality matrix (r x r);
#'@param max_iter   maximum number of iterations of the EM algorithm; suggestion 10 000;
#'@param prec       precision (stopping criteria); suggestion 0.0001;
#'@param n          length of the time series.
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
#'data <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,100)$x
#'out <- bootstrapfun(mu,sigma,Q,max_iter=10000,prec=0.0001,n=100)
#'
#'@keywords internal
#'
#'@export
#'
bootstrapfun <- function(mu,sigma,Q,max_iter,prec,n){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }


  eta0 = sample(reg,1);  # initial regime chosen at random

  if(reg==1){y1 = mu+sigma*rnorm(n)}else{
    y1 = Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0,n)$x;}

  out2 = EstHMM1d(y1,reg,max_iter, prec);
  # mu1       = out2$mu;
  # sigma1    = out2$sigma;
  # eta1      = out2$eta;
  # Q1        = out2$Q;
  cvm_sim   = out2$cvm;
  #out = list(mu1=mu1 , sigma1=sigma1, Q1=Q1, cvm_sim=cvm_sim)
  #out = list(cvm_sim=cvm_sim)
  return(cvm_sim)

}
