#'@title Goodness-of-fit test of a univariate Gaussian Hidden Markov Model
#'
#'@description This function performs a goodness-of-fit test of a Gaussian HMM based on a Cramér-von Mises statistic
#'using parametric bootstrap.
#'
#'
#'@param       y         (n x 1) data vector
#'@param       reg       number of regimes
#'@param       max_iter  maxmimum number of iterations of the EM algorithm; suggestion 10 000
#'@param       prec      precision (stopping criteria); suggestion 0.0001
#'@param       n_sample  number of bootstrap samples; suggestion 1000
#'@param       n_cores   number of cores to use in the parallel computing
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{pvalue}{pvalue of the Cram\'er-von Mises statistic in percent}
#'@return \item{mu}{estimated mean for each regime}
#'@return \item{sigma}{estimated standard deviation for each regime}
#'@return \item{Q}{(reg x reg) estimated transition matrix}
#'@return \item{eta}{(n x reg) conditional probabilities of being in regime k at time t given observations up to time t}
#'@return \item{lambda}{(n x reg) probabilities of being in regime k at time t given all observations}
#'@return \item{cvm}{Cramér-von Mises statistic for the goodness-of-fit test}
#'@return \item{W}{Pseudo-observations that should be uniformly distributed under the null hypothesis of a Gaussian HMM}
#'@return \item{LL}{Log-likelihood}
#'
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'@examples
#'\donttest{
#'Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05)
#'data <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,100)$x
#'gof <- GofHMM1d(data, 2, max_iter=10000, prec=0.0001, n_sample=100,n_cores=2)
#'}
#'@importFrom doParallel registerDoParallel
#'@importFrom parallel makeCluster stopCluster
#'@import     foreach
#'@importFrom stats dnorm pnorm  rnorm qnorm sd
#'
#'@export
#'


GofHMM1d <-function(y, reg, max_iter=10000, prec=0.0001, n_sample=1000,n_cores){
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)


  #n_cores = detectCores()-2;

  out0 = EstHMM1d(y,reg,max_iter, prec);
  print("End of estimation");

  mu       = out0$mu;
  sigma    = out0$sigma;
  Q        = out0$Q;
  cvm      = out0$cvm;
  W        = out0$W;
  eta      = out0$eta;
  lambda   = out0$lambda;
  nu       = out0$nu;
  LL       = out0$LL;

  n= length(y);
 fun = c('Sim.Markov.Chain','Sim.HMM.Gaussian.1d','EstHMM1d','Sn','bootstrapfun','em.step')

 # result <- foreach(i=1:n_sample, .packages='GaussianHMM1d') %dopar% bootstrapfun(mu,sigma,Q,max_iter,prec,n)
  result <- foreach(i=1:n_sample, .export=fun) %dopar% bootstrapfun(mu,sigma,Q,max_iter,prec,n)

  cvm_sim = rep(0,n_sample)
  for (i in 1:n_sample){
    cvm_sim[i] = result[[i]]
  }
  stopCluster(cl)

  pvalue = 100*mean( cvm_sim > cvm)

  out = list(pvalue =  pvalue, mu=mu, sigma=sigma, Q=Q, eta = eta, lambda=lambda, nu = nu, cvm = cvm, W = W)
  return(out)
}


