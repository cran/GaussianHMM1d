#'@title Estimated Regimes for the univariate Gaussian HMM
#'
#'@description This function computes and plots the most likely regime for univariate Gaussian HMM using
#'probabilities of being in regime k at time t given all observations (lambda)
#'and probabilities of being in regime k at time t given observations up to time t (eta).
#'@param t     (nx1) vector of dates (years, ...); if no dates then t=[1:length(y)]
#'@param  y     (nx1) vector of data;
#'@param lambda (nxreg) probabilities of being in regime k at time t given all observations;
#'@param eta    (nxreg) probabilities of being in regime k at time t given observations up to time t;
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{A}{Estimated Regime using lambda}
#'@return \item{B}{Estimated Regime using eta}
#'@return \item{runsA}{Estimated number of runs using lambda}
#'@return \item{runsB}{Estimated number of runs using eta}
#'@return \item{pA}{Graph for the estimated regime for each observation using lambda}
#'@return \item{pB}{Graph for the estimated regime for each observation using eta}
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'@examples  Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05);
#'data <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,100)$x
#'t=c(1:100);
#'est <- EstHMM1d(data, 2)
#'EstRegime(t,data,est$lambda, est$eta)
#'
#'@importFrom graphics plot
#'@export
#

EstRegime=function(t,y,lambda,eta){

  n=dim(lambda)[1]

  A=matrix(0,n,1)
  B=matrix(0,n,1)

for (i in 1:n){
A[i,]=which.max(lambda[i,])}

  runsA=length(rle(c(A))$length)

  for (i in 1:n){
    B[i,]=which.max(eta[i,])}

  runsB=length(rle(c(B))$length)

  pA=plot(t,y,col=A, main = "Estimated regime for univariate Gaussian HMM using lambda")
  pB=plot(t,y,col=B, main = "Estimated regime for univariate Gaussian HMM using eta")
  return(list(A=A, B=B, runsA=runsA, runsB=runsB))
  return(pA)
  return(pB)


}
