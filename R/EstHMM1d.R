#'@title Estimation of a univariate Gaussian Hidden Markov Model (HMM)
#'
#' @description This function estimates parameters (mu, sigma, Q) of a univariate Hidden Markov Model.
#' It computes also the probability of being in each regime, given the past observations (eta)
#' and the whole series (lambda). The conditional distribution given past observations is applied to
#' obtains pseudo-observations W that should be uniformly distributed under the null hypothesis.
#' A Cramér-von Mises test statistic is then computed.
#'
#'@param y   (nx1) vector of data
#'@param reg    number of regimes
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10 000
#'@param eps precision (stopping criteria); suggestion 0.0001.
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{mu}{estimated mean for each regime}
#'@return \item{sigma}{stimated standard deviation for each regime}
#'@return \item{Q}{(reg x reg) estimated transition matrix}
#'@return \item{eta}{(n x reg) probabilities of being in regime k at time t given observations up to time t}
#'@return \item{lambda}{(n x reg) probabilities of being in regime k at time t given all observations}
#'@return \item{cvm}{Cramér-von Mises statistic for the goodness-of-fit test}
#'@return \item{U}{Pseudo-observations that should be uniformly distributed under the null hypothesis of a Gaussian HMM}
#'@return \item{LL}{Log-likelihood}
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05)
#'data <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,100)$x
#'est <- EstHMM1d(data, 2, max_iter=10000, eps=0.0001)
#'
#'@importFrom stats dnorm pnorm   sd
#'@export
#'
EstHMM1d <- function(y,reg,max_iter=10000,eps=0.0001){



  n = length(y)

  out0=.C("est_hmm_1d_new",
                  as.double(y),
                  as.integer(reg),
                  as.integer(n),
                  as.double(eps),
                  as.double(max_iter),
                  mu=double(reg),
                  sigma = double(reg),
                  nu = double(reg),
                  Qvec = double(reg*reg),
                  etavec = double(n*reg),
                  lambdavec = double(n*reg),
                  U = double(n),
                  cvm = double(1),
                  L  = double(n))

eta = matrix(out0$etavec,ncol=reg)
lambda = matrix(out0$lambdavec,ncol=reg)
Q= matrix(out0$Qvec,ncol=reg)

mu <- out0$mu
sigma <- out0$sigma
nu     <- out0$nu
U     <- out0$U
cvm = out0$cvm
LL = sum(log(out0$L))
if(is.na(sum(U))){out=NULL;
cat("Singular values. Number of regimes too large.")}else{
out=list(mu=mu,sigma=sigma,Q=Q,eta=eta,lambda=lambda,nu=nu,U=U,cvm=cvm,LL=LL)}
return(out)
}
