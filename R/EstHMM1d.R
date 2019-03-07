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
#'@param prec precision (stopping criteria); suggestion 0.0001.
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{mu}{estimated mean for each regime}
#'@return \item{sigma}{stimated standard deviation for each regime}
#'@return \item{Q}{(reg x reg) estimated transition matrix}
#'@return \item{eta}{(n x reg) probabilities of being in regime k at time t given observations up to time t}
#'@return \item{lambda}{(n x reg) probabilities of being in regime k at time t given all observations}
#'@return \item{cvm}{Cramér-von Mises statistic for the goodness-of-fit test}
#'@return \item{W}{Pseudo-observations that should be uniformly distributed under the null hypothesis of a Gaussian HMM}
#'@return \item{LL}{Log-likelihood}
#'
#'@references Chapter 10.2 of B. Rémillard (2013). Statistical Methods for Financial Engineering,
#'Chapman and Hall/CRC Financial Mathematics Series, Taylor & Francis.
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2); mu <- c(-0.3 ,0.7) ; sigma <- c(0.15,0.05)
#'data <- Sim.HMM.Gaussian.1d(mu,sigma,Q,eta0=1,100)$x
#'est <- EstHMM1d(data, 2, max_iter=10000, prec=0.0001)
#'
#'@importFrom stats dnorm pnorm   sd
#'@export
#'
EstHMM1d <- function(y,reg,max_iter=10000,prec=0.0001){
n <- length(y)
r <- reg


if(r==1){
  mu     <- mean(y)
  sigma  <- sd(y)
  Q      <- 1
  eta    <- rep(1,n)
  lambda <- eta

  Z = (y-mu)/sigma
  L = dnorm(Z)/sigma
  W = pnorm(Z);

  LL     <- sum(log(L))
  cvm    <- Sn(W)
  nu    <- 1} else{

n0 <- floor(n/r)

x <- matrix(0,n0,r)

ind0 <- (1:n0)

for(j in 1:r){
    ind <- (j-1)*n0 + ind0
    x[,j] <- y[ind]
}

mu0    <- apply(x,2,mean)
sigma0 <- apply(x,2,sd)
Q0     <- matrix(1,r,r)/r

for(k in 1:100){
  curr.step <- em.step(y,mu0,sigma0,Q0)
  #eta <- curr.step$eta
  mu0    <- curr.step$mu
  Q0     <- curr.step$Q
  sigma0 <- curr.step$sigma
}


rprec= r*prec;


for(k in 1:max_iter){
  curr.step <- em.step(y,mu0,sigma0,Q0)
  eta <- curr.step$eta

  sum1 <- sum(abs(mu0))
  sum2 <- sum(abs(curr.step$mu-mu0))


  if ( sum2 < (sum1*rprec) ){
    break
  }
  mu0    <- curr.step$mu
  Q0     <- curr.step$Q
  sigma0 <- curr.step$sigma


}


mu <- curr.step$mu
Q     <- curr.step$Q
sigma <- curr.step$sigma
LL<-curr.step$LL
lambda<-curr.step$lambda
nu=curr.step$nu

Z=matrix(0,n,r)

for (i in 1:r){Z[,i]=(y-mu[i])/sigma[i]}

U = pnorm(Z);

eta00 = matrix(1,1,r)/r;

w00 = rbind(eta00, eta) %*% Q;

dd=dim(w00)[1]

w   = w00[-dd,];

W = as.vector(rowSums( w * U));

cvm = Sn(W);
}

return(list(mu=mu,sigma=sigma,Q=Q,eta=eta,nu=nu,lambda=lambda,LL=LL,cvm=cvm, W=W))
}

