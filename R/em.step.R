#'@title Function to perform the E-M steps for the estimation of the paramters
#'
#'@description This function perform the E-M steps for the estimation of the parameters
#'of a univariate Gaussian HMM.
#'
#'@param y       points at which the density function is comptuted (mx1);
#'@param mu      vector of means for each regime (r x 1);
#'@param sigma   vector of standard deviations for each regime (r x 1);
#'@param  Q      transition probality matrix (r x r);
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
#'out <- em.step(data,mu,sigma,Q)
#'@importFrom stats dnorm
#'
#'@export
#'
em.step   <- function(y,mu,sigma,Q){

  n <- length(y)
  r <- length(mu)
  gamma    <- matrix(0,n,r)
  gammabar <- matrix(0,n,r)
  lambda   <- matrix(0,n,r)
  f        <- matrix(0,n,r)
  Lambda   <- array(0, dim=c(r,r,n))
  ##
  for(j in 1:r){
    z      <- (y-mu[j])/sigma[j]
    f[,j] <- dnorm(z)/sigma[j]
  }

  ##
  gammabar[n,]<-1/r

  for(k in 1: (n-1)){
    i <- n-k
    j=i+1
    v =  ( gammabar[j,] * f[j,] ) %*% t(Q)
    gammabar[i,] = v/sum(v)
    #gammabar[i,] <- ( gammabar[i+1,] * f[i+1,] ) * t(Q)
  }

  ##
  L=matrix(0,1,n)
  eta=matrix(0,n,r)
  eta0 <- matrix(1,1,r)/r

  v <- ( eta0 %*% Q) * f[1,]
  L[1]=sum(v)
  eta[1,] <- v/sum(v)

  for(i in 2:n){
    v        <- ( eta[i-1,] %*% Q) * f[i,]
    L[i]=sum(v)
    eta[i,] <- v/L[i]
  }

  LL=sum(log(L))
  ##
  v      <- eta * gammabar
  sv0=rowSums(v)
  for(j in 1:(r)){lambda[,j]=v[,j]/sv0}

  #sv     <- matrix(1,1,r) %*% apply(v,1,sum)
  #lambda <- v / sv

  ##
  #Lambda[,,n] <- matrix(1,r,1) %*% lambda[n,] * Q

  #for (j in (1:r)){  Lambda[j,,n]=lambda[n,j] %*% Q[j,]  }
  for (j in (1:r)){  Lambda[j,,n]=lambda[n,j]*Q[j,]  }

  gf=gammabar*f


  for(i in 1:(n-1)){
    M <- Q*( matrix(1,r,1) *(eta[i,]))%*%( matrix( 1,1,r) * gf[i+1,] )
    c <- sum(apply(M,2,sum))
    Lambda[,,i] <- M/c
  }

  ##
  nu <- apply(lambda,2,mean)

  munew=matrix(0,1,r)
  sigmanew=matrix(0,1,r)
  Qnew=matrix(0,r,r)


  for (j in 1:r)
  {
    w=lambda[,j]/nu[j]/n
    #yy <- matrix(1,1,r) %*% y
    munew[j]<- sum( y * w )
    sigmanew[j] <- sqrt( sum( y^2 * w) - munew[j]^2 )
    Qnew[j,] <- apply((Lambda[j,,]),1,mean) /  nu[j]
  }



  # w  <- lambda / matrix(1,n,1)%*%nu /n
  #yy <- matrix(1,1,r) %*% y
  #munew <- sum( yy * w )
  #sigmanew <- sqrt( apply( yy^2 * w,2,sum) - munew^2 )
  #Qnew <- apply(Lambda,3,mean) / matrix(1,1,r)%*% t(nu)

  return(list(nu=nu, munew=munew, sigmanew=sigmanew, Qnew=Qnew, eta=eta, gammabar=gammabar, lambda=lambda, Lambda=Lambda, LL=LL))
}
