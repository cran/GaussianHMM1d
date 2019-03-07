#'@title Simulation of a finite Markov chain
#'
#'@description This function generates a Markov chain X(1), ..., X(n) with transition matrix Q,
#' starting from a state eta0.
#'
#'@param  Q Transition probality matrix (r x r);
#'@param  n length of series;
#'@param  eta0 inital value in {1,...,r}.
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{x}{Simulated Markov chain}
#'
#'@examples Q <- matrix(c(0.8, 0.3, 0.2, 0.7),2,2) ;
#'  sim <- Sim.Markov.Chain(Q,eta0=1,n=100)
#'
#'@export
#'
Sim.Markov.Chain = function(Q,n,eta0)
{

dd = dim(Q);


r = dd[2];  #number of regimes


x = 0 * c(1:n);
x0 = matrix(0,nrow=n,ncol=r);


ind = eta0;
for (k in 1:r)
  {
     x0[,k] = sample(r,n,replace= TRUE, prob=Q[k,]);
  }

for (i in 1:n)
{
   x[i]=  x0[i,ind];
   ind = x[i];
}

x
}
