#'@title Cramer-von Mises statistic for goodness-of-fit of the null hypothesis of a univariate uniform distrubtion over [0,1]
#'
#'@description This function computes the Cramér-von Mises statistic Sn for goodness-of-fit of the null hypothesis of a univariate uniform distrubtion over [0,1]
#'
#'@param U vector of pseudos-observations (apprimating uniform variates)
#'
#'@author Bouchra R Nasri  and Bruno N Rémillard, January 31, 2019
#'
#'@return \item{Sn}{Cramér-von Mises statistic}
#'
#'@export
Sn = function(U)
{

n   = length(U)
l   = (2*seq(1,n)-1)/n
U0  = sort(U)

stat=n/3+sum((U0-l)*U0)

return(stat)

}
