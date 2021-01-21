#'@title Cramer-von Mises statistic for the goodness-of-fit test of the null hypothesis of a univariate uniform distribution over [0,1]
#'
#'@description This function computes the Cramer-von Mises statistic Sn for goodness-of-fit of the null hypothesis of a univariate uniform distrubtion over [0,1]
#'
#'@param U vector of pseudos-observations (approximating uniform)
#'
#'@return \item{sta}{Cramer-von Mises statistic}
#'
#'@export


Snd1 <- function(U){

  n = length(U)

  u = sort(U)

  t = t((-0.5+rep(1:n))) / n

  sta = 1/(12*n) + sum( (u-t)^2 ) ;

  return(sta)

}
