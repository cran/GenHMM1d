#'@title Predicted probabilities of regimes of a univariate HMM for a new observation
#'
#'@description This function computes the predicted probabilities of the regimes for a new observation of a univariate HMM, given observations up to time n
#'
#'@param ynew   new observations
#'@param ZI     1 if zero-inflated, 0 otherwise (default)
#'@param family distribution name; run the function distributions() for help
#'@param theta  parameters; (r  x p)
#'@param Q      probability transition  matrix for the regimes; (r  x r)
#'@param eta    vector of the estimated probability of each regime at time n; (1  x r)
#'
#'
#'@return \item{etanew}{predicted probabilities of the regimes}
#'
#'@examples
#'family = "gaussian"
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96, 0.04)
#'ForecastHMMeta(1.5, 0, family, theta, Q, eta)
#'
#'@export


ForecastHMMeta<-function(ynew, ZI=0, family, theta, Q, eta){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }

  un = 1+ZI

  etanew = matrix(0, nrow=length(ynew), ncol=r)

  for (j in un:r){
    for (l in 1:r){
      etanew[,j] = etanew[,j] + eta[l] * Q[l, j]
    }
    etanew[,j] = etanew[,j] * PDF(family, ynew, theta[j,])
  }

  if(ZI==1)
  {
    for (l in 1:r){
      etanew[,1] = etanew[,1] + eta[l] * Q[l, 1]
    }
    etanew[,1] = etanew[,1] * (ynew==0)
  }

  etanew = etanew / rowSums(etanew)

  return(etanew)


}


