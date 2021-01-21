#'@title Predicted probabilities of regimes of a univariate HMM given a new observation
#'
#'@description This function computes the predicted probabilities of the regimes for a new observation of a univariate HMM, given observations up to time n
#'
#'@param ynew  the new observations
#'@param family    distribution name; run the function distributions() for help
#'@param theta  parameters; (r  x p)
#'@param Q    probability transition  matrix; (r  x r)
#'@param eta    vector of the estimated probability of each regime at time n; (1  x r)
#'
#'
#'@return \item{etanew}{predicted probabilities of the regimes}
#'
#'@examples
#'family = "gaussian"
#'
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96091218, 0.03908782)
#'
#'forecastedhmmeta = ForecastHMMeta(c(1.5), family, theta=theta, Q=Q, eta=eta)
#'print('Forecasted regime probabilities : ')
#'print(forecastedhmmeta)
#'
#'@export


ForecastHMMeta<-function(ynew, family, theta, Q, eta){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }

  etanew = matrix(0, nrow=length(ynew), ncol=r)

  for (j in 1:r){
    for (l in 1:r){
      etanew[,j] = etanew[,j] + eta[l] * Q[l, j]
    }
    etanew[,j] = etanew[,j] * PDF(family, ynew, theta[j,])
  }

  etanew = etanew / rowSums(etanew)

  return(etanew)


}


