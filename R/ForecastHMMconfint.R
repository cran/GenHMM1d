#'@title Forecasted confidence interval of a univariate HMM at times n+k1, n+k2,....
#'
#'@description This function computes the forecasted confidence interval of a univariate HMM for multiple horizons, given observations up to time n
#'
#'@param U  values between 0 and 1
#'@param family    distribution name; run the function distributions() for help
#'@param theta  parameters; (r  x p)
#'@param Q   probability transition  matrix; (r  x r)
#'@param eta    vector of the estimated probability of each regime at time n; (1 x r)
#'@param k   prediction times (may be a vector of integers).
#'
#'
#'@return \item{qlow}{lower bound of the forecasted confidence interval}
#'@return \item{qhigh}{upper bound of the forecasted confidence interval}
#'
#'@examples
#'family = "gaussian"
#'
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96091218, 0.03908782)
#'
#'forecastedhmmconfint = ForecastHMMconfint(U=c(0.1, 0.9), family, theta=theta, Q=Q,
#' eta=eta, k=c(1,2,3,4,5))
#'print('Forecasted confidence interval : ')
#'print(forecastedhmmconfint)
#'
#'@export


ForecastHMMconfint<-function(U, family, theta, Q, eta, k=1){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }


  qlow = matrix(0, nrow=length(k), ncol=1)
  qhigh = matrix(0, nrow=length(k), ncol=1)

  for (d in 1:length(k)){
    Q_prime = matrixcalc::matrix.power(Q, as.integer(k[d]))
    for (l in 1:r){
      for (j in 1:r){
        qlow[d] = qlow[d] + eta[j] * Q_prime[j, l] * QUANTILE(family, U[1], theta[l,])

        qhigh[d] = qhigh[d] + eta[j] * Q_prime[j, l] * QUANTILE(family, U[2], theta[l,])
      }
    }
  }

  out = list(qlow=qlow, qhigh=qhigh)

  return(out)


}


