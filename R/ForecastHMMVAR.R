#'@title Value at risk (VAR) of a univariate HMM at time n+k1, n+k2, ...
#'
#'@description This function computes the VAR of a univariate HMM for multiple horizons, given observations up to time n
#'
#'@param U  value (1  x 1) between 0 and 1
#'@param family    distribution name; run the function distributions() for help
#'@param theta  parameters; (r  x p)
#'@param Q    probability transition  matrix; (r  x r)
#'@param eta    vector of the estimated probability of each regime at time n; (1  x r)
#'@param k  prediction times (may be a vector of integers).
#'
#'@return \item{var}{values at risk (1 x horizon)}
#'
#'
#'@examples
#'\donttest{
#'family = "gaussian"
#'
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96091218, 0.03908782)
#'
#'
#'forecastedVAR = ForecastHMMVAR(U=c(0.01), family, theta=theta, Q=Q, eta=eta, k=c(1,2,3,4,5))
#'print('Forecasted VAR : ')
#'print(forecastedVAR)
#'}
#'
#'
#'
#'@export


ForecastHMMVAR<-function(U, family, theta, Q, eta, k=1){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }


  var = matrix(0, nrow=length(U), ncol=length(k))

  for (d in 1:length(k)){
    Q_prime = matrixcalc::matrix.power(Q, as.integer(k[d]))
    for (l in 1:r){
      for (j in 1:r){
        var[,d] = var[,d] + eta[j] * Q_prime[j, l] * QUANTILE(family, U, theta[l,])
      }
    }
  }



  return(var)


}


