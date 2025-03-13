#'@title Value at risk (VAR) of a univariate HMM at time n+k1, n+k2, ...
#'
#'@description This function computes the VAR of a univariate HMM for multiple horizons, given observations up to time n
#'
#'@param U       values (n  x 1) between 0 and 1
#'@param ZI      1 if zero-inflated, 0 otherwise (default)
#'@param family  distribution name; run the function distributions() for help
#'@param theta   parameters; (r  x p)
#'@param Q       probability transition  matrix for the regimes; (r  x r)
#'@param eta     vector of the estimated probability of each regime at time n; (1  x r)
#'@param k      prediction times (may be a vector of integers).
#'
#'@return \item{var}{values at risk (1 x horizon)}
#'
#'
#'@examples
#'\donttest{
#'family = "gaussian"
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96, 0.04)
#'U=c(0.01,0.05)
#' k=c(1,2,3,4,5)
#'ForecastHMMVAR(U, 0, family, theta, Q, eta=eta,k)
#'}
#'
#'
#'
#'@export


ForecastHMMVAR<-function(U, ZI=0, family, theta, Q, eta, k=1){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }

  un=1+ZI;

  n = length(U)
  var  = matrix(0, nrow=n, ncol=length(k))
  var0 = matrix(0, nrow=n, ncol=r)



  for (l in un:r){
         var0[,l] = QUANTILE( U, theta[l,],family)
  }
  if(ZI==1) var0[,1] = 0

  ub = apply(var0,1,max)
  lb = apply(var0,1,min)-ZI

  for (d in 1:length(k)){

    for(i in 1:n)
    {
       p=U[i];
       a = lb[i]
       b = ub[i]

      x0 = (a+b)/2;

      for (j in 1:20){
        u0 = ForecastHMMCdf(x0, ZI, family, theta, Q, eta, k[d])
        if(u0<p){a = x0}else {b=x0}
        x0 = (a+b)/2;
      }

      var[i,d] = x0
    }


  }



  return(var)


}
