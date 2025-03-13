#'@title Bootstrap for the univariate distributions
#'
#'@param theta     parameters
#'@param Q         transition matrix for the regimes
#'@param ZI        1 if zero-inflated, 0 otherwise (default)
#'@param family    distribution name; (run the command distributions() for help)
#'@param n         number of simulated observations
#'@param size      additional parameter for some discrete distributions; run the command distributions() for help
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10000
#'@param eps       precision (stopping criteria); suggestion 0.0001
#'@param useFest   TRUE (default) to use the first estimated parameters as starting value for the bootstrap, FALSE otherwise
#'
#'
#'@return Internal function used for the parametric bootstrap
#'
#'@export
#'@keywords internal

bfun <- function(theta,Q, ZI, family, n, size=0, max_iter=10000, eps=1e-4, useFest=TRUE){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }


   y = SimHMMGen(theta, size, Q, ZI, family, n)$SimData

  if (useFest){
      esthmmg = EstHMMGen(y, ZI, reg, family, start=0, max_iter=max_iter, eps=eps, size=size, theta0=theta)
  } else {
      esthmmg = EstHMMGen(y, ZI, reg, family, max_iter=max_iter, eps=eps, size=size)

  }

  out = list(theta1=esthmmg$theta, Q1=esthmmg$Q, eta1=esthmmg$eta, nu1=esthmmg$nu, U1=esthmmg$U,
             cvm_sim=esthmmg$cvm, W1=esthmmg$W,lambda1=esthmmg$lambda, LL1=esthmmg$LL, AIC1=esthmmg$AIC,
             BIC1=esthmmg$BIC, stats1=esthmmg$stats, pred_l1=esthmmg$pred_l, pred_e1=esthmmg$pred_e,
             runs_e1=esthmmg$runs_e, runs_l1=esthmmg$runs_l)

  return(out)

}
