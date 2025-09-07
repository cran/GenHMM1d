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
#'
#'
#'@return Internal function used for the parametric bootstrap
#'
#'@export
#'@keywords internal

bfun <- function(theta,Q, ZI, family, start, n, size, max_iter, eps){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }


   ysim = SimHMMGen(theta, size, Q, ZI, family, n)$SimData


      esthmmg = EstHMMGen(ysim, ZI, reg, family, start, max_iter, eps, size)



  out = list(theta1=esthmmg$theta, Q1=esthmmg$Q, eta1=esthmmg$eta, nu1=esthmmg$nu, U1=esthmmg$U,
             cvm_sim=esthmmg$cvm, W1=esthmmg$W,lambda1=esthmmg$lambda, LL1=esthmmg$LL, AIC1=esthmmg$AIC,
             BIC1=esthmmg$BIC, stats1=esthmmg$stats, pred_l1=esthmmg$pred_l, pred_e1=esthmmg$pred_e,
             runs_e1=esthmmg$runs_e, runs_l1=esthmmg$runs_l)

  return(out)

}
