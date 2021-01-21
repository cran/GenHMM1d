#'@title Bootstrap for the univariate distributions
#'
#'@param Q transition matrix;
#'@param family    distribution name; (run the command distributions() for help)
#'@param theta   parameters
#'@param n    number of simulated observations
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10000
#'@param size additional parameter for some discrete distributions; run the command distributions() for help
#'@param eps       precision (stopping criteria); suggestion 0.001
#'@param n_sample    number of bootstrap samples; suggestion 1000
#'@param UseFest   1 (default) to use the first estimated parameters as starting value for the bootstrap, 0 otherwise
#'
#'
#'@return Internal function used for the parametric bootstrap
#'
#'@export
#'@keywords internal

bootstrapfun <- function(Q, family, theta, n, size=0, max_iter=10000, eps=10e-4, useFest=1){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }


  if (family == "negativebinomial" || family != "binomial" || family != "betanegativebinomial" || family != "betabinomial") {
    theta = cbind(matrix(size, nrow=dim(theta)[1], ncol=1), theta)
    useFest = 0
  }


  y = SimHMMGen(Q=Q, family=family, theta=theta, n=n)$SimData

  if (useFest==1){
      esthmmg = EstHMMGen(y=y, reg=reg, family=family, start=0, max_iter=max_iter, eps=eps, graph=0, size=size, theta0=theta)
  } else {
      esthmmg = EstHMMGen(y=y, reg=reg, family=family, max_iter=max_iter, eps=eps, size=size)

  }
  print("Bootstrap estimation done +1")
  out = list(theta1=esthmmg$theta, Q1=esthmmg$Q, eta1=esthmmg$eta, nu1=esthmmg$nu, U1=esthmmg$U,
             cvm_sim=esthmmg$cvm, W1=esthmmg$W,lambda1=esthmmg$lambda, LL1=esthmmg$LL, AIC1=esthmmg$AIC,
             BIC1=esthmmg$BIC, stats1=esthmmg$stats, pred_l1=esthmmg$pred_l, pred_e1=esthmmg$pred_e,
             runs_e1=esthmmg$runs_e, runs_l1=esthmmg$runs_l)

  return(out)

}
