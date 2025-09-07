#'@title Goodness-of-fit of univariate hidden Markov model
#'
#'@description This function performs a goodness-of-fit test for a univariate hidden Markov model
#'
#'
#'@param y         observations
#'@param ZI        1 if zero-inflated, 0 otherwise (default)
#'@param reg       number of regimes
#'@param family    distribution name; run the function distributions() for help
#'@param start     starting parameter for the estimation
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10000
#'@param eps       precision (stopping criteria); suggestion 0.0001.
#'@param size      additional parameter for some discrete distributions; run the command distributions() for help
#'@param n_samples number of bootstrap samples; suggestion 1000
#'@param n_cores   number of cores to use in the parallel computing
#'
#'
#'@return \item{pvalue}{pvalue of the Cramer-von Mises statistic in percent}
#'@return \item{theta}{Estimated parameters; (r x p) }
#'@return \item{Q}{estimated transition matrix; ; (r x r)}
#'@return \item{eta}{(conditional probabilities of being in regime k at time t given observations up to time t; (n x r)}
#'@return \item{lambda}{conditional probabilities of being in regime k at time t given all observations; (n x r)}
#'@return \item{U}{pseudo-observations that should be uniformly distributed under the null hypothesis}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
#'@return \item{W}{matrix of Rosenblatt transforms; (n x r)}
#'@return \item{LL}{log-likelihood}
#'@return \item{nu}{stationary distribution}
#'@return \item{AIC}{Akaike information criterion}
#'@return \item{BIC}{bayesian information criterion}
#'@return \item{CAIC}{consistent Akaike information criterion}
#'@return \item{AICcorrected}{Akaike information criterion corrected}
#'@return \item{HQC}{Hannan-Quinn information criterion}
#'@return \item{stats}{Empirical means and standard deviation of each regimes using lambda}
#'@return \item{pred_l}{Estimated regime using lambda}
#'@return \item{pred_e}{Estimated regime using eta}
#'@return \item{runs_l}{Estimated number of runs using lambda}
#'@return \item{runs_e}{Estimated number of runs using eta}
#'
#'@importFrom foreach %dopar%
#'@import doParallel
#'
#'
#'@examples
#'family = "gaussian"
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2) ; theta = matrix(c(0, 1.7, 0, 1),2,2) ;
#'y = SimHMMGen(theta, size=0, Q, ZI=1, family,  100)$SimData
#'out=GofHMMGen(y,1,2,family,n_samples=10)
#'
#'
#'@export
#'
#'
#'
GofHMMGen <-function(y, ZI=0, reg, family, start=0, max_iter=10000, eps=1e-4, size=0, n_samples=1000, n_cores=1){
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)

  if( (family %in% c("negativebinomial","binomial", "betanegativebinomial", "betabinomial")) == FALSE){size=0}

  esthmmgen = EstHMMGen(y, ZI, reg, family, start, max_iter, eps, size)
  theta = esthmmgen$theta
  Q = esthmmgen$Q
  eta = esthmmgen$eta
  nu = esthmmgen$nu
  U = esthmmgen$U
  cvm_est = esthmmgen$cvm
  W = esthmmgen$W
  lambda = esthmmgen$lambda
  LL = esthmmgen$LL
  AIC =esthmmgen$AIC
  BIC = esthmmgen$BIC
  CAIC =esthmmgen$CAIC
  AICcorrected = esthmmgen$AICcorrected
  HQC = esthmmgen$HQC
  stats = esthmmgen$stats
  pred_l = esthmmgen$pred_l
  pred_e = esthmmgen$pred_e
  runs_e = esthmmgen$runs_e
  runs_l = esthmmgen$runs_l

  eta0 = sample(1:reg, n_samples, replace = T)
  n = length(y);



  result <- foreach::foreach(i=1:n_samples, .packages="GenHMM1d") %dopar% bfun(theta, Q, ZI, family, start, n, size, max_iter, eps)


  cvm_sim1 = rep(0,n_samples)
  for (i in 1:n_samples){
    cvm_sim1[i] = result[[i]]$cvm_sim
  }

  pvalue = 100*mean( stats::na.omit(cvm_sim1 > cvm_est))

  out = list(  pvalue=pvalue, theta=theta, Q=Q, eta=eta, nu=nu, U=U, cvm=cvm_est, cvms = cvm_sim1, W=W,
               lambda=lambda, LL=LL, AIC=AIC, BIC=BIC, stats=stats, pred_l=pred_l, pred_e=pred_e,
               runs_e=runs_e, runs_l=runs_l)


  parallel::stopCluster(cl)


  return(out)
}
