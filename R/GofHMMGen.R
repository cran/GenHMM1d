#'@title Goodness-of-fit of univariate hidden Markov model
#'
#'@description This function performs goodness-of-fit test of an univariate hidden Markov model
#'
#'
#'@param y   observations
#'@param reg    number of regimes
#'@param family   distribution name; run the function distributions() for help
#'@param start  starting parameter for the estimation
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10000
#'@param eps precision (stopping criteria); suggestion 0.0001.
#'@param graph 1 for a graph, 0 otherwise (default); only for continuous distributions
#'@param size additional parameter for some discrete distributions; run the command distributions() for help
#'@param n_sample number of bootstrap samples; suggestion 1000
#'@param n_cores number of cores to use in the parallel computing
#'@param useFest 1 (default) to use the first estimated parameters as starting value for the bootstrap, 0 otherwise
#'
#'
#'@return \item{pvalue}{pvalue of the Cramer-von Mises statistic in percent}
#'@return \item{theta}{Estimated parameters; (r x p) }
#'@return \item{Q}{estimated transition matrix; ; (r x r)}
#'@return \item{eta}{(conditional probabilities of being in regime k at time t given observations up to time t; (n x r)}
#'@return \item{lambda}{conditional probabilities of being in regime k at time t given all observations; (n x r)}
#'@return \item{U}{matrix of Rosenblatt transforms; (n x r)}
#'@return \item{cvm}{Cramer-von-Mises statistic for goodness-of-fit}
#'@return \item{W}{pseudo-observations that should be uniformly distributed under the null hypothesis}
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
#'@export
#'
#'
#'
GofHMMGen <-function(y, reg, family, start=0, max_iter=10000, eps=10e-4, graph=0, size=0, n_sample=100, n_cores=1, useFest=1){
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)

  if (family != "negativebinomial" && family != "binomial" && family != "betanegativebinomial" && family != "betabinomial") {
    size=0
  }

  esthmmgen = EstHMMGen(y, reg=reg, family=family, start=start, max_iter=max_iter, eps=eps, size=size)
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

  eta0 = sample(1:reg, n_sample, replace = T)
  n = length(y);

  print("First estimation done")

  result <- foreach::foreach(i=1:n_sample, .packages="GenHMM1d") %dopar% bootstrapfun(Q, family, theta, n, size, max_iter, eps, useFest)


  cvm_sim1 = rep(0,n_sample)
  for (i in 1:n_sample){
    cvm_sim1[i] = result[[i]]$cvm_sim
  }

  pvalue = 100*mean( stats::na.omit(cvm_sim1 > cvm_est))

  out = list(  pvalue=pvalue, theta=theta, Q=Q, eta=eta, nu=nu, U=U, cvm=cvm_est, cvms = cvm_sim1, W=W,
               lambda=lambda, LL=LL, AIC=AIC, BIC=BIC, stats=stats, pred_l=pred_l, pred_e=pred_e,
               runs_e=runs_e, runs_l=runs_l)


  parallel::stopCluster(cl)


  return(out)
}
