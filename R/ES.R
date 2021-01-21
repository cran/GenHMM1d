#'@title Expected shortfall function
#'
#'@description This function compute the expected shortfall of an univariate distribution
#'
#'@param family   distribution name; run the function distributions() for help
#'@param p  value (1  x 1) at which the expected shortfall needs to be computed; between 0 and 1; (e.g 0.01, 0.05)
#'@param param  parameters of the distribution; (1 x p)
#'@param size additional parameter for some discrete distributions; run the command distributions() for help
#'@param Nsim number of simulations
#'
#'@return \item{es}{expected shortfall}
#'
#'@examples
#'family = "gaussian"
#'
#'theta = matrix(c(-1.5, 1.7),1,2) ;
#'es = ES(family, (0.01), theta)
#'print('Expected shortfall : ')
#'print(es$es)
#'
#'
#'@export



ES<-function(family, p, param, size=0, Nsim=25000){

  var = QUANTILE(family, p, param, size)

  param_sim = matrix(param, nrow = 1, ncol = length(param))

  sim = SimHMMGen(1, family, param_sim, Nsim)

  diffp = sim$SimData-var

  diffp[diffp>0] = 0

  es = var + (1/(1-p)) * base::mean(replace(diffp, diffp == 0, NA), na.rm = TRUE)

  out = list(es=es, var=var, sim=sim)
  return(out)


}
