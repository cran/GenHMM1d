#'@title Tail conditional median function
#'
#'@description This function computes the tail conditional median of a univariate distribution
#'
#'@param family   distribution name; run the command distributions() for help
#'@param p (1  x 1) values  between 0 and 1
#'@param param  parameters of the distribution
#'@param size additional parameter for some discrete distributions
#'@param Nsim number of simulations
#'
#'@return \item{tcm}{tail conditional median}
#'
#'@examples
#'family = "gaussian"
#'
#'Q = 1 ; theta = matrix(c(-1.5, 1.7),1,2) ;
#'tcm = TCM(family, (0.01), theta)
#'print('Tail conditional mean : ')
#'print(tcm$tcm)
#'
#'
#'@export



TCM<-function(family, p, param, size=0, Nsim=25000){

  var = QUANTILE(family, p, param, size)

  param_sim = matrix(param, nrow = 1, ncol = length(param))

  sim = SimHMMGen(1, family, param_sim, Nsim)

  diffp = sim$SimData-var

  diffp[diffp>0] = 0

  tcm = var + (1/(1-p)) * stats::median(replace(diffp, diffp == 0, NA), na.rm = TRUE)

  out = list(tcm=tcm, var=var, sim=sim)
  return(out)


}
