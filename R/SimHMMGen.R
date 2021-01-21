#'HMM observations simulation
#'
#'@title Simulation of univariate hidden Markov model
#'
#'@description This function simulates observation from a univariate hidden Markov model
#'
#'@param Q  transition probability matrix; (r x r)
#'@param family   distribution name; run the function distributions() for help
#'@param theta   parameters; (r x p)
#'@param n    number of simulated observations
#'@param graph 1 for a graph, 0 otherwise (default); only for continuous distributions
#'
#'@return \item{SimData}{Simulated data}
#'@return \item{MC}{Simulated Markov chain}
#'@return \item{Sim}{Simulated Data for each regime}
#'
#'@examples
#'family = "gaussian"
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2) ; theta = matrix(c(-1.5, 1.7, 1, 1),2,2) ;
#'sim = SimHMMGen(Q, family, theta, 500, 0)
#'
#'
#'
#'family = "binomial"
#'size = 5
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2) ; thetaB = matrix(c(size, size, 0.2, 0.7),2,2) ;
#'simB = SimHMMGen(Q, family, thetaB, 500, graph=0)$SimData
#
#'
#'
#'@import ggplot2
#'
#'@export
#'
SimHMMGen<-function(Q, family, theta, n, graph=0){
  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    reg = dim(QQ0)[1]
  } else {
    reg = dim(Q)[2]
  }

  #if(dim(Q)[2] <= 1){
  # warning("Transition matrix not assigned, will simulate one regime copula.")
  #}else
  if(nargs()<=3 ){
    stop("Requires at least four input arguments.")
  }




  if(reg >=2){
    MC = SimMarkovChain(Q,n)
  }else  MC = rep(1,n+1)
  Sim   = matrix(0,n,reg)
  SimData = matrix(0,n,1)


  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]

           for( k in 1:reg){
             Sim[,k] = VaRES::varaep(stats::runif(n), q1 = theta[k,1], q2 = theta[k,2], alpha = theta[k,3])
           }
         } ,


         "asymlaplace" = {    ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::ralap(n, location = theta[k,1], scale = theta[k,2], kappa = theta[k,3])
           }
         } ,


         "asympower" = {    ## [01, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varasypower(stats::runif(n), a = theta[k,1], lambda = theta[k,2], delta = theta[k,3])
           }
         } ,


         "asymt" = {    ## [R+, R+, 01, R]

           for( k in 1:reg){
             Sim[,k] =  theta[k,4] + VaRES::varast(stats::runif(n), nu1 = theta[k,1], nu2 = theta[k,2], alpha = theta[k,3])
           }
         } ,


         "beard" = {    ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbeard(stats::runif(n), a = theta[k,1], b = theta[k,2], rho = theta[k,3])
           }
         } ,


         "benini" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rbenini(n, y0 = theta[k,1], shape = theta[k,2])
           }
         } ,


         "benford" = {     ## [1 ou 2]

           for( k in 1:reg){
             Sim[,k] = VGAM::rbenf(n, ndigits = theta[k,1])
           }
         } ,


         "bernoulli" = {     ## [01]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rbern(n, prob = theta[k,1])
           }
         } ,


         "beta" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rbeta(n, shape1 = theta[k,1], shape2 = theta[k,2])
           }
         } ,


         "betabinomial" = {     ## [N+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rbbinom(n, size = theta[k,1], alpha = theta[k,2], beta = theta[k,3])
           }
         } ,


         "betageometric" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rbetageom(n, shape1 = theta[k,1], shape2 = theta[k,2])
           }
         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rbnbinom(n, size = theta[k,1], alpha = theta[k,2], beta = theta[k,3])
           }
         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetaburr(stats::runif(n), a = theta[k,1], b = theta[k,2], c = theta[k,3], d = theta[k,4])
           }
         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetaburr7(stats::runif(n), a = theta[k,1], b = theta[k,2], c = theta[k,3], k = theta[k,4])
           }
         } ,


         "betaexponential" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetaexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2], b = theta[k,3])
           }
         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetafrechet(stats::runif(n), a = theta[k,1], b = theta[k,2], alpha = theta[k,3],
                                             sigma = theta[k,4])
           }
         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetagompertz(stats::runif(n), b = theta[k,1], c = theta[k,2], d = theta[k,3],
                                              eta = theta[k,4])
           }
         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetagumbel(stats::runif(n), a = theta[k,1], b = theta[k,2], mu = theta[k,3],
                                            sigma = theta[k,4])
           }
         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetagumbel2(stats::runif(n), a = theta[k,1], b = theta[k,2], c = theta[k,3],
                                             d = theta[k,4])
           }
         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetalognorm(stats::runif(n), a = theta[k,1], b = theta[k,2], mu = theta[k,3],
                                             sigma = theta[k,4])
           }
         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetalomax(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,



         "betanormal" = {     ## [R+, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rbetanorm(n, shape1 = theta[k,1], shape2 = theta[k,2],
                                 mean = theta[k,3], sd = theta[k,4])
           }
         } ,


         "betaprime" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rbetapr(n, shape1 = theta[k,1], shape2 = theta[k,2], scale = theta[k,3])
           }
         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varbetaweibull(stats::runif(n), a = theta[k,1], b = theta[k,2], alpha = theta[k,3],
                                             sigma = theta[k,4])
           }
         } ,


         "bhattacharjee" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rbhatt(n, mu = theta[k,1], sigma = theta[k,2], a = theta[k,3])
           }
         } ,


         "binomial" = {     ## [N+, 01]

           for( k in 1:reg){
             Sim[,k] = stats::rbinom(n, size = theta[k,1], prob = theta[k,2])
           }
         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rfatigue(n, alpha = theta[k,1], beta = theta[k,2], mu = theta[k,3])
           }
         } ,


         "boxcox" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rboxcox(n, m = theta[k,1], s = theta[k,2], f = theta[k,3])
           }
         } ,


         "burr" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rburr(n, shape1 = theta[k,1], shape2 = theta[k,2], scale = theta[k,3])
           }
         } ,


         "burr2param" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varburr(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,


         "cauchy" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rcauchy(n, location = theta[k,1], scale = theta[k,2])
           }
         } ,


         "chen" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varchen(stats::runif(n), b = theta[k,1], lambda = theta[k,2])
           }
         } ,


         "chi" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = EnvStats::rchi(n, df = theta[k,1])
           }
         } ,


         "chisquared" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = stats::rchisq(n, df = theta[k,1])
           }
         } ,


         "clg" = {     ## [R+, R+, R]
            thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varclg(stats::runif(n), a = thetaS[k,1], b = thetaS[k,2], thetaS[k,3])
           }
         } ,


         "complementarybeta" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varcompbeta(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,



         "dagum" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rdagum(n, scale = theta[k,1], shape1.a = theta[k,2], shape2.p = theta[k,3])
           }
         } ,


         "diffzeta" = {     ## [R+, >1]

           for( k in 1:reg){
             Sim[,k] = VGAM::rdiffzeta(n, shape = theta[k,1], start = theta[k,2])
           }
         } ,


         "discretegamma" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rdgamma(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "discretelaplace" = {     ## [R, 01]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rdlaplace(n, location = theta[k,1], scale = theta[k,2])
           }
         } ,


         "discretenormal" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rdnorm(n, mean = theta[k,1], sd = theta[k,2])
           }
         } ,


         "discreteweibull" = {     ## [01, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rdweibull(n, shape1 = theta[k,1], shape2 = theta[k,2])
           }
         } ,



         "doubleweibull" = {     ## [R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::vardweibull(stats::runif(n), c = theta[k,1], mu = theta[k,2], sigma = theta[k,3])
           }
         } ,


         "ev" = {

           for( k in 1:reg){     ## [R, R+]
             Sim[,k] = VGAM::rgev(n, location = theta[k,1], scale = theta[k,2], shape = 0)
           }
         } ,


         "exponential" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = stats::rexp(n, rate = theta[k,1])
           }
         } ,


         "exponentialextension" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexpext(stats::runif(n), lambda = theta[k,1], a = theta[k,2])
           }
         } ,



         "exponentialgeometric" = {     ## [R+, 01]

           for( k in 1:reg){
             Sim[,k] = VGAM::rexpgeom(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]

           for( k in 1:reg){
             Sim[,k] = VGAM::rexplog(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "exponentialpoisson" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexppois(stats::runif(n), b = theta[k,1], lambda = theta[k,2])
           }
         } ,


         "exponentialpower" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexppower(stats::runif(n), mu = theta[k,1], sigma = theta[k,2], a = theta[k,3])
           }
         } ,


         "exponentiatedexponential" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexpexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2])
           }
         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexplogis(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varexpweibull(stats::runif(n), a = theta[k,1], alpha = theta[k,2], sigma = theta[k,3])
           }
         } ,


         "F" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rf(n, df1 = theta[k,1], df2 = theta[k,2])
           }
         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rfpareto(n, min = theta[k,1], shape1 = theta[k,2],
                                        shape2 = theta[k,3], shape3 = theta[k,4],
                                        scale = theta[k,5])
           }
         } ,


         "fisk" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rfisk(n, scale = theta[k,1], shape1.a = theta[k,2])
           }
         } ,


         "foldednormal" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rfoldnorm(n, mean = theta[k,1], sd = theta[k,2])
           }
         } ,


         "frechet" = {     ## [R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rfrechet(n, shape = theta[k,1], location = theta[k,2], scale = theta[k,3])
           }
         } ,


         "gamma" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rgamma(n, shape =  theta[k,1], scale = theta[k,2])
           }
         } ,


         "gammapoisson" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rgpois(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "gaussian" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rnorm(n, mean = theta[k,1], sd = theta[k,2])
           }
         } ,


         "gev" = {     ## [R, R+, R]

           for( k in 1:reg){
             Sim[,k] = VGAM::rgev(n, location = theta[k,1], scale = theta[k,2], shape = theta[k,3])
           }
         } ,


         "geninvbeta" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::vargeninvbeta(stats::runif(n), a = theta[k,1], c = theta[k,2], d = theta[k,3])
           }
         } ,


         "genlogis" = {     ## [R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::vargenlogis(stats::runif(n), a = theta[k,1], mu = theta[k,2], sigma = theta[k,3])
           }
         } ,


         "genlogis3" = {     ## [R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::vargenlogis3(stats::runif(n), a = theta[k,1], mu = theta[k,2], sigma = theta[k,3])
           }
         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::vargenlogis4(stats::runif(n), a = theta[k,1], alpha = theta[k,2], mu = theta[k,3], sigma = theta[k,4])
           }
         } ,


         "genpowerweibull" = {     ## [R+, R+]
            thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::vargenpowerweibull(stats::runif(n), a = thetaS[k,1], theta = thetaS[k,2])
           }
         } ,



         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

           for( k in 1:reg){
             Sim[,k] = GeneralizedHyperbolic::rghyp(n, mu = theta[k,1], delta = theta[k,2],
                                                    alpha = theta[k,3], beta = theta[k,4],
                                                    lambda = theta[k,5])
           }
         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           for( k in 1:reg){
             Sim[,k] = GLDEX::rgl(n, lambda1 = theta[k,1], lambda2 = theta[k,2], lambda3 = theta[k,3],
                                  lambda4 = theta[k,4])
           }
         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = gamlss.dist::rGT(n, mu = theta[k,1], sigma = theta[k,2], nu = theta[k,3], tau = theta[k,4])
           }
         } ,


         "geometric" = {     ## [01]

           for( k in 1:reg){
             Sim[,k] = stats::rgeom(n, prob = theta[k,1])
           }
         } ,


         "gompertz" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = ssdtools::rgompertz(n, lscale = theta[k,1], lshape = theta[k,2])
           }
         } ,


         "gpd" = {     ## [R, R+, R]

           for( k in 1:reg){
             Sim[,k] = VGAM::rgpd(n, location = theta[k,1], scale = theta[k,2], shape = theta[k,3])
           }
         } ,


         "gumbel" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rgumbel(n, location = theta[k,1], scale = theta[k,2])
           }
         } ,


         "gumbel2" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rgumbelII(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "halfcauchy" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rhcauchy(n, sigma = theta[k,1])
           }
         } ,


         "halflogistic" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varhalflogis(stats::runif(n), lambda = theta[k,1])
           }
         } ,


         "halfnormal" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rhnorm(n, sigma = theta[k,1])
           }
         } ,


         "halft" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rht(n, nu = theta[k,1], sigma = theta[k,2])
           }
         } ,


         "hjorth" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rhjorth(n, m = theta[k,1], s = theta[k,2], f = theta[k,3])
           }
         } ,


         "hblaplace" = {     ## [01, R, R+]
            thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varHBlaplace(stats::runif(n), a = thetaS[k,1], theta = thetaS[k,2], phi = thetaS[k,3])
           }
         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]

           for( k in 1:reg){
             Sim[,k] = GeneralizedHyperbolic::rhyperb(n, mu = theta[k,1], delta = theta[k,2],
                                                    alpha = theta[k,3], beta = theta[k,4])
           }
         } ,


         "huber" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rhuber(n, mu = theta[k,1], sigma = theta[k,2])
           }
         } ,


         "hzeta" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rhzeta(n, shape = theta[k,1])
           }
         } ,


         "inversebeta" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varinvbeta(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,


         "inverseburr" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvburr(n, shape1 = theta[k,1], shape2 = theta[k,2], scale = theta[k,3])
           }
         } ,


         "inversechisquared" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rinvchisq(n, nu = theta[k,1])
           }
         } ,


         "inverseexponential" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvexp(n, scale = theta[k,1])
           }
         } ,


         "inverseexpexponential" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varinvexpexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2])
           }
         } ,


         "inversegamma" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rinvgamma(n, alpha = theta[k,1], beta = theta[k,2])
           }
         } ,


         "inverselomax" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rinv.lomax(n, scale = theta[k,1], shape2.p = theta[k,2])
           }
         } ,


         "inverseparalogistic" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvparalogis(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "inversepareto" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvpareto(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvtrgamma(n, shape1 = theta[k,1], shape2 = theta[k,2], scale = theta[k,3])
           }
         } ,


         "inverseweibull" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rinvweibull(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "kumaraswamy" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rkumar(n, shape1 = theta[k,1], shape2 = theta[k,2])
           }
         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2], b = theta[k,3])
           }
         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumgamma(stats::runif(n), a = theta[k,1], b = theta[k,2], c = theta[k,3], d = theta[k,4])
           }
         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumgumbel(stats::runif(n), a = theta[k,1], b = theta[k,2], mu = theta[k,3],
                                    sigma = theta[k,4])
           }
         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumhalfnorm(stats::runif(n), sigma = theta[k,1], a = theta[k,2], b = theta[k,3])
           }
         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumloglogis(stats::runif(n), a = theta[k,1], b = theta[k,2], alpha = theta[k,3],
                                      beta = theta[k,4])
           }
         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumnormal(stats::runif(n), mu = theta[k,1], sigma = theta[k,2], a = theta[k,3],
                                           b = theta[k,4])
           }
         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varkumweibull(stats::runif(n), a = theta[k,1], b = theta[k,2], alpha = theta[k,3],
                                      sigma = theta[k,4])
           }
         } ,



         "laplace" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rlaplace(n, mu = theta[k,1], sigma = theta[k,2])
           }
         } ,


         "levy" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rlevy(n, m = theta[k,1], s = theta[k,2])
           }
         } ,


         "linearfailurerate" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varlfr(stats::runif(n), a = theta[k,1], b = theta[k,2])
           }
         } ,


         "lindley" = {     ## [R+]
           thetaS=theta
           for( k in 1:reg){
             Sim[,k] = VGAM::rlind(n, theta = thetaS[k,1])
           }
         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varLNbeta(stats::runif(n), lambda = theta[k,1], a = theta[k,2], b = theta[k,3])
           }
         } ,


         "logcauchy" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varlogcauchy(stats::runif(n), mu = theta[k,1], sigma = theta[k,2])
           }
         } ,



         "loggamma" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rlgamma(n, location = theta[k,1], scale = theta[k,2], shape = theta[k,3])
           }
         } ,


         "loggumbel" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = ssdtools::rlgumbel(n, llocation = theta[k,1], lscale = theta[k,2])
           }
         } ,


         "loglaplace" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rloglap(n, location.ald = theta[k,1], scale.ald = theta[k,2], kappa = theta[k,3])
           }
         } ,


         "loglog" = {     ## [R+, >1]

           for( k in 1:reg){
             Sim[,k] = VaRES::varloglog(stats::runif(n), a = theta[k,1], lambda = theta[k,2])
           }
         } ,


         "loglogistic" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rllogis(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "lognormal" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rlnorm(n, meanlog = theta[k,1], sdlog = theta[k,2])
           }
         } ,


         "lognormal3" = {     ## [R, R+, R]

           for( k in 1:reg){
             Sim[,k] = EnvStats::rlnorm3(n, meanlog = theta[k,1], sdlog = theta[k,2], threshold = theta[k,3])
           }
         } ,


         "logistic" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rlogis(n, location = theta[k,1], scale = theta[k,2])
           }
         } ,


         "logisticexponential" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varlogisexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2])
           }
         } ,


         "logisticrayleigh" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varlogisrayleigh(stats::runif(n), a = theta[k,1], lambda = theta[k,2])
           }
         } ,


         "logseries" = {     ## [01]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rlgser(n, theta = theta[k,1])
           }
         } ,


         "lomax" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rlomax(n, scale = theta[k,1], shape3.q = theta[k,2])
           }
         } ,


         "makeham" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rmakeham(n, scale = theta[k,1], shape = theta[k,2], epsilon = theta[k,3])
           }
         } ,


         "maxwell" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rmaxwell(n, rate = theta[k,1])
           }
         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]
           thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varMlaplace(stats::runif(n), theta = thetaS[k,1], phi = thetaS[k,2], psi = thetaS[k,3])
           }
         } ,


         "moexponential" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varmoexp(stats::runif(n), lambda = theta[k,1], a = theta[k,2])
           }
         } ,


         "moweibull" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varmoweibull(stats::runif(n), a = theta[k,1], b = theta[k,2], lambda = theta[k,3])
           }
         } ,


         "nakagami" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rnaka(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "ncchisquared" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rchisq(n, df = theta[k,1], ncp = theta[k,2])
           }
         } ,


         "ncF" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rf(n, df1 = theta[k,1], df2 = theta[k,2], ncp = theta[k,3])
           }
         } ,


         "negativebinomial" = {     ## [N+, 01]

           for( k in 1:reg){
             Sim[,k] = stats::rnbinom(n, size = theta[k,1], prob = theta[k,2])
           }
         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]

           for( k in 1:reg){
             Sim[,k] = GeneralizedHyperbolic::rnig(n, mu = theta[k,1], delta = theta[k,2],
                                                   alpha = theta[k,3], beta = theta[k,4])
           }
         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rnsbeta(n, shape1 = theta[k,1], shape2 = theta[k,2],
                                           min = theta[k,3], max = theta[k,4])
           }
         } ,


         "paralogistic" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rparalogistic(n, scale = theta[k,1], shape1.a = theta[k,2])
           }
         } ,


         "pareto" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rpareto(n, a = theta[k,1], b = theta[k,2])
           }
         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varparetostable(stats::runif(n), lambda = theta[k,1], nu = theta[k,2], sigma = theta[k,3])
           }
         } ,


         "pareto1" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rparetoI(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "pareto2" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rparetoII(n, location = theta[k,1], scale = theta[k,2], shape = theta[k,3])
           }
         } ,


         "pareto3" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rparetoIII(n, location = theta[k,1], scale = theta[k,2], inequality = theta[k,3])
           }
         } ,


         "pareto4" = {     ## [R, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rparetoIV(n, location = theta[k,1], scale = theta[k,2], inequality = theta[k,3], shape = theta[k,4])
           }
         } ,


         "perks" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rperks(n, scale = theta[k,1], shape = theta[k,2])
           }
         } ,


         "pctalaplace" = {     ## [R+, R]
           thetaS=theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varPCTAlaplace(stats::runif(n), a = thetaS[k,1], theta = thetaS[k,2])
           }
         } ,


         "poisson" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = stats::rpois(n, lambda = theta[k,1])
           }
         } ,


         "power1" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varpower1(stats::runif(n), a = theta[k,1])
           }
         } ,


         "power2" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varpower2(stats::runif(n), b = theta[k,1])
           }
         } ,


         "powerdistribution" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rpower(n, alpha = theta[k,1], beta = theta[k,2])
           }
         } ,


         "powerexponential" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rpowexp(n, m = theta[k,1], s = theta[k,2], f = theta[k,3])
           }
         } ,


         "rayleigh" = {     ## [R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rrayleigh(n, scale = theta[k,1])
           }
         } ,


         "reflectedgamma" = {     ## [R+, R, R+]
           thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varrgamma(stats::runif(n), a = thetaS[k,1], theta = thetaS[k,2], phi = thetaS[k,3])
           }
         } ,



         "rice" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rrice(n, sigma = theta[k,1], vee = theta[k,2])
           }
         } ,


         "scaledchisquared" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rinvchisq(n, nu = theta[k,1], tau = theta[k,2])
           }
         } ,


         "schabe" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varschabe(stats::runif(n), gamma = theta[k,1], theta = theta[k,2])
           }
         } ,



         "simplex" = {     ## [01, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rsimplex(n, m = theta[k,1], s = theta[k,2])
           }
         } ,


         "skewedlaplace" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rskewlaplace(n, m = theta[k,1], s = theta[k,2], f = theta[k,3])
           }
         } ,


         "skewedt" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = skewt::rskt(n, df = theta[k,1], gamma = theta[k,2])
           }
         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

           for( k in 1:reg){
             Sim[,k] = sn::rst(n, xi = theta[k,1], omega = theta[k,2], alpha = theta[k,3], nu = theta[k,4])
           }
         } ,


         "skewednormal" = {     ## [R, R+, R, R]

           for( k in 1:reg){
             Sim[,k] = sn::rsn(n, xi = theta[k,1], omega = theta[k,2], alpha = theta[k,3])
           }
         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           for( k in 1:reg){
             Sim[,k] = gamlss.dist::rSEP(n, mu = theta[k,1], sigma = theta[k,2], nu = theta[k,3], tau = theta[k,4])
           }
         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           for( k in 1:reg){
             Sim[,k] = sgt::rsgt(n, mu = theta[k,1], sigma = theta[k,2], lambda = theta[k,3], p = theta[k,4])
           }
         } ,


         "slash" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rslash(n, mu = theta[k,1], sigma = theta[k,2])
           }
         } ,


         "stable" = {     ## [02, -1+1, R+, R]

           for( k in 1:reg){
             Sim[,k] = stabledist::rstable(n, alpha = theta[k,1], beta = theta[k,2], gamma = theta[k,3],
                               delta = theta[k,4])
           }
         } ,


         "stacy" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = VaRES::varstacygamma(stats::runif(n), gamma = theta[k,1], c = theta[k,2], theta = theta[k,3])
           }
         } ,


         "t" = {     ## [R, R+, R+]

           for( k in 1:reg){
             Sim[,k] = theta[k,1] + theta[k,2]*stats::rt(n, df = theta[k,3])
           }
         } ,



         "tobit" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rtobit(n, mean = theta[k,1], sd = theta[k,2])
           }
         } ,


         "topple" = {     ## [01]

           for( k in 1:reg){
             Sim[,k] = VGAM::rtopple(n, shape = theta[k,1])
           }
         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rtrbeta(n, shape1 = theta[k,1], shape2 = theta[k,2], shape3 = theta[k,3],
                                       scale = theta[k,4])
           }
         } ,


         "transformedgamma" = {     ## [R+, R+, R+]

           for( k in 1:reg){
             Sim[,k] = actuar::rtrgamma(n, shape1 = theta[k,1], shape2 = theta[k,2], scale = theta[k,3])
           }
         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rtnorm(n, mean = theta[k,1], sd = theta[k,2], a = theta[k,3],
                                          b = theta[k,4])
           }
         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]

           for( k in 1:reg){
             Sim[,k] = VGAM::rtruncpareto(n, lower = theta[k,1], upper = theta[k,2], shape = theta[k,3])
           }
         } ,


         "twosidedpower" = {     ## [01, R+]

           for( k in 1:reg){
             Sim[,k] = rmutil::rtwosidedpower(n, m = theta[k,1], s =  theta[k,2])
           }
         } ,


         "wald" = {     ## [R, R+]

           for( k in 1:reg){
             Sim[,k] = extraDistr::rwald(n, mu = theta[k,1], lambda = theta[k,2])
           }
         } ,


         "weibull" = {     ## [R+, R+]

           for( k in 1:reg){
             Sim[,k] = stats::rweibull(n, shape = theta[k,1], scale = theta[k,2])
           }
         } ,


         "xie" = {     ## [R+, R+, R+]
           thetaS = theta
           for( k in 1:reg){
             Sim[,k] = VaRES::varxie(stats::runif(n), a = thetaS[k,1], b = thetaS[k,2], lambda = thetaS[k,3])
           }
         } ,


         "yules" = {     ## [R+] mais > 0.5

           for( k in 1:reg){
             Sim[,k] = VGAM::ryules(n, shape = theta[k,1])
           }
         } ,


  )


  for(i in 1:n){
    k = MC[i]
    SimData[i,] = Sim[i,k]
  }

  if (graph != 0){
    ddd = data.frame("observations" = SimData)
    ddd$regimes = as.factor(MC)
    ddd$index = c(1:n)
    print(ggplot2::ggplot(ddd, ggplot2::aes(ddd$index, ddd$observations), color="black") +
            ggplot2::geom_line(size=0.5) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = "white"),
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                  plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title="Simulated data")
    )

    print(ggplot2::ggplot(ddd, ggplot2::aes(ddd$regimes, ddd$observations, color=ddd$regimes)) +
            ggplot2::geom_boxplot() +
            ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2)) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = "white"),
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                  legend.position = c(.02, .95),
                  legend.justification = c("left", "top"),
                  legend.box = "horizontal",
                  legend.margin = ggplot2::margin(6, 6, 6, 6),
                  legend.box.background = ggplot2::element_rect(color="black", size=1),
                  legend.box.margin = ggplot2::margin(1, 1, 1, 1),
                  legend.title = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title="Simulated regimes of each observation", x="regimes")
    )


    print(ggplot2::ggplot(ddd, ggplot2::aes(ddd$regimes, ddd$observations, color=ddd$regimes, fill = ddd$regimes)) +
            ggplot2::geom_violin() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = "white"),
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                  legend.position = c(.05, .95),
                  legend.justification = c("left", "top"),
                  legend.box = "horizontal",
                  legend.margin = ggplot2::margin(6, 6, 6, 6),
                  legend.box.background = ggplot2::element_rect(color="black", size=1),
                  legend.box.margin = ggplot2::margin(1, 1, 1, 1),
                  legend.title = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title="Simulated regimes of each observation", x="regimes")
    )


    print(ggplot2::ggplot(ddd, ggplot2::aes(ddd$regimes, ddd$observations, color=ddd$regimes, fill = ddd$regimes)) +
            ggplot2::geom_point() +
            ggplot2::geom_boxplot() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = "white"),
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                  legend.position = c(.02, .95),
                  legend.justification = c("left", "top"),
                  legend.box = "horizontal",
                  legend.margin = ggplot2::margin(6, 6, 6, 6),
                  legend.box.background = ggplot2::element_rect(color="black", size=1),
                  legend.box.margin = ggplot2::margin(1, 1, 1, 1),
                  legend.title = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title="Simulated regimes of each observation", x="regimes")
    )


  }


  out=list(SimData=SimData,MC=MC,Sim=Sim)
  return(out)
}

