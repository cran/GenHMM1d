#'@title Quantile function
#'
#'@description This function computes the quantile function of a univariate distribution
#'
#'@param family    distribution name; run the function distributions() for help
#'@param p  values at which the quantile needs to be computed; between 0 and 1; (e.g 0.01, 0.05)
#'@param param  parameters of the distribution; (1 x p)
#'@param size additional parameter for some discrete distributions; run the command distributions() for help
#'
#'@return \item{q}{quantile/VAR}
#'
#'
#'@examples
#'family = "gaussian"
#'
#'Q = 1 ; theta = matrix(c(-1.5, 1.7),1,2) ;
#'quantile = QUANTILE(family, (0.01), theta)
#'print('Quantile : ')
#'print(quantile)
#'
#'@export


QUANTILE<-function(family,p,param,size=0){



  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]


           q = VaRES::varaep(p, q1 = param[1], q2 = param[2], alpha = param[3])

         } ,


         "asymlaplace" = {    ## [R, R+, R+]


           q = VGAM::qalap(p, location = param[1], scale = param[2], kappa = param[3])

         } ,


         "asympower" = {    ## [01, R+, R+]


           q = VaRES::varasypower(p, a = param[1], lambda = param[2], delta = param[3])

         } ,


         "asymt" = {    ## [R+, R+, 01, R]


           q = param[4] + VaRES::varast(p, nu1 = param[1], nu2 = param[2], alpha = param[3])

         } ,


         "beard" = {    ## [R+, R+, R+]


           q = VaRES::varbeard(p, a = param[1], b = param[2], rho = param[3])

         } ,


         "benini" = {     ## [R, R+]


           q = VGAM::qbenini(p, y0 = param[1], shape = param[2])

         } ,


         "benford" = {     ## [1 ou 2]


           q = VGAM::qbenf(p, ndigits = param[1])

         } ,


         "bernoulli" = {     ## [01]


           q = extraDistr::qbern(p, prob = param[1])

         } ,


         "beta" = {     ## [R+, R+]


           q = stats::qbeta(p, shape1 = param[1], shape2 = param[2])

         } ,



         # "betabinomial" = {     ## [N+, R+, R+]
         #
         #   q = extraDistr::qbbinom(p, size = size, alpha = param[1], beta = param[2])
         #
         # } ,


         # "betageometric" = {     ## [R+, R+]
         #
         #   q = VGAM::qbetageom(p, shape1 = param[1], shape2 = param[2])
         #
         # } ,


         # "betanegativebinomial" = {     ## [N+, R+, R+]
         #
         #   q = extraDistr::qbnbinom(p, size = size, alpha = param[1], beta = param[2])
         #
         # } ,


         "betaburr" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetaburr(p, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetaburr7(p, a = param[1], b = param[2], c = param[3], k = param[4])

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


           q = VaRES::varbetaexp(p, lambda = param[1], a = param[2], b = param[3])

         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetafrechet(p, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetagompertz(p, b = param[1], c = param[2], d = param[3], eta = param[4])

         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]


           q = VaRES::varbetagumbel(p, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetagumbel2(p, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]


           q = VaRES::varbetalognorm(p, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetalomax(p, a = param[1], b = param[2], alpha = param[3], lambda = param[4])

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


           q = VGAM::qbetanorm(p, shape1 = param[1], shape2 = param[2],
                               mean = param[3], sd = param[4])

         } ,


         "betaprime" = {     ## [R+, R+, R+]


           q = extraDistr::qbetapr(p, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


           q = VaRES::varbetaweibull(p, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         # "bhattacharjee" = {     ## [R, R+, R+]
         #
         #
         #   q = extraDistr::qbhatt(p, mu = param[1], sigma = param[2], a = param[3])
         #
         # } ,


         "binomial" = {     ## [N+, 01]

           q = stats::qbinom(p, size = size, prob = param[1])

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]


           q = extraDistr::qfatigue(p, alpha = param[1], beta = param[2], mu = param[3])

         } ,


         "boxcox" = {     ## [R+, R+, R+]


           q = rmutil::qboxcox(p, m = param[1], s = param[2], q = param[3])

         } ,


         "burr" = {     ## [R+, R+, R+]


           q = actuar::qburr(p, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "burr2param" = {     ## [R+, R+]


           q = VaRES::varburr(p, a = param[1], b = param[2])

         } ,


         "cauchy" = {     ## [R, R+]


           q = stats::qcauchy(p, location = param[1], scale = param[2])

         } ,


         "chen" = {     ## [R+, R+]


           q = VaRES::varchen(p, b = param[1], lambda = param[2])

         } ,


         "chi" = {     ## [R+]


           q = EnvStats::qchi(p, dq = param[1])

         } ,


         "chisquared" = {     ## [R+]


           q = stats::qchisq(p, dq = param[1])

         } ,


         "clg" = {     ## [R+, R+, R]


           q = VaRES::varclg(p, a = param[1], b = param[2], param[3])

         } ,


         "complementarybeta" = {     ## [R+, R+]


           q = VaRES::varcompbeta(p, a = param[1], b = param[2])

         } ,



         "dagum" = {     ## [R+, R+, R+]


           q = VGAM::qdagum(p, scale = param[1], shape1.a = param[2], shape2.p = param[3])

         } ,


         "diffzeta" = {     ## [R+, >1]


           q = VGAM::qdiffzeta(p, shape = param[1], start = param[2])

         } ,


         # "discretegamma" = {     ## [R+, R+]
         #
         #   q = extraDistr::qdgamma(p, shape = param[1], scale = param[2])
         #
         # } ,


         # "discretelaplace" = {     ## [R, 01]
         #
         #   q = extraDistr::qdlaplace(p, location = param[1], scale = param[2])
         #
         # } ,


         # "discretenormal" = {     ## [R, R+]
         #
         #   q = extraDistr::qdnorm(p, mean = param[1], sd = param[2])
         #
         # } ,


         "discreteweibull" = {     ## [01, R+]

           q = extraDistr::qdweibull(p, shape1 = param[1], shape2 = param[2])

         } ,



         "doubleweibull" = {     ## [R+, R, R+]

           q = VaRES::vardweibull(p, c = param[1], mu = param[2], sigma = param[3])

         } ,


         "ev" = {

           ## [R, R+]
           q = VGAM::qgev(p, location = param[1], scale = param[2], shape = 0)

         } ,


         "exponential" = {     ## [R+]


           q = stats::qexp(p, rate = param[1])

         } ,


         "exponentialextension" = {     ## [R+, R+]


           q = VaRES::varexpext(p, lambda = param[1], a = param[2])

         } ,



         "exponentialgeometric" = {     ## [R+, 01]


           q = VGAM::qexpgeom(p, scale = param[1], shape = param[2])

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


           q = VGAM::qexplog(p, scale = param[1], shape = param[2])

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


           q = VaRES::varexppois(p, b = param[1], lambda = param[2])

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


           q = VaRES::varexppower(p, mu = param[1], sigma = param[2], a = param[3])

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


           q = VaRES::varexpexp(p, lambda = param[1], a = param[2])

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


           q = VaRES::varexplogis(p, a = param[1], b = param[2])

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


           q = VaRES::varexpweibull(p, a = param[1], alpha = param[2], sigma = param[3])

         } ,


         "F" = {     ## [R+, R+]


           q = stats::qf(p, df1 = param[1], df2 = param[2])

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]


           q = actuar::qfpareto(p, min = param[1], shape1 = param[2],
                                shape2 = param[3], shape3 = param[4],
                                scale = param[5])

         } ,


         "fisk" = {     ## [R+, R+]

           q = VGAM::qfisk(p, scale = param[1], shape1.a = param[2])

         } ,


         "foldednormal" = {     ## [R, R+]


           q = VGAM::qfoldnorm(p, mean = param[1], sd = param[2])

         } ,


         "frechet" = {     ## [R+, R, R+]


           q = VGAM::qfrechet(p, shape = param[1], location = param[2], scale = param[3])

         } ,


         "gamma" = {     ## [R+, R+]


           q = stats::qgamma(p, shape =  param[1], scale = param[2])

         } ,


         "gaussian" = {     ## [R, R+]


           q = stats::qnorm(p, mean = param[1], sd = param[2])

         } ,


         "gev" = {     ## [R, R+, R]


           q = VGAM::qgev(p, location = param[1], scale = param[2], shape = param[3])

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


           q = VaRES::vargeninvbeta(p, a = param[1], c = param[2], d = param[3])

         } ,


         "genlogis" = {     ## [R+, R, R+]


           q = VaRES::vargenlogis(p, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis3" = {     ## [R+, R, R+]


           q = VaRES::vargenlogis3(p, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


           q = VaRES::vargenlogis4(p, a = param[1], alpha = param[2], mu = param[3], sigma = param[4])

         } ,


         "genpowerweibull" = {     ## [R+, R+]


           q = VaRES::vargenpowerweibull(p, a = param[1], theta = param[2])

         } ,


         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)


           q = GeneralizedHyperbolic::qghyp(p, mu = param[1], delta = param[2],
                                            alpha = param[3], beta = param[4],
                                            lambda = param[5])

         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           q = GLDEX::qgl(p, lambda1 = param[1], lambda2 = param[2], lambda3 = param[3], lambda4 = param[4])

         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

           q = gamlss.dist::qGT(p, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,


         "geometric" = {     ## [01]

           q = stats::qgeom(p, prob = param[1])

         } ,


         "gompertz" = {     ## [R+, R+]


           q = ssdtools::qgompertz(p, lscale = param[1], lshape = param[2])

         } ,


         "gpd" = {     ## [R, R+, R]


           q = VGAM::qgpd(p, location = param[1], scale = param[2], shape = param[3])

         } ,


         "gumbel" = {     ## [R, R+]


           q = VGAM::qgumbel(p, location = param[1], scale = param[2])

         } ,


         "gumbel2" = {     ## [R+, R+]


           q = VGAM::qgumbelII(p, scale = param[1], shape = param[2])

         } ,


         "halfcauchy" = {     ## [R+]


           q = extraDistr::qhcauchy(p, sigma = param[1])

         } ,


         "halflogistic" = {     ## [R+]


           q = VaRES::varhalflogis(p, lambda = param[1])

         } ,


         "halfnormal" = {     ## [R+]


           q = extraDistr::qhnorm(p, sigma = param[1])

         } ,


         "halft" = {     ## [R+, R+]


           q = extraDistr::qht(p, nu = param[1], sigma = param[2])

         } ,


         "hjorth" = {     ## [R+, R+, R+]


           q = rmutil::qhjorth(p, m = param[1], s = param[2], q = param[3])

         } ,


         "hblaplace" = {     ## [01, R, R+]


           q = VaRES::varHBlaplace(p, a = param[1], theta = param[2], phi = param[3])

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]


           q = GeneralizedHyperbolic::qhyperb(p, mu = param[1], delta = param[2],
                                              alpha = param[3], beta = param[4])

         } ,


         "huber" = {     ## [R, R+]


           q = extraDistr::qhuber(p, mu = param[1], sigma = param[2])

         } ,


         "hzeta" = {     ## [R+]

           q = VGAM::qhzeta(p, shape = param[1])

         } ,


         "inversebeta" = {     ## [R+, R+]


           q = VaRES::varinvbeta(p, a = param[1], b = param[2])

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


           q = actuar::qinvburr(p, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inversechisquared" = {     ## [R+]


           q = extraDistr::qinvchisq(p, nu = param[1])

         } ,


         "inverseexponential" = {     ## [R+]


           q = actuar::qinvexp(p, scale = param[1])

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


           q = VaRES::varinvexpexp(p, lambda = param[1], a = param[2])

         } ,


         "inversegamma" = {     ## [R+, R+]


           q = extraDistr::qinvgamma(p, alpha = param[1], beta = param[2])

         } ,


         "inverselomax" = {     ## [R+, R+]


           q = VGAM::qinv.lomax(p, scale = param[1], shape2.p = param[2])

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


           q = actuar::qinvparalogis(p, shape = param[1], scale = param[2])

         } ,


         "inversepareto" = {     ## [R+, R+]


           q = actuar::qinvpareto(p, shape = param[1], scale = param[2])

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


           q = actuar::qinvtrgamma(p, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inverseweibull" = {     ## [R+, R+]


           q = actuar::qinvweibull(p, shape = param[1], scale = param[2])

         } ,


         "kumaraswamy" = {     ## [R+, R+]


           q = VGAM::qkumar(p, shape1 = param[1], shape2 = param[2])

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


           q = VaRES::varkumexp(p, lambda = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


           q = VaRES::varkumgamma(p, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


           q = VaRES::varkumgumbel(p, a = param[1], b = param[2], mu = param[3],
                                 sigma = param[4])

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


           q = VaRES::varkumhalfnorm(p, sigma = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]


           q = VaRES::varkumloglogis(p, a = param[1], b = param[2], alpha = param[3],
                                   beta = param[4])

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]


           q = VaRES::varkumnormal(p, mu = param[1], sigma = param[2], a = param[3],
                                 b = param[4])

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


           q = VaRES::varkumweibull(p, a = param[1], b = param[2], alpha = param[3],
                                  sigma = param[4])

         } ,



         "laplace" = {     ## [R, R+]


           q = extraDistr::qlaplace(p, mu = param[1], sigma = param[2])

         } ,


         "levy" = {     ## [R, R+]


           q = rmutil::qlevy(p, m = param[1], s = param[2])

         } ,


         "linearfailurerate" = {     ## [R+, R+]


           q = VaRES::varlfr(p, a = param[1], b = param[2])

         } ,


         # "lindley" = {     ## [R+]
         #
         #   q = VGAM::qlind(p, theta = param[1])
         #
         # } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


           q = VaRES::varLNbeta(p, lambda = param[1], a = param[2], b = param[3])

         } ,


         "logcauchy" = {     ## [R, R+]


           q = VaRES::varlogcauchy(p, mu = param[1], sigma = param[2])

         } ,



         "loggamma" = {     ## [R, R+, R+]


           q = VGAM::qlgamma(p, location = param[1], scale = param[2], shape = param[3])

         } ,


         "loggumbel" = {     ## [R, R+]


           q = ssdtools::qlgumbel(p, llocation = param[1], lscale = param[2])

         } ,


         "loglaplace" = {     ## [R, R+, R+]


           q = VGAM::qloglap(p, location.ald = param[1], scale.ald = param[2], kappa = param[3])

         } ,


         "loglog" = {     ## [R+, >1]


           q = VaRES::varloglog(p, a = param[1], lambda = param[2])

         } ,


         "loglogistic" = {     ## [R+, R+]


           q = actuar::qllogis(p, shape = param[1], scale = param[2])

         } ,


         "lognormal" = {     ## [R, R+]


           q = stats::qlnorm(p, meanlog = param[1], sdlog = param[2])

         } ,


         "lognormal3" = {     ## [R, R+, R]


           q = EnvStats::qlnorm3(p, meanlog = param[1], sdlog = param[2], threshold = param[3])

         } ,


         "logistic" = {     ## [R, R+]


           q = stats::qlogis(p, location = param[1], scale = param[2])

         } ,


         "logisticexponential" = {     ## [R+, R+]


           q = VaRES::varlogisexp(p, lambda = param[1], a = param[2])

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


           q = VaRES::varlogisrayleigh(p, a = param[1], lambda = param[2])

         } ,


         "logseries" = {     ## [01]

           q = extraDistr::qlgser(p, theta = param[1])

         } ,



         "lomax" = {     ## [R+, R+]


           q = VGAM::qlomax(p, scale = param[1], shape3.q = param[2])

         } ,


         "makeham" = {     ## [R+, R+, R+]


           q = VGAM::qmakeham(p, scale = param[1], shape = param[2], epsilon = param[3])

         } ,


         "maxwell" = {     ## [R+]


           q = VGAM::qmaxwell(p, rate = param[1])

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


           q = VaRES::varMlaplace(p, theta = param[1], phi = param[2], psi = param[3])

         } ,


         "moexponential" = {     ## [R+, R+]


           q = VaRES::varmoexp(p, lambda = param[1], a = param[2])

         } ,


         "moweibull" = {     ## [R+, R+, R+]


           q = VaRES::varmoweibull(p, a = param[1], b = param[2], lambda = param[3])

         } ,


         "nakagami" = {     ## [R+, R+]


           q = VGAM::qnaka(p, scale = param[1], shape = param[2])

         } ,


         "ncchisquared" = {     ## [R+, R+]


           q = stats::qchisq(p, dq = param[1], ncp = param[2])

         } ,


         "ncF" = {     ## [R+, R+, R+]


           q = stats::qf(p, df1 = param[1], df2 = param[2], ncp = param[3])

         } ,


         "negativebinomial" = {     ## [N+, 01]

           q = stats::qnbinom(p, size = size, prob = param[1])

         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]


           q = GeneralizedHyperbolic::qnig(p, mu = param[1], delta = param[2],
                                           alpha = param[3], beta = param[4])

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


           q = extraDistr::qnsbeta(p, shape1 = param[1], shape2 = param[2],
                                   min = param[3], max = param[4])

         } ,


         "paralogistic" = {     ## [R+, R+]


           q = VGAM::qparalogistic(p, scale = param[1], shape1.a = param[2])

         } ,


         "pareto" = {     ## [R+, R+]


           q = extraDistr::qpareto(p, a = param[1], b = param[2])

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


           q = VaRES::varparetostable(p, lambda = param[1], nu = param[2], sigma = param[3])

         } ,


         "pareto1" = {     ## [R+, R+]


           q = VGAM::qparetoI(p, scale = param[1], shape = param[2])

         } ,


         "pareto2" = {     ## [R, R+, R+]


           q = VGAM::qparetoII(p, location = param[1], scale = param[2], shape = param[3])

         } ,


         "pareto3" = {     ## [R, R+, R+]


           q = VGAM::qparetoIII(p, location = param[1], scale = param[2], inequality = param[3])

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


           q = VGAM::qparetoIV(p, location = param[1], scale = param[2], inequality = param[3], shape = param[4])

         } ,


         "perks" = {     ## [R+, R+]


           q = VGAM::qperks(p, scale = param[1], shape = param[2])

         } ,


         "pctalaplace" = {     ## [R+, R]


           q = VaRES::varPCTAlaplace(p, a = param[1], theta = param[2])

         } ,


         "poisson" = {     ## [R+]

           q = stats::qpois(p, lambda = param[1])

         } ,


         "power1" = {     ## [R+]


           q = VaRES::varpower1(p, a = param[1])

         } ,


         "power2" = {     ## [R+]


           q = VaRES::varpower2(p, b = param[1])

         } ,


         "powerdistribution" = {     ## [R+, R+]


           q = extraDistr::qpower(p, alpha = param[1], beta = param[2])

         } ,


         "powerexponential" = {     ## [R, R+, R+]


           q = rmutil::qpowexp(p, m = param[1], s = param[2], q = param[3])

         } ,


         "rayleigh" = {     ## [R+]


           q = VGAM::qrayleigh(p, scale = param[1])

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


           q = VaRES::varrgamma(p, a = param[1], theta = param[2], phi = param[3])

         } ,



         "rice" = {     ## [R+, R+]


           q = VGAM::qrice(p, sigma = param[1], vee = param[2])

         } ,


         "scaledchisquared" = {     ## [R+, R+]


           q = extraDistr::qinvchisq(p, nu = param[1], tau = param[2])

         } ,


         "schabe" = {     ## [R+, R+]


           q = VaRES::varschabe(p, gamma = param[1], theta = param[2])

         } ,



         "simplex" = {     ## [01, R+]


           q = rmutil::qsimplex(p, m = param[1], s = param[2])

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


           q = rmutil::qskewlaplace(p, m = param[1], s = param[2], q = param[3])

         } ,


         "skewedt" = {     ## [R+, R+]


           q = skewt::qskt(p, dq = param[1], gamma = param[2])

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

           f = sn::qst(p, xi = param[1], omega = param[2], alpha = param[3], nu = param[4])

         } ,


         "skewednormal" = {     ## [R, R+, R, R]

           f = sn::qsn(p, xi = param[1], omega = param[2], alpha = param[3])

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           f = sgt::qsgt(p, mu = param[1], sigma = param[2], lambda = param[3], p = param[4])

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           f = gamlss.dist::qSEP(p, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,


         # "slash" = {     ## [R, R+]
         #
         #   q = extraDistr::qslash(p, mu = param[1], sigma = param[2])
         #
         # } ,


         "stable" = {     ## [02, -1+1, R+, R]


           q = stabledist::qstable(p, alpha = param[1], beta = param[2], gamma = param[3],
                                   delta = param[4])

         } ,


         "stacy" = {     ## [R+, R+, R+]


           q = VaRES::varstacygamma(p, gamma = param[1], c = param[2], theta = param[3])

         } ,


         "t" = {     ## [R, R+, R+]


           q =  param[1] + param[2]*stats::qt(p, df = param[3])

         } ,



         "tobit" = {     ## [R, R+]


           q = VGAM::qtobit(p, mean = param[1], sd = param[2])

         } ,


         "topple" = {     ## [01]


           q = VGAM::qtopple(p, shape = param[1])

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]


           q = actuar::qtrbeta(p, shape1 = param[1], shape2 = param[2], shape3 = param[3],
                               scale = param[4])

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]


           q = actuar::qtrgamma(p, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


           q = extraDistr::qtnorm(p, mean = param[1], sd = param[2], a = param[3],
                                  b = param[4])

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


           q = VGAM::qtruncpareto(p, lower = param[1], upper = param[2], shape = param[3])

         } ,


         "twosidedpower" = {     ## [01, R+]


           q = rmutil::qtwosidedpower(p, m = param[1], s =  param[2])

         } ,


         # "wald" = {     ## [R, R+]
         #
         #   q = extraDistr::qwald(p, mu = param[1], lambda = param[2])
         #
         # } ,


         "weibull" = {     ## [R+, R+]


           q = stats::qweibull(p, shape = param[1], scale = param[2])

         } ,



         "xie" = {     ## [R+, R+, R+]


           q = VaRES::varxie(p, a = param[1], b = param[2], lambda = param[3])

         } ,


         "yules" = {     ## [R+] mais > 0.5

           q = VGAM::qyules(p, shape = param[1])

         } ,


  )


  return(q)


}
