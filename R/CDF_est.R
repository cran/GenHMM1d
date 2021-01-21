#'@title Cumulative density function
#'
#'@description This function compute the cumulative density function of an univariate distribution (used) in the CVM calculation
#'
#'@param family   distribution name; run the command distributions() for help
#'@param y  observations
#'@param param  parameters of the distribution
#'@param size additional parameter for some discrete distributions
#'@param u_Ros used internally for the computation of the Rosenblatt transform for univariate discrete distribution
#'
#'@return \item{f}{cdf}
#'
#'
#'@export
#'@keywords internal


CDF_est<-function(family,y,param,size=0,u_Ros=0){



  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]

           f = VaRES::paep(y, q1 = param[1], q2 = param[2], alpha = param[3])

         } ,


         "asymlaplace" = {    ## [R, R+, R+]

           f = VGAM::palap(y, location = param[1], scale = param[2], kappa = param[3])

         } ,


         "asympower" = {    ## [01, R+, R+]

           f = VaRES::pasypower(y, a = param[1], lambda = param[2], delta = param[3])

         } ,


         "asymt" = {    ## [R+, R+, 01, R]

           f = VaRES::past(y-param[4], nu1 = param[1], nu2 = param[2], alpha = param[3])

         } ,


         "beard" = {    ## [R+, R+, R+]

           f = VaRES::pbeard(y, a = param[1], b = param[2], rho = param[3])

         } ,


         "benini" = {     ## [R, R+]


           f = VGAM::pbenini(y, y0 = param[1], shape = param[2])

         } ,


         "benford" = {     ## [1 ou 2]


           f = VGAM::pbenf(y, ndigits = param[1])

         } ,


         "bernoulli" = {     ## [01]


           f = extraDistr::pbern(y, prob = param[1])

         } ,


         "beta" = {     ## [R+, R+]


           f = stats::pbeta(y, shape1 = param[1], shape2 = param[2])

         } ,


         "betabinomial" = {     ## [N+, R+, R+]

           f_1 = extraDistr::pbbinom(y-0.7, size = size, alpha = param[1], beta = param[2])
           f_2 = extraDistr::pbbinom(y, size = size, alpha = param[1], beta = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "betageometric" = {     ## [R+, R+]

           f_1 = VGAM::pbetageom(y-0.7, shape1 = param[1], shape2 = param[2])
           f_2 = VGAM::pbetageom(y, shape1 = param[1], shape2 = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]

           f_1 = extraDistr::pbnbinom(y-0.7, size = size, alpha = param[1], beta = param[2])
           f_2 = extraDistr::pbnbinom(y, size = size, alpha = param[1], beta = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetaburr(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetaburr7(y, a = param[1], b = param[2], c = param[3], k = param[4])

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


           f = VaRES::pbetaexp(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetafrechet(y, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetagompertz(y, b = param[1], c = param[2], d = param[3], eta = param[4])

         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]


           f = VaRES::pbetagumbel(y, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetagumbel2(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]


           f = VaRES::pbetalognorm(y, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetalomax(y, a = param[1], b = param[2], alpha = param[3], lambda = param[4])

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


           f = VGAM::pbetanorm(y, shape1 = param[1], shape2 = param[2],
                               mean = param[3], sd = param[4])

         } ,


         "betaprime" = {     ## [R+, R+, R+]


           f = extraDistr::pbetapr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


           f = VaRES::pbetaweibull(y, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         "bhattacharjee" = {     ## [R, R+, R+]


           f = extraDistr::pbhatt(y, mu = param[1], sigma = param[2], a = param[3])

         } ,


         "binomial" = {     ## [N+, 01]

           f_1 = stats::pbinom(y-0.7, size = size, prob = param[1])
           f_2 = stats::pbinom(y, size = size, prob = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]


           f = extraDistr::pfatigue(y, alpha = param[1], beta = param[2], mu = param[3])

         } ,


         "boxcox" = {     ## [R+, R+, R+]


           f = rmutil::pboxcox(y, m = param[1], s = param[2], f = param[3])

         } ,


         "burr" = {     ## [R+, R+, R+]


           f = actuar::pburr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "burr2param" = {     ## [R+, R+]


           f = VaRES::pburr(y, a = param[1], b = param[2])

         } ,


         "cauchy" = {     ## [R, R+]


           f = stats::pcauchy(y, location = param[1], scale = param[2])

         } ,


         "chen" = {     ## [R+, R+]


           f = VaRES::pchen(y, b = param[1], lambda = param[2])

         } ,


         "chi" = {     ## [R+]


           f = EnvStats::pchi(y, df = param[1])

         } ,


         "chisquared" = {     ## [R+]


           f = stats::pchisq(y, df = param[1])

         } ,


         "clg" = {     ## [R+, R+, R]


           f = VaRES::pclg(y, a = param[1], b = param[2], param[3])

         } ,


         "complementarybeta" = {     ## [R+, R+]


           f = VaRES::pcompbeta(y, a = param[1], b = param[2])

         } ,



         "dagum" = {     ## [R+, R+, R+]


           f = VGAM::pdagum(y, scale = param[1], shape1.a = param[2], shape2.p = param[3])

         } ,


         "diffzeta" = {     ## [R+, >1]


           f = VGAM::pdiffzeta(y, shape = param[1], start = param[2])

         } ,


         "discretegamma" = {     ## [R+, R+]

           f_1 = extraDistr::pdgamma(y-0.7, shape = param[1], scale = param[2])
           f_2 = extraDistr::pdgamma(y, shape = param[1], scale = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "discretelaplace" = {     ## [R, 01]

           f_1 = extraDistr::pdlaplace(y-0.7, location = param[1], scale = param[2])
           f_2 = extraDistr::pdlaplace(y, location = param[1], scale = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "discretenormal" = {     ## [R, R+]

           f_1 = extraDistr::pdnorm(y-0.7, mean = param[1], sd = param[2])
           f_2 = extraDistr::pdnorm(y, mean = param[1], sd = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "discreteweibull" = {     ## [01, R+]

           f_1 = extraDistr::pdweibull(y-0.7, shape1 = param[1], shape2 = param[2])
           f_2 = extraDistr::pdweibull(y, shape1 = param[1], shape2 = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,



         "doubleweibull" = {     ## [R+, R, R+]

           f = VaRES::pdweibull(y, c = param[1], mu = param[2], sigma = param[3])

         } ,


         "ev" = {

           ## [R, R+]
           f = VGAM::pgev(y, location = param[1], scale = param[2], shape = 0)

         } ,


         "exponential" = {     ## [R+]


           f = stats::pexp(y, rate = param[1])

         } ,


         "exponentialextension" = {     ## [R+, R+]


           f = VaRES::pexpext(y, lambda = param[1], a = param[2])

         } ,



         "exponentialgeometric" = {     ## [R+, 01]


           f = VGAM::pexpgeom(y, scale = param[1], shape = param[2])

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


           f = VGAM::pexplog(y, scale = param[1], shape = param[2])

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


           f = VaRES::pexppois(y, b = param[1], lambda = param[2])

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


           f = VaRES::pexppower(y, mu = param[1], sigma = param[2], a = param[3])

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


           f = VaRES::pexpexp(y, lambda = param[1], a = param[2])

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


           f = VaRES::pexplogis(y, a = param[1], b = param[2])

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


           f = VaRES::pexpweibull(y, a = param[1], alpha = param[2], sigma = param[3])

         } ,


         "F" = {     ## [R+, R+]


           f = stats::pf(y, df1 = param[1], df2 = param[2])

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]


           f = actuar::pfpareto(y, min = param[1], shape1 = param[2],
                                shape2 = param[3], shape3 = param[4],
                                scale = param[5])

         } ,


         "fisk" = {     ## [R+, R+]

           f = VGAM::pfisk(y, scale = param[1], shape1.a = param[2])

         } ,


         "foldednormal" = {     ## [R, R+]


           f = VGAM::pfoldnorm(y, mean = param[1], sd = param[2])

         } ,


         "frechet" = {     ## [R+, R, R+]


           f = VGAM::pfrechet(y, shape = param[1], location = param[2], scale = param[3])

         } ,


         "gamma" = {     ## [R+, R+]


           f = stats::pgamma(y, shape =  param[1], scale = param[2])

         } ,



         "gammapoisson" = {     ## [R+, R+]

           f_1 = extraDistr::pgpois(y-0.7, shape = param[1], scale = param[2])
           f_2 = extraDistr::pgpois(y, shape = param[1], scale = param[2])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "gaussian" = {     ## [R, R+]


           f = stats::pnorm(y, mean = param[1], sd = param[2])

         } ,


         "gev" = {     ## [R, R+, R]


           f = VGAM::pgev(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


           f = VaRES::pgeninvbeta(y, a = param[1], c = param[2], d = param[3])

         } ,


         "genlogis" = {     ## [R+, R, R+]


           f = VaRES::pgenlogis(y, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis3" = {     ## [R+, R, R+]


           f = VaRES::pgenlogis3(y, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


           f = VaRES::pgenlogis4(y, a = param[1], alpha = param[2], mu = param[3], sigma = param[4])

         } ,


         "genpowerweibull" = {     ## [R+, R+]


           f = VaRES::pgenpowerweibull(y, a = param[1], theta = param[2])

         } ,



         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

           f = GeneralizedHyperbolic::pghyp(y, mu = param[1], delta = param[2],
                                            alpha = param[3], beta = param[4],
                                            lambda = param[5])

         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           f = GLDEX::pgl(y, lambda1 = param[1], lambda2 = param[2], lambda3 = param[3], lambda4 = param[4])

         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

           f = gamlss.dist::pGT(y, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,



         "geometric" = {     ## [01]

           f_1 = stats::pgeom(y-0.7, prob = param[1])
           f_2 = stats::pgeom(y, prob = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "gompertz" = {     ## [R+, R+]


           f = ssdtools::pgompertz(y, lscale = param[1], lshape = param[2])

         } ,


         "gpd" = {     ## [R, R+, R]


           f = VGAM::pgpd(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "gumbel" = {     ## [R, R+]


           f = VGAM::pgumbel(y, location = param[1], scale = param[2])

         } ,


         "gumbel2" = {     ## [R+, R+]


           f = VGAM::pgumbelII(y, scale = param[1], shape = param[2])

         } ,


         "halfcauchy" = {     ## [R+]


           f = extraDistr::phcauchy(y, sigma = param[1])

         } ,


         "halflogistic" = {     ## [R+]


           f = VaRES::phalflogis(y, lambda = param[1])

         } ,


         "halfnormal" = {     ## [R+]


           f = extraDistr::phnorm(y, sigma = param[1])

         } ,


         "halft" = {     ## [R+, R+]


           f = extraDistr::pht(y, nu = param[1], sigma = param[2])

         } ,


         "hjorth" = {     ## [R+, R+, R+]


           f = rmutil::phjorth(y, m = param[1], s = param[2], f = param[3])

         } ,


         "hblaplace" = {     ## [01, R, R+]


           f = VaRES::pHBlaplace(y, a = param[1], theta = param[2], phi = param[3])

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]  [mu, delta, alpha, beta] (avec alpha^2 > beta^2)


           f = GeneralizedHyperbolic::phyperb(y, mu = param[1], delta = param[2],
                                              alpha = param[3], beta = param[4])

         } ,


         "huber" = {     ## [R, R+]


           f = extraDistr::phuber(y, mu = param[1], sigma = param[2])

         } ,


         "hzeta" = {     ## [R+]

           f_1 = VGAM::phzeta(y-0.7, shape = param[1])
           f_2 = VGAM::phzeta(y, shape = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "inversebeta" = {     ## [R+, R+]


           f = VaRES::pinvbeta(y, a = param[1], b = param[2])

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


           f = actuar::pinvburr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inversechisquared" = {     ## [R+]


           f = extraDistr::pinvchisq(y, nu = param[1])

         } ,


         "inverseexponential" = {     ## [R+]


           f = actuar::pinvexp(y, scale = param[1])

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


           f = VaRES::pinvexpexp(y, lambda = param[1], a = param[2])

         } ,


         "inversegamma" = {     ## [R+, R+]


           f = extraDistr::pinvgamma(y, alpha = param[1], beta = param[2])

         } ,


         "inverselomax" = {     ## [R+, R+]


           f = VGAM::pinv.lomax(y, scale = param[1], shape2.p = param[2])

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


           f = actuar::pinvparalogis(y, shape = param[1], scale = param[2])

         } ,


         "inversepareto" = {     ## [R+, R+]


           f = actuar::pinvpareto(y, shape = param[1], scale = param[2])

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


           f = actuar::pinvtrgamma(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inverseweibull" = {     ## [R+, R+]


           f = actuar::pinvweibull(y, shape = param[1], scale = param[2])

         } ,


         "kumaraswamy" = {     ## [R+, R+]


           f = VGAM::pkumar(y, shape1 = param[1], shape2 = param[2])

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


           f = VaRES::pkumexp(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


           f = VaRES::pkumgamma(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


           f = VaRES::pkumgumbel(y, a = param[1], b = param[2], mu = param[3],
                                 sigma = param[4])

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


           f = VaRES::pkumhalfnorm(y, sigma = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]


           f = VaRES::pkumloglogis(y, a = param[1], b = param[2], alpha = param[3],
                                   beta = param[4])

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]


           f = VaRES::pkumnormal(y, mu = param[1], sigma = param[2], a = param[3],
                                 b = param[4])

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


           f = VaRES::pkumweibull(y, a = param[1], b = param[2], alpha = param[3],
                                  sigma = param[4])

         } ,



         "laplace" = {     ## [R, R+]


           f = extraDistr::plaplace(y, mu = param[1], sigma = param[2])

         } ,


         "levy" = {     ## [R, R+]


           f = rmutil::plevy(y, m = param[1], s = param[2])

         } ,


         "linearfailurerate" = {     ## [R+, R+]


           f = VaRES::plfr(y, a = param[1], b = param[2])

         } ,


         "lindley" = {     ## [R+]


           f = VGAM::plind(y, theta = param[1])

         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


           f = VaRES::pLNbeta(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "logcauchy" = {     ## [R, R+]


           f = VaRES::plogcauchy(y, mu = param[1], sigma = param[2])

         } ,



         "loggamma" = {     ## [R, R+, R+]


           f = VGAM::plgamma(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "loggumbel" = {     ## [R, R+]


           f = ssdtools::plgumbel(y, llocation = param[1], lscale = param[2])

         } ,


         "loglaplace" = {     ## [R, R+, R+]


           f = VGAM::ploglap(y, location.ald = param[1], scale.ald = param[2], kappa = param[3])

         } ,


         "loglog" = {     ## [R+, >1]


           f = VaRES::ploglog(y, a = param[1], lambda = param[2])

         } ,


         "loglogistic" = {     ## [R+, R+]


           f = actuar::pllogis(y, shape = param[1], scale = param[2])

         } ,


         "lognormal" = {     ## [R, R+]


           f = stats::plnorm(y, meanlog = param[1], sdlog = param[2])

         } ,


         "lognormal3" = {     ## [R, R+, R]


           f = EnvStats::plnorm3(y, meanlog = param[1], sdlog = param[2], threshold = param[3])

         } ,


         "logistic" = {     ## [R, R+]


           f = stats::plogis(y, location = param[1], scale = param[2])

         } ,


         "logisticexponential" = {     ## [R+, R+]


           f = VaRES::plogisexp(y, lambda = param[1], a = param[2])

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


           f = VaRES::plogisrayleigh(y, a = param[1], lambda = param[2])

         } ,


         "logseries" = {     ## [01]

           f_1 = extraDistr::plgser(y-0.7, theta = param[1])
           f_2 = extraDistr::plgser(y, theta = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,



         "lomax" = {     ## [R+, R+]


           f = VGAM::plomax(y, scale = param[1], shape3.q = param[2])

         } ,


         "makeham" = {     ## [R+, R+, R+]


           f = VGAM::pmakeham(y, scale = param[1], shape = param[2], epsilon = param[3])

         } ,


         "maxwell" = {     ## [R+]


           f = VGAM::pmaxwell(y, rate = param[1])

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


           f = VaRES::pMlaplace(y, theta = param[1], phi = param[2], psi = param[3])

         } ,


         "moexponential" = {     ## [R+, R+]


           f = VaRES::pmoexp(y, lambda = param[1], a = param[2])

         } ,


         "moweibull" = {     ## [R+, R+, R+]


           f = VaRES::pmoweibull(y, a = param[1], b = param[2], lambda = param[3])

         } ,


         "nakagami" = {     ## [R+, R+]


           f = VGAM::pnaka(y, scale = param[1], shape = param[2])

         } ,


         "ncchisquared" = {     ## [R+, R+]


           f = stats::pchisq(y, df = param[1], ncp = param[2])

         } ,


         "ncF" = {     ## [R+, R+, R+]


           f = stats::pf(y, df1 = param[1], df2 = param[2], ncp = param[3])

         } ,


         "negativebinomial" = {     ## [N+, 01]

           f_1 = stats::pnbinom(y-0.7, size = size, prob = param[1])
           f_2 = stats::pnbinom(y, size = size, prob = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "normalinversegaussian" = {     ## [R, R, R+, R]


           f = GeneralizedHyperbolic::pnig(y, mu = param[1], delta = param[2],
                                           alpha = param[3], beta = param[4])

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


           f = extraDistr::pnsbeta(y, shape1 = param[1], shape2 = param[2],
                                   min = param[3], max = param[4])

         } ,


         "paralogistic" = {     ## [R+, R+]


           f = VGAM::pparalogistic(y, scale = param[1], shape1.a = param[2])

         } ,


         "pareto" = {     ## [R+, R+]


           f = extraDistr::ppareto(y, a = param[1], b = min(y))

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


           f = VaRES::pparetostable(y, lambda = param[1], nu = param[2], sigma = param[3])

         } ,


         "pareto1" = {     ## [R+, R+]


           f = VGAM::pparetoI(y, scale = param[1], shape = param[2])

         } ,


         "pareto2" = {     ## [R, R+, R+]


           f = VGAM::pparetoII(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "pareto3" = {     ## [R, R+, R+]


           f = VGAM::pparetoIII(y, location = param[1], scale = param[2], inequality = param[3])

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


           f = VGAM::pparetoIV(y, location = param[1], scale = param[2], inequality = param[3], shape = param[4])

         } ,


         "perks" = {     ## [R+, R+]


           f = VGAM::pperks(y, scale = param[1], shape = param[2])

         } ,


         "pctalaplace" = {     ## [R+, R]


           f = VaRES::pPCTAlaplace(y, a = param[1], theta = param[2])

         } ,


         "poisson" = {     ## [R+]

           f_1 = stats::ppois(y-0.7, lambda = param[1])
           f_2 = stats::ppois(y, lambda = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


         "power1" = {     ## [R+]


           f = VaRES::ppower1(y, a = param[1])

         } ,


         "power2" = {     ## [R+]


           f = VaRES::ppower2(y, b = param[1])

         } ,


         "powerdistribution" = {     ## [R+, R+]


           f = extraDistr::ppower(y, alpha = param[1], beta = param[2])

         } ,


         "powerexponential" = {     ## [R, R+, R+]


           f = rmutil::ppowexp(y, m = param[1], s = param[2], f = param[3])

         } ,


         "rayleigh" = {     ## [R+]


           f = VGAM::prayleigh(y, scale = param[1])

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


           f = VaRES::prgamma(y, a = param[1], theta = param[2], phi = param[3])

         } ,



         "rice" = {     ## [R+, R+]


           f = VGAM::price(y, sigma = param[1], vee = param[2])

         } ,


         "scaledchisquared" = {     ## [R+, R+]


           f = extraDistr::pinvchisq(y, nu = param[1], tau = param[2])

         } ,


         "schabe" = {     ## [R+, R+]


           f = VaRES::pschabe(y, gamma = param[1], theta = param[2])

         } ,



         "simplex" = {     ## [01, R+]


           f = rmutil::psimplex(y, m = param[1], s = param[2])

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


           f = rmutil::pskewlaplace(y, m = param[1], s = param[2], f = param[3])

         } ,


         "skewedt" = {     ## [R+, R+]


           f = skewt::pskt(y, df = param[1], gamma = param[2])

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

           f = sn::pst(y, xi = param[1], omega = param[2], alpha = param[3], nu = param[4])

         } ,


         "skewednormal" = {     ## [R, R+, R, R]

           f = sn::psn(y, xi = param[1], omega = param[2], alpha = param[3])

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           f = sgt::psgt(y, mu = param[1], sigma = param[2], lambda = param[3], p = param[4])

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           f = gamlss.dist::pSEP(y, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,


         "slash" = {     ## [R, R+]


           f = extraDistr::pslash(y, mu = param[1], sigma = param[2])

         } ,


         "stable" = {     ## [02, -1+1, R+, R]


           f = stabledist::pstable(y, alpha = param[1], beta = param[2], gamma = param[3],
                                   delta = param[4])

         } ,


         "stacy" = {     ## [R+, R+, R+]


           f = VaRES::pstacygamma(y, gamma = param[1], c = param[2], theta = param[3])

         } ,


         "t" = {     ## [R, R+, R+]


           f =  stats::pt((y-param[1])/param[2], df = param[3])

         } ,



         "tobit" = {     ## [R, R+]


           f = VGAM::ptobit(y, mean = param[1], sd = param[2])

         } ,


         "topple" = {     ## [01]


           f = VGAM::ptopple(y, shape = param[1])

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]


           f = actuar::ptrbeta(y, shape1 = param[1], shape2 = param[2], shape3 = param[3],
                               scale = param[4])

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]


           f = actuar::ptrgamma(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


           f = extraDistr::ptnorm(y, mean = param[1], sd = param[2], a = param[3],
                                  b = param[4])

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


           f = VGAM::ptruncpareto(y, lower = param[1], upper = param[2], shape = param[3])

         } ,


         "twosidedpower" = {     ## [01, R+]


           f = rmutil::ptwosidedpower(y, m = param[1], s =  param[2])

         } ,


         "wald" = {     ## [R, R+]


           f = extraDistr::pwald(y, mu = param[1], lambda = param[2])

         } ,


         "weibull" = {     ## [R+, R+]


           f = stats::pweibull(y, shape = param[1], scale = param[2])

         } ,


         "xie" = {     ## [R+, R+, R+]


           f = VaRES::pxie(y, a = param[1], b = param[2], lambda = param[3])

         } ,


         "yules" = {     ## [R+] mais > 0.5

           f_1 = VGAM::pyules(y, shape = param[1])
           f_2 = VGAM::pyules(y, shape = param[1])
           f = (1-u_Ros)*f_1 + u_Ros*f_2

         } ,


  )


  f[is.nan(f)] <- 0
  return(f)


}



pgeometricprog<-function(y,prob){

  f = matrix(0, nrow = length(y), ncol = 1)

  for (i in 1:length(y)){
    if (y[1]>=0){
      f[i] = 1 - (1-prob)^(y[i]+1)
    }
  }

  return(f)

}
