#'@title Probability density function
#'
#'@description This function computes the probability density function (pdf) of a univariate distribution
#'
#'@param family   distribution name; run the function distributions() for help
#'@param y  observations
#'@param param  parameters of the distribution; (1 x p)
#'
#'@return \item{f}{pdf}
#'
#'@export


PDF<-function(family,y,param){


  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]


             f = VaRES::daep(y, q1 = param[1], q2 = param[2], alpha = param[3])

         } ,


         "asymlaplace" = {    ## [R, R+, R+]


             f = VGAM::dalap(y, location = param[1], scale = param[2], kappa = param[3])

         } ,


         "asympower" = {    ## [01, R+, R+]


             f = VaRES::dasypower(y, a = param[1], lambda = param[2], delta = param[3])

         } ,


         "asymt" = {    ## [R+, R+, 01, R]

             f = VaRES::dast(y-param[4], nu1 = param[1], nu2 = param[2], alpha = param[3])

         } ,


         "beard" = {    ## [R+, R+, R+]


             f = VaRES::dbeard(y, a = param[1], b = param[2], rho = param[3])

         } ,


         "benini" = {     ## [R, R+]


             f = VGAM::dbenini(y, y0 = param[1], shape = param[2])

         } ,


         "benford" = {     ## [1 ou 2]


             f = VGAM::dbenf(y, ndigits = param[1])

         } ,


         "bernoulli" = {     ## [01]

             f = extraDistr::dbern(y, prob = param[1])

         } ,


         "beta" = {     ## [R+, R+]

             f = stats::dbeta(y, shape1 = param[1], shape2 = param[2])

         } ,



         "betabinomial" = {     ## [N+, R+, R+]


             f = extraDistr::dbbinom(y, size = param[1], alpha = param[2], beta = param[3])

         } ,


         "betageometric" = {     ## [R+, R+]


             f = VGAM::dbetageom(y, shape1 = param[1], shape2 = param[2])

         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]


             f = extraDistr::dbnbinom(y, size = param[1], alpha = param[2], beta = param[3])

         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetaburr(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetaburr7(y, a = param[1], b = param[2], c = param[3], k = param[4])

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


             f = VaRES::dbetaexp(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetafrechet(y, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetagompertz(y, b = param[1], c = param[2], d = param[3], eta = param[4])

         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]


             f = VaRES::dbetagumbel(y, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetagumbel2(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]


             f = VaRES::dbetalognorm(y, a = param[1], b = param[2], mu = param[3], sigma = param[4])

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetalomax(y, a = param[1], b = param[2], alpha = param[3], lambda = param[4])

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


             f = VGAM::dbetanorm(y, shape1 = param[1], shape2 = param[2],
                                       mean = param[3], sd = param[4])

         } ,


         "betaprime" = {     ## [R+, R+, R+]


             f = extraDistr::dbetapr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


             f = VaRES::dbetaweibull(y, a = param[1], b = param[2], alpha = param[3], sigma = param[4])

         } ,


         "bhattacharjee" = {     ## [R, R+, R+]


             f = extraDistr::dbhatt(y, mu = param[1], sigma = param[2], a = param[3])

         } ,


         "binomial" = {     ## [N+, 01]

             f = stats::dbinom(y, size = param[1], prob = param[2])

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]


             f = extraDistr::dfatigue(y, alpha = param[1], beta = param[2], mu = param[3])

         } ,


         "boxcox" = {     ## [R+, R+, R+]


             f = rmutil::dboxcox(y, m = param[1], s = param[2], f = param[3])

         } ,


         "burr" = {     ## [R+, R+, R+]


             f = actuar::dburr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "burr2param" = {     ## [R+, R+]


             f = VaRES::dburr(y, a = param[1], b = param[2])

         } ,


         "cauchy" = {     ## [R, R+]


             f = stats::dcauchy(y, location = param[1], scale = param[2])

         } ,


         "chen" = {     ## [R+, R+]


             f = VaRES::dchen(y, b = param[1], lambda = param[2])

         } ,


         "chi" = {     ## [R+]


             f = EnvStats::dchi(y, df = param[1])

         } ,


         "chisquared" = {     ## [R+]


             f = stats::dchisq(y, df = param[1])

         } ,


         "clg" = {     ## [R+, R+, R]


             f = VaRES::dclg(y, a = param[1], b = param[2], param[3])

         } ,


         "complementarybeta" = {     ## [R+, R+]


             f = VaRES::dcompbeta(y, a = param[1], b = param[2])

         } ,



         "dagum" = {     ## [R+, R+, R+]


             f = VGAM::ddagum(y, scale = param[1], shape1.a = param[2], shape2.p = param[3])

         } ,


         "diffzeta" = {     ## [R+, >1]


             f = VGAM::ddiffzeta(y, shape = param[1], start = param[2])

         } ,


         "discretegamma" = {     ## [R+, R+]


             f = extraDistr::ddgamma(y, shape = param[1], scale = param[2])

         } ,


         "discretelaplace" = {     ## [R, 01]


             f = extraDistr::ddlaplace(y, location = param[1], scale = param[2])

         } ,


         "discretenormal" = {     ## [R, R+]


             f = extraDistr::ddnorm(y, mean = param[1], sd = param[2])

         } ,


         "discreteweibull" = {     ## [01, R+]


             f = extraDistr::ddweibull(y, shape1 = param[1], shape2 = param[2])

         } ,



         "doubleweibull" = {     ## [R+, R, R+]

             f = VaRES::ddweibull(y, c = param[1], mu = param[2], sigma = param[3])

         } ,


         "ev" = {

                 ## [R, R+]
             f = VGAM::dgev(y, location = param[1], scale = param[2], shape = 0)

         } ,


         "exponential" = {     ## [R+]


             f = stats::dexp(y, rate = param[1])

         } ,


         "exponentialextension" = {     ## [R+, R+]


             f = VaRES::dexpext(y, lambda = param[1], a = param[2])

         } ,



         "exponentialgeometric" = {     ## [R+, 01]


             f = VGAM::dexpgeom(y, scale = param[1], shape = param[2])

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


             f = VGAM::dexplog(y, scale = param[1], shape = param[2])

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


             f = VaRES::dexppois(y, b = param[1], lambda = param[2])

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


             f = VaRES::dexppower(y, mu = param[1], sigma = param[2], a = param[3])

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


             f = VaRES::dexpexp(y, lambda = param[1], a = param[2])

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


             f = VaRES::dexplogis(y, a = param[1], b = param[2])

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


             f = VaRES::dexpweibull(y, a = param[1], alpha = param[2], sigma = param[3])

         } ,



         "F" = {     ## [R+, R+]


             f = stats::df(y, df1 = param[1], df2 = param[2])

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]


             f = actuar::dfpareto(y, min = param[1], shape1 = param[2],
                                        shape2 = param[3], shape3 = param[4],
                                        scale = param[5])

         } ,


         "fisk" = {     ## [R+, R+]

             f = VGAM::dfisk(y, scale = param[1], shape1.a = param[2])

         } ,


         "foldednormal" = {     ## [R, R+]


             f = VGAM::dfoldnorm(y, mean = param[1], sd = param[2])

         } ,


         "frechet" = {     ## [R+, R, R+]


             f = VGAM::dfrechet(y, shape = param[1], location = param[2], scale = param[3])

         } ,


         "gamma" = {     ## [R+, R+]


             f = stats::dgamma(y, shape =  param[1], scale = param[2])

         } ,


         "gammapoisson" = {     ## [R+, R+]

             f = extraDistr::dgpois(y, shape = param[1], scale = param[2])

         } ,


         "gaussian" = {     ## [R, R+]


             f = stats::dnorm(y, mean = param[1], sd = param[2])

         } ,


         "gev" = {     ## [R, R+, R]


             f = VGAM::dgev(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


             f = VaRES::dgeninvbeta(y, a = param[1], c = param[2], d = param[3])

         } ,


         "genlogis" = {     ## [R+, R, R+]


             f = VaRES::dgenlogis(y, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis3" = {     ## [R+, R, R+]


             f = VaRES::dgenlogis3(y, a = param[1], mu = param[2], sigma = param[3])

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


             f = VaRES::dgenlogis4(y, a = param[1], alpha = param[2], mu = param[3], sigma = param[4])

         } ,


         "genpowerweibull" = {     ## [R+, R+]


             f = VaRES::dgenpowerweibull(y, a = param[1], theta = param[2])

         } ,



         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]   [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

           f = GeneralizedHyperbolic::dghyp(y, mu = param[1], delta = param[2],
                                            alpha = param[3], beta = param[4],
                                            lambda = param[5])

         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           f = GLDEX::dgl(y, lambda1 = param[1], lambda2 = param[2], lambda3 = param[3], lambda4 = param[4])

         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

            f = gamlss.dist::dGT(y, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,


         "geometric" = {     ## [01]


             f = stats::dgeom(y, prob = param[1])

         } ,


         "gompertz" = {     ## [R+, R+]


             f = ssdtools::dgompertz(y, lscale = param[1], lshape = param[2])

         } ,


         "gpd" = {     ## [R, R+, R]


             f = VGAM::dgpd(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "gumbel" = {     ## [R, R+]


             f = VGAM::dgumbel(y, location = param[1], scale = param[2])

         } ,


         "gumbel2" = {     ## [R+, R+]


             f = VGAM::dgumbelII(y, scale = param[1], shape = param[2])

         } ,


         "halfcauchy" = {     ## [R+]


             f = extraDistr::dhcauchy(y, sigma = param[1])

         } ,


         "halflogistic" = {     ## [R+]


             f = VaRES::dhalflogis(y, lambda = param[1])

         } ,


         "halfnormal" = {     ## [R+]


             f = extraDistr::dhnorm(y, sigma = param[1])

         } ,


         "halft" = {     ## [R+, R+]


             f = extraDistr::dht(y, nu = param[1], sigma = param[2])

         } ,


         "hjorth" = {     ## [R+, R+, R+]


             f = rmutil::dhjorth(y, m = param[1], s = param[2], f = param[3])

         } ,


         "hblaplace" = {     ## [01, R, R+]


             f = VaRES::dHBlaplace(y, a = param[1], theta = param[2], phi = param[3])

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]


             f = GeneralizedHyperbolic::dhyperb(y, mu = param[1], delta = param[2],
                                                alpha = param[3], beta = param[4])

         } ,


         "huber" = {     ## [R, R+]


             f = extraDistr::dhuber(y, mu = param[1], sigma = param[2])

         } ,


         "hzeta" = {     ## [R+]


             f = VGAM::dhzeta(y, shape = param[1])

         } ,


         "inversebeta" = {     ## [R+, R+]


             f = VaRES::dinvbeta(y, a = param[1], b = param[2])

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


             f = actuar::dinvburr(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inversechisquared" = {     ## [R+]


             f = extraDistr::dinvchisq(y, nu = param[1])

         } ,


         "inverseexponential" = {     ## [R+]


             f = actuar::dinvexp(y, scale = param[1])

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


             f = VaRES::dinvexpexp(y, lambda = param[1], a = param[2])

         } ,


         "inversegamma" = {     ## [R+, R+]


             f = extraDistr::dinvgamma(y, alpha = param[1], beta = param[2])

         } ,


         "inverselomax" = {     ## [R+, R+]


             f = VGAM::dinv.lomax(y, scale = param[1], shape2.p = param[2])

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


             f = actuar::dinvparalogis(y, shape = param[1], scale = param[2])

         } ,


         "inversepareto" = {     ## [R+, R+]


             f = actuar::dinvpareto(y, shape = param[1], scale = param[2])

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


             f = actuar::dinvtrgamma(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "inverseweibull" = {     ## [R+, R+]


             f = actuar::dinvweibull(y, shape = param[1], scale = param[2])

         } ,


         "kumaraswamy" = {     ## [R+, R+]


             f = VGAM::dkumar(y, shape1 = param[1], shape2 = param[2])

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


             f = VaRES::dkumexp(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


             f = VaRES::dkumgamma(y, a = param[1], b = param[2], c = param[3], d = param[4])

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


             f = VaRES::dkumgumbel(y, a = param[1], b = param[2], mu = param[3],
                                           sigma = param[4])

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


             f = VaRES::dkumhalfnorm(y, sigma = param[1], a = param[2], b = param[3])

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]


             f = VaRES::dkumloglogis(y, a = param[1], b = param[2], alpha = param[3],
                                             beta = param[4])

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]


             f = VaRES::dkumnormal(y, mu = param[1], sigma = param[2], a = param[3],
                                           b = param[4])

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


             f = VaRES::dkumweibull(y, a = param[1], b = param[2], alpha = param[3],
                                            sigma = param[4])

         } ,



         "laplace" = {     ## [R, R+]


             f = extraDistr::dlaplace(y, mu = param[1], sigma = param[2])

         } ,


         "levy" = {     ## [R, R+]


             f = rmutil::dlevy(y, m = param[1], s = param[2])

         } ,


         "linearfailurerate" = {     ## [R+, R+]


             f = VaRES::dlfr(y, a = param[1], b = param[2])

         } ,


         "lindley" = {     ## [R+]


             f = VGAM::dlind(y, theta = param[1])

         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


             f = VaRES::dLNbeta(y, lambda = param[1], a = param[2], b = param[3])

         } ,


         "logcauchy" = {     ## [R, R+]


             f = VaRES::dlogcauchy(y, mu = param[1], sigma = param[2])

         } ,



         "loggamma" = {     ## [R, R+, R+]


             f = VGAM::dlgamma(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "loggumbel" = {     ## [R, R+]


             f = ssdtools::dlgumbel(y, llocation = param[1], lscale = param[2])

         } ,


         "loglaplace" = {     ## [R, R+, R+]


             f = VGAM::dloglap(y, location.ald = param[1], scale.ald = param[2], kappa = param[3])

         } ,


         "loglog" = {     ## [R+, >1]


             f = VaRES::dloglog(y, a = param[1], lambda = param[2])

         } ,


         "loglogistic" = {     ## [R+, R+]


             f = actuar::dllogis(y, shape = param[1], scale = param[2])

         } ,


         "lognormal" = {     ## [R, R+]


             f = stats::dlnorm(y, meanlog = param[1], sdlog = param[2])

         } ,


         "lognormal3" = {     ## [R, R+, R]


             f = EnvStats::dlnorm3(y, meanlog = param[1], sdlog = param[2], threshold = param[3])

         } ,


         "logistic" = {     ## [R, R+]


             f = stats::dlogis(y, location = param[1], scale = param[2])

         } ,


         "logisticexponential" = {     ## [R+, R+]


             f = VaRES::dlogisexp(y, lambda = param[1], a = param[2])

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


             f = VaRES::dlogisrayleigh(y, a = param[1], lambda = param[2])

         } ,


         "logseries" = {     ## [01]


             f = extraDistr::dlgser(y, theta = param[1])

         } ,



         "lomax" = {     ## [R+, R+]


             f = VGAM::dlomax(y, scale = param[1], shape3.q = param[2])

         } ,


         "makeham" = {     ## [R+, R+, R+]


             f = VGAM::dmakeham(y, scale = param[1], shape = param[2], epsilon = param[3])

         } ,


         "maxwell" = {     ## [R+]


             f = VGAM::dmaxwell(y, rate = param[1])

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


             f = VaRES::dMlaplace(y, theta = param[1], phi = param[2], psi = param[3])

         } ,



         "moexponential" = {     ## [R+, R+]


             f = VaRES::dmoexp(y, lambda = param[1], a = param[2])

         } ,


         "moweibull" = {     ## [R+, R+, R+]


             f = VaRES::dmoweibull(y, a = param[1], b = param[2], lambda = param[3])

         } ,


         "nakagami" = {     ## [R+, R+]


             f = VGAM::dnaka(y, scale = param[1], shape = param[2])

         } ,


         "ncchisquared" = {     ## [R+, R+]


             f = stats::dchisq(y, df = param[1], ncp = param[2])

         } ,


         "ncF" = {     ## [R+, R+, R+]


             f = stats::df(y, df1 = param[1], df2 = param[2], ncp = param[3])

         } ,


         "negativebinomial" = {     ## [N+, 01]

             f = stats::dnbinom(y, size = param[1], prob = param[2])

         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]

           f = GeneralizedHyperbolic::dnig(y, mu = param[1], delta = param[2],
                                           alpha = param[3], beta = param[4])

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


             f = extraDistr::dnsbeta(y, shape1 = param[1], shape2 = param[2],
                                           min = param[3], max = param[4])

         } ,



         "paralogistic" = {     ## [R+, R+]


             f = VGAM::dparalogistic(y, scale = param[1], shape1.a = param[2])

         } ,


         "pareto" = {     ## [R+, R+]


             f = extraDistr::dpareto(y, a = param[1], b = min(y))

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


             f = VaRES::dparetostable(y, lambda = param[1], nu = param[2], sigma = param[3])

         } ,


         "pareto1" = {     ## [R+, R+]


             f = VGAM::dparetoI(y, scale = param[1], shape = param[2])

         } ,


         "pareto2" = {     ## [R, R+, R+]


             f = VGAM::dparetoII(y, location = param[1], scale = param[2], shape = param[3])

         } ,


         "pareto3" = {     ## [R, R+, R+]


             f = VGAM::dparetoIII(y, location = param[1], scale = param[2], inequality = param[3])

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


             f = VGAM::dparetoIV(y, location = param[1], scale = param[2], inequality = param[3], shape = param[4])

         } ,


         "perks" = {     ## [R+, R+]


             f = VGAM::dperks(y, scale = param[1], shape = param[2])

         } ,


         "pctalaplace" = {     ## [R+, R]


             f = VaRES::dPCTAlaplace(y, a = param[1], theta = param[2])

         } ,


         "poisson" = {     ## [R+]


             f = stats::dpois(y, lambda = param[1])

         } ,


         "power1" = {     ## [R+]


             f = VaRES::dpower1(y, a = param[1])

         } ,


         "power2" = {     ## [R+]


             f = VaRES::dpower2(y, b = param[1])

         } ,


         "powerdistribution" = {     ## [R+, R+]


             f = extraDistr::dpower(y, alpha = param[1], beta = param[2])

         } ,


         "powerexponential" = {     ## [R, R+, R+]


             f = rmutil::dpowexp(y, m = param[1], s = param[2], f = param[3])

         } ,


         "rayleigh" = {     ## [R+]


             f = VGAM::drayleigh(y, scale = param[1])

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


             f = VaRES::drgamma(y, a = param[1], theta = param[2], phi = param[3])

         } ,



         "rice" = {     ## [R+, R+]


             f = VGAM::drice(y, sigma = param[1], vee = param[2])

         } ,


         "scaledchisquared" = {     ## [R+, R+]


             f = extraDistr::dinvchisq(y, nu = param[1], tau = param[2])

         } ,


         "schabe" = {     ## [R+, R+]


             f = VaRES::dschabe(y, gamma = param[1], theta = param[2])

         } ,



         "simplex" = {     ## [01, R+]


             f = rmutil::dsimplex(y, m = param[1], s = param[2])

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


             f = rmutil::dskewlaplace(y, m = param[1], s = param[2], f = param[3])

         } ,


         "skewedt" = {     ## [R+, R+]


             f = skewt::dskt(y, df = param[1], gamma = param[2])

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

             f = sn::dst(y, xi = param[1], omega = param[2], alpha = param[3], nu = param[4])

         } ,


         "skewednormal" = {     ## [R, R+, R, R]

             f = sn::dsn(y, xi = param[1], omega = param[2], alpha = param[3])

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

            f = sgt::dsgt(y, mu = param[1], sigma = param[2], lambda = param[3], p = param[4])

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           f = gamlss.dist::dSEP(y, mu = param[1], sigma = param[2], nu = param[3], tau = param[4])

         } ,


         "slash" = {     ## [R, R+]


             f = extraDistr::dslash(y, mu = param[1], sigma = param[2])

         } ,


         "stable" = {     ## [02, -1+1, R+, R]


             f = stabledist::dstable(y, alpha = param[1], beta = param[2], gamma = param[3],
                                           delta = param[4])

         } ,


         "stacy" = {     ## [R+, R+, R+]


             f = VaRES::dstacygamma(y, gamma = param[1], c = param[2], theta = param[3])

         } ,


         "t" = {     ## [R, R+, R+]


             f =  stats::dt((y-param[1])/param[2], df = param[3])

         } ,



         "tobit" = {     ## [R, R+]


             f = VGAM::dtobit(y, mean = param[1], sd = param[2])

         } ,


         "topple" = {     ## [01]


             f = VGAM::dtopple(y, shape = param[1])

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]


             f = actuar::dtrbeta(y, shape1 = param[1], shape2 = param[2], shape3 = param[3],
                                       scale = param[4])

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]


             f = actuar::dtrgamma(y, shape1 = param[1], shape2 = param[2], scale = param[3])

         } ,


         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


             f = extraDistr::dtnorm(y, mean = param[1], sd = param[2], a = param[3],
                                          b = param[4])

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


             f = VGAM::dtruncpareto(y, lower = param[1], upper = param[2], shape = param[3])

         } ,


         "twosidedpower" = {     ## [01, R+]


             f = rmutil::dtwosidedpower(y, m = param[1], s =  param[2])

         } ,


         "wald" = {     ## [R, R+]


             f = extraDistr::dwald(y, mu = param[1], lambda = param[2])

         } ,


         "weibull" = {     ## [R+, R+]


             f = stats::dweibull(y, shape = param[1], scale = param[2])

         } ,



         "xie" = {     ## [R+, R+, R+]


             f = VaRES::dxie(y, a = param[1], b = param[2], lambda = param[3])

         } ,


         "yules" = {     ## [R+] mais > 0.5


             f = VGAM::dyules(y, shape = param[1])

         } ,


  )

  f[is.nan(f)] <- 0
  return(f)


}





