#'@title Probability density function
#'
#'@description This function computes the probability density function of a univariate distribution with unconstrained parameters
#'
#'@param family   distribution name; run the command distributions() for help
#'@param y  observations
#'@param param  parameters of the distribution (unconstrained,  in R)
#'@param size additional parameter for some discrete distributions
#'
#'@return \item{f}{pdf}
#'
#'@export
#'@keywords internal


PDF_unc<-function(family,y,param,size=0){


  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]


           f = VaRES::daep(y, q1 = exp(param[1]), q2 = exp(param[2]),
                           alpha = 1/(1+exp(-param[3])))+1e-10

         } ,





         "asympower" = {    ## [01, R+, R+]


           f = VaRES::dasypower(y, a = 1/(1+exp(-param[1])), lambda = exp(param[2]),
                                delta = exp(param[3]))+1e-10

         } ,


         "asymt" = {    ## [R+, R+, 01, R]


           f = VaRES::dast(y-param[4], nu1 = exp(param[1]), nu2 = exp(param[2]),
                           alpha = 1/(1+exp(-param[3])) )+1e-10

         } ,


         "beard" = {    ## [R+, R+, R+]


           f = VaRES::dbeard(y, a = exp(param[1]), b = exp(param[2]), rho = exp(param[3]))+1e-10

         } ,


         "benini" = {     ## [R, R+]


           f = VGAM::dbenini(y, y0 = param[1], shape = exp(param[2]))+1e-10

         } ,


         "benford" = {     ## [1 ou 2]

           if (exp(param[1])>2) {
             f = VGAM::dbenf(y, ndigits = 2)+1e-10
           } else if (exp(param[1])<1) {
             f = VGAM::dbenf(y, ndigits = 1)+1e-10
           } else {
             f = VGAM::dbenf(y, ndigits = round(param[1]))+1e-10
           }
         } ,


         "bernoulli" = {     ## [01]

           f = extraDistr::dbern(y, prob = 1/(1+exp(-param[1])) )+1e-10

         } ,


         "beta" = {     ## [R+, R+]


           f = stats::dbeta(y, shape1 = exp(param[1]), shape2 = exp(param[2]))+1e-10

         } ,


         "betabinomial" = {     ## [N+, R+, R+]


           f = extraDistr::dbbinom(y, size = size, alpha = exp(param[1]), beta = exp(param[2]))+1e-10

         } ,


         "betageometric" = {     ## [R+, R+]


           f = VGAM::dbetageom(y, shape1 = exp(param[1]), shape2 = exp(param[2]))+1e-10

         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]


           f = extraDistr::dbnbinom(y, size = size, alpha = exp(param[1]), beta = exp(param[2]))+1e-10

         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]


           f = VaRES::dbetaburr(y, a = exp(param[1]), b = exp(param[2]), c = exp(param[3]),
                                d = exp(param[4]))+1e-10

         } ,


         "betaburr7" = {     ## [R+, R+]


           f = VaRES::dbetaburr7(y, a = exp(param[1]), b = exp(param[2]), c = exp(param[3]),
                                 k = exp(param[4]))+1e-10

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


           f = VaRES::dbetaexp(y, lambda = exp(param[1]), a = exp(param[2]), b = exp(param[3]))+1e-10

         } ,


         "betafrechet" = {     ## [R+, R+]


           f = VaRES::dbetafrechet(y, a = exp(param[1]), b = exp(param[2]), alpha = exp(param[3]),
                                   sigma = exp(param[4]))+1e-10

         } ,


         "betagompertz" = {     ## [R+, R+]


           f = VaRES::dbetagompertz(y, b = exp(param[1]), c = exp(param[2]), d = exp(param[3]),
                                    eta = exp(param[4]))+1e-10

         } ,


         "betagumbel" = {     ## [R+, R+]


           f = VaRES::dbetagumbel(y, a = exp(param[1]), b = exp(param[2]), mu = param[3],
                                  sigma = exp(param[4]))+1e-10

         } ,


         "betagumbel2" = {     ## [R+, R+]


           f = VaRES::dbetagumbel2(y, a = exp(param[1]), b = exp(param[2]), c = exp(param[3]),
                                   d = exp(param[4]))+1e-10

         } ,


         "betalognormal" = {     ## [R+, R+]


           f = VaRES::dbetalognorm(y, a = exp(param[1]), b = exp(param[2]), mu = param[3],
                                   sigma = exp(param[4]))+1e-10

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


           f = VaRES::dbetalomax(y, a = exp(param[1]), b = exp(param[2]), alpha = exp(param[3]),
                                 lambda = exp(param[4]))+1e-10

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


           f = VGAM::dbetanorm(y, shape1 = exp(param[1]), shape2 = exp(param[2]),
                               mean = param[3], sd = exp(param[4]))+1e-10

         } ,


         "betaprime" = {     ## [R+, R+, R+]


           f = extraDistr::dbetapr(y, shape1 = exp(param[1]), shape2 = exp(param[2]), scale = exp(param[3]))+1e-10

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


           f = VaRES::dbetaweibull(y, a = exp(param[1]), b = exp(param[2]), alpha = exp(param[3]),
                                   sigma = exp(param[4]))+1e-10

         } ,


         "bhattacharjee" = {     ## [R, R+, R+]


           f = extraDistr::dbhatt(y, mu = param[1], sigma = exp(param[2]), a = exp(param[3]))+1e-10

         } ,


         "binomial" = {     ## [N+, 01]

           f = stats::dbinom(y, size = size, prob = 1/(1+exp(-param[1])) )+1e-10

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]


           f = extraDistr::dfatigue(y, alpha = exp(param[1]), beta = exp(param[2]), mu = param[3])+1e-10

         } ,


         "boxcox" = {     ## [R+, R+, R+]


           f = rmutil::dboxcox(y, m = exp(param[1]), s = exp(param[2]), f = exp(param[3]))+1e-10

         } ,


         "burr" = {     ## [R+, R+, R+]


           f = actuar::dburr(y, shape1 = exp(param[1]), shape2 = exp(param[2]), scale = exp(param[3]))+1e-10

         } ,


         "burr2param" = {     ## [R+, R+]


           f = VaRES::dburr(y, a = exp(param[1]), b = exp(param[2]))+1e-10

         } ,


         "cauchy" = {     ## [R, R+]


           f = stats::dcauchy(y, location = param[1], scale = exp(param[2]))+1e-10

         } ,


         "chen" = {     ## [R+, R+]


           f = VaRES::dchen(y, b = exp(param[1]), lambda = exp(param[2]))+1e-10

         } ,


         "chi" = {     ## [R+]


           f = EnvStats::dchi(y, df = exp(param[1]))+1e-10

         } ,


         "chisquared" = {     ## [R+]


           f = stats::dchisq(y, df = exp(param[1]))+1e-10

         } ,


         "clg" = {     ## [R+, R+, R]


           f = VaRES::dclg(y, a = exp(param[1]), b = exp(param[2]), param[3])+1e-10

         } ,


         "complementarybeta" = {     ## [R+, R+]


           f = VaRES::dcompbeta(y, a = exp(param[1]), b = exp(param[2]))+1e-10

         } ,



         "dagum" = {     ## [R+, R+, R+]


           f = VGAM::ddagum(y, scale = exp(param[1]), shape1.a = exp(param[2]),
                            shape2.p = exp(param[3]))+1e-10

         } ,


         "diffzeta" = {     ## [R+, >1]


           f = VGAM::ddiffzeta(y, shape = exp(param[1]), start = 1+exp(-param[2]) )+1e-10

         } ,


         "discretegamma" = {     ## [R+, R+]


           f = extraDistr::ddgamma(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "discretelaplace" = {     ## [R, 01]


           f = extraDistr::ddlaplace(y, location = param[1], scale = 1/(1+exp(-param[2])) )+1e-10

         } ,


         "discretenormal" = {     ## [R, R+]


           f = extraDistr::ddnorm(y, mean = param[1], sd = exp(param[2]))+1e-10

         } ,


         "discreteweibull" = {     ## [01, R+]


           f = extraDistr::ddweibull(y, shape1 = 1/(1+exp(-param[1])), shape2 = exp(param[2]))+1e-10

         } ,



         "doubleweibull" = {     ## [R+, R, R+]

             f = VaRES::ddweibull(y, c = exp(param[1]), mu = param[2], sigma = exp(param[3]))+1e-10

         } ,


         "ev" = {

           ## [R, R+]
           f = VGAM::dgev(y, location = param[1], scale = exp(param[2]), shape = 0)+1e-10

         } ,


         "exponential" = {     ## [R+]


           f = stats::dexp(y, rate = exp(param[1]))+1e-10

         } ,


         "exponentialextension" = {     ## [R+, R+]


           f = VaRES::dexpext(y, lambda = exp(param[1]), a = exp(param[2]))+1e-10

         } ,



         "exponentialgeometric" = {     ## [R+, 01]


           f = VGAM::dexpgeom(y, scale = exp(param[1]), shape = 1/(1+exp(-param[2])) )+1e-10

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


           f = VGAM::dexplog(y, scale = exp(param[1]), shape = 1/(1+exp(-param[2])) )+1e-10

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


           f = VaRES::dexppois(y, b = exp(param[1]), lambda = exp(param[2]))+1e-10

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


           f = VaRES::dexppower(y, mu = param[1], sigma = exp(param[2]), a = exp(param[3]))+1e-10

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


           f = VaRES::dexpexp(y, lambda = exp(param[1]), a = exp(param[2]))+1e-10

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


           f = VaRES::dexplogis(y, a = exp(param[1]), b = exp(param[2]))+1e-10

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


           f = VaRES::dexpweibull(y, a = exp(param[1]), alpha = exp(param[2]),
                                  sigma = exp(param[3]))+1e-10

         } ,


         "F" = {     ## [R+, R+]


           f = stats::df(y, df1 = exp(param[1]), df2 = exp(param[2]))+1e-10

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]


           f = actuar::dfpareto(y, min = param[1], shape1 = exp(param[2]),
                                shape2 = exp(param[3]), shape3 = exp(param[4]),
                                scale = exp(param[5]))+1e-10

         } ,


         "fisk" = {     ## [R+, R+]

           f = VGAM::dfisk(y, scale = exp(param[1]), shape1.a = exp(param[2]))+1e-10

         } ,


         "foldednormal" = {     ## [R, R+]


           f = VGAM::dfoldnorm(y, mean = param[1], sd = exp(param[2]))+1e-10

         } ,


         "frechet" = {     ## [R+, R, R+]


           f = VGAM::dfrechet(y, shape = exp(param[1]), location = param[2],
                              scale = exp(param[3]))+1e-10

         } ,


         "gamma" = {     ## [R+, R+]


           f = stats::dgamma(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,



         "gaussian" = {     ## [R, R+]


           f = stats::dnorm(y, mean = param[1], sd = exp(param[2]))+1e-10

         } ,


         "gev" = {     ## [R, R+, R]


           f = VGAM::dgev(y, location = param[1], scale = exp(param[2]), shape = param[3])+1e-10

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


           f = VaRES::dgeninvbeta(y, a = exp(param[1]), c = exp(param[2]),
                                  d = exp(param[3]))+1e-10

         } ,


         "genlogis" = {     ## [R+, R, R+]


           f = VaRES::dgenlogis(y, a = exp(param[1]), mu = param[2], sigma = exp(param[3]))+1e-10

         } ,


         "genlogis3" = {     ## [R+, R, R+]


           f = VaRES::dgenlogis3(y, a = exp(param[1]), mu = param[2], sigma = exp(param[3]))+1e-10

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


           f = VaRES::dgenlogis4(y, a = exp(param[1]), alpha = exp(param[2]),
                                 mu = param[3], sigma = exp(param[4]))+1e-10

         } ,


         "genpowerweibull" = {     ## [R+, R+]


           f = VaRES::dgenpowerweibull(y, a = exp(param[1]), theta = exp(param[2]))+1e-10

         } ,



         "geometric" = {     ## [01]

           f = stats::dgeom(y, prob = 1/(1+exp(-param[1])) )+1e-10

         } ,




         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)


            f = GeneralizedHyperbolic::dghyp(y, mu = param[1], delta = exp(param[2]),
                                             alpha = exp(param[3]),
                                             beta = exp(param[3]-10e-15) * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 ),
                                             lambda = param[5])+1e-10
         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           f = GLDEX::dgl(y, lambda1 = param[1], lambda2 = exp(param[2]),
                          lambda3 = param[3], lambda4 = param[4])+1e-10

         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

           f = gamlss.dist::dGT(y, mu = param[1], sigma = exp(param[2]),
                                nu = exp(param[3]), tau = exp(param[4]))+1e-10

         } ,



         "gompertz" = {     ## [R+, R+]


           f = ssdtools::dgompertz(y, lscale = exp(param[1]), lshape = exp(param[2]))+1e-10

         } ,


         "gpd" = {     ## [R, R+, R]


           f = VGAM::dgpd(y, location = param[1], scale = exp(param[2]), shape = param[3])+1e-10

         } ,


         "gumbel" = {     ## [R, R+]


           f = VGAM::dgumbel(y, location = param[1], scale = exp(param[2]))+1e-10

         } ,


         "gumbel2" = {     ## [R+, R+]


           f = VGAM::dgumbelII(y, scale = exp(param[1]), shape = exp(param[2]))+1e-10

         } ,


         "halfcauchy" = {     ## [R+]


           f = extraDistr::dhcauchy(y, sigma = exp(param[1]))+1e-10

         } ,


         "halflogistic" = {     ## [R+]


           f = VaRES::dhalflogis(y, lambda = exp(param[1]))+1e-10

         } ,


         "halfnormal" = {     ## [R+]


           f = extraDistr::dhnorm(y, sigma = exp(param[1]))+1e-10

         } ,


         "halft" = {     ## [R+, R+]


           f = extraDistr::dht(y, nu = exp(param[1]), sigma = exp(param[2]))+1e-10

         } ,


         "hjorth" = {     ## [R+, R+, R+]


           f = rmutil::dhjorth(y, m = exp(param[1]), s = exp(param[2]), f = exp(param[3]))+1e-10

         } ,


         "hblaplace" = {     ## [01, R, R+]


           f = VaRES::dHBlaplace(y, a = 1/(1+exp(-param[1])), theta = param[2], phi = exp(param[3]))+1e-10

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]


             f = GeneralizedHyperbolic::dhyperb(y, mu = param[1], delta = exp(param[2]),
                                                alpha = exp(param[3]),
                                                beta = exp(param[3]) * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 ))+1e-10
         } ,


         "huber" = {     ## [R, R+]


           f = extraDistr::dhuber(y, mu = param[1], sigma = exp(param[2]))+1e-10

         } ,


         "hzeta" = {     ## [R+]

           f = VGAM::dhzeta(y, shape = exp(param[1]))+1e-10

         } ,


         "inversebeta" = {     ## [R+, R+]


           f = VaRES::dinvbeta(y, a = exp(param[1]), b = exp(param[2]))+1e-10

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


           f = actuar::dinvburr(y, shape1 = exp(param[1]), shape2 = exp(param[2]),
                                scale = exp(param[3]))+1e-10

         } ,


         "inversechisquared" = {     ## [R+]


           f = extraDistr::dinvchisq(y, nu = exp(param[1]))+1e-10

         } ,


         "inverseexponential" = {     ## [R+]


           f = actuar::dinvexp(y, scale = exp(param[1]))+1e-10

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


           f = VaRES::dinvexpexp(y, lambda = exp(param[1]), a = exp(param[2]))+1e-10

         } ,


         "inversegamma" = {     ## [R+, R+]


           f = extraDistr::dinvgamma(y, alpha = exp(param[1]), beta = exp(param[2]))+1e-10

         } ,


         "inverselomax" = {     ## [R+, R+]


           f = VGAM::dinv.lomax(y, scale = exp(param[1]), shape2.p = exp(param[2]))+1e-10

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


           f = actuar::dinvparalogis(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "inversepareto" = {     ## [R+, R+]


           f = actuar::dinvpareto(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


           f = actuar::dinvtrgamma(y, shape1 = exp(param[1]), shape2 = exp(param[2]),
                                   scale = exp(param[3]))+1e-10

         } ,


         "inverseweibull" = {     ## [R+, R+]


           f = actuar::dinvweibull(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "kumaraswamy" = {     ## [R+, R+]


           f = VGAM::dkumar(y, shape1 = exp(param[1]), shape2 = exp(param[2]))+1e-10

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


           f = VaRES::dkumexp(y, lambda = exp(param[1]), a = exp(param[2]),
                              b = exp(param[3]))+1e-10

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


           f = VaRES::dkumgamma(y, a = exp(param[1]), b = exp(param[2]),
                                c = exp(param[3]), d = exp(param[4]))+1e-10

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


           f = VaRES::dkumgumbel(y, a = exp(param[1]), b = exp(param[2]), mu = param[3],
                                 sigma = exp(param[4]))+1e-10

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


           f = VaRES::dkumhalfnorm(y, sigma = exp(param[1]), a = exp(param[2]),
                                   b = exp(param[3]))+1e-10

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]


           f = VaRES::dkumloglogis(y, a = exp(param[1]), b = exp(param[2]), alpha = exp(param[3]),
                                   beta = exp(param[4]))+1e-10

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]


           f = VaRES::dkumnormal(y, mu = param[1], sigma = exp(param[2]), a = exp(param[3]),
                                 b = exp(param[4]))+1e-10

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


           f = VaRES::dkumweibull(y, a = exp(param[1]), b = exp(param[2]), alpha = exp(param[3]),
                                  sigma = exp(param[4]))+1e-10

         } ,



         "laplace" = {     ## [R, R+]


           f = extraDistr::dlaplace(y, mu = param[1], sigma = exp(param[2]))+1e-10

         } ,


         "levy" = {     ## [R, R+]


           f = rmutil::dlevy(y, m = param[1], s = exp(param[2]))+1e-10

         } ,


         "linearfailurerate" = {     ## [R+, R+]


           f = VaRES::dlfr(y, a = exp(param[1]), b = exp(param[2]))+1e-10

         } ,


         "lindley" = {     ## [R+]


           f = VGAM::dlind(y, theta = exp(param[1]))+1e-10

         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


           f = VaRES::dLNbeta(y, lambda = exp(param[1]), a = exp(param[2]),
                              b = exp(param[3]))+1e-10

         } ,


         "logcauchy" = {     ## [R, R+]


           f = VaRES::dlogcauchy(y, mu = param[1], sigma = exp(param[2]))+1e-10

         } ,



         "loggamma" = {     ## [R, R+, R+]


           f = VGAM::dlgamma(y, location = param[1], scale = exp(param[2]),
                             shape = exp(param[3]))+1e-10

         } ,



         "loggumbel" = {     ## [R, R+]


           f = ssdtools::dlgumbel(y, llocation = param[1], lscale = exp(param[2]))+1e-10

         } ,



         "loglog" = {     ## [R+, >1]


           f = VaRES::dloglog(y, a = exp(param[1]), lambda = 1 + exp(-param[2]) )+1e-10

         } ,


         "loglogistic" = {     ## [R+, R+]


           f = actuar::dllogis(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "lognormal" = {     ## [R, R+]


           f = stats::dlnorm(y, meanlog = param[1], sdlog = exp(param[2]))+1e-10

         } ,


         "lognormal3" = {     ## [R, R+, R]


           f = EnvStats::dlnorm3(y, meanlog = param[1], sdlog = exp(param[2]),
                                 threshold = param[3])+1e-10

         } ,


         "logistic" = {     ## [R, R+]


           f = stats::dlogis(y, location = param[1], scale = exp(param[2]))+1e-10

         } ,


         "logisticexponential" = {     ## [R+, R+]


           f = VaRES::dlogisexp(y, lambda = exp(param[1]), a = exp(param[2]))+1e-10

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


           f = VaRES::dlogisrayleigh(y, a = exp(param[1]), lambda = exp(param[2]))+1e-10

         } ,


         "logseries" = {     ## [01]

           f = extraDistr::dlgser(y, theta = 1/(1+exp(-param[1])) )+1e-10

         } ,


         "lomax" = {     ## [R+, R+]


           f = VGAM::dlomax(y, scale = exp(param[1]), shape3.q = exp(param[2]))+1e-10

         } ,


         "makeham" = {     ## [R+, R+, R+]


           f = VGAM::dmakeham(y, scale = exp(param[1]), shape = exp(param[2]),
                              epsilon = exp(param[3]))+1e-10

         } ,


         "maxwell" = {     ## [R+]


           f = VGAM::dmaxwell(y, rate = exp(param[1]))+1e-10

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


           f = VaRES::dMlaplace(y, theta = param[1], phi = exp(param[2]),
                                psi = exp(param[3]))+1e-10

         } ,



         "moexponential" = {     ## [R+, R+]


           f = VaRES::dmoexp(y, lambda = exp(param[1]), a = exp(param[2]))+1e-10

         } ,


         "moweibull" = {     ## [R+, R+, R+]


           f = VaRES::dmoweibull(y, a = exp(param[1]), b = exp(param[2]), lambda = exp(param[3]))+1e-10

         } ,


         "nakagami" = {     ## [R+, R+]


           f = VGAM::dnaka(y, scale = exp(param[1]), shape = exp(param[2]))+1e-10

         } ,


         "ncchisquared" = {     ## [R+, R+]


           f = stats::dchisq(y, df = exp(param[1]), ncp = exp(param[2]))+1e-10

         } ,


         "ncF" = {     ## [R+, R+, R+]


           f = stats::df(y, df1 = exp(param[1]), df2 = exp(param[2]), ncp = exp(param[3]))+1e-10

         } ,


         "negativebinomial" = {     ## [N+, 01]

           f = stats::dnbinom(y, size = size, prob = 1/(1+exp(-param[1])) )+1e-10

         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]

           f = GeneralizedHyperbolic::dnig(y, mu = param[1], delta = exp(param[2]),
                                                alpha = exp(param[3]),
                                                beta = exp(param[3]) * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 ))+1e-10

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


           f = extraDistr::dnsbeta(y, shape1 = exp(param[1]), shape2 = exp(param[2]),
                                   min = param[3], max = param[4])+1e-10

         } ,


         "paralogistic" = {     ## [R+, R+]


           f = VGAM::dparalogistic(y, scale = exp(param[1]), shape1.a = exp(param[2]))+1e-10

         } ,


         "pareto" = {     ## [R+]


           f = extraDistr::dpareto(y, a = exp(param[1]), b = min(y))+1e-10

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


           f = VaRES::dparetostable(y, lambda = exp(param[1]), nu = exp(param[2]),
                                    sigma = exp(param[3]))+1e-10

         } ,


         "pareto1" = {     ## [R+, R+]


           f = VGAM::dparetoI(y, scale = exp(param[1]), shape = exp(param[2]))+1e-10

         } ,


         "pareto2" = {     ## [R, R+, R+]


           f = VGAM::dparetoII(y, location = param[1], scale = exp(param[2]),
                               shape = exp(param[3]))+1e-10

         } ,


         "pareto3" = {     ## [R, R+, R+]


           f = VGAM::dparetoIII(y, location = param[1], scale = exp(param[2]),
                                inequality = exp(param[3]))+1e-10

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


           f = VGAM::dparetoIV(y, location = param[1], scale = exp(param[2]), inequality = exp(param[3]),
                               shape = exp(param[4]))+1e-10

         } ,


         "perks" = {     ## [R+, R+]


           f = VGAM::dperks(y, scale = exp(param[1]), shape = exp(param[2]))+1e-10

         } ,


         "pctalaplace" = {     ## [R+, R]


           f = VaRES::dPCTAlaplace(y, a = exp(param[1]), theta = param[2])+1e-10

         } ,


         "poisson" = {     ## [R+]


           f = stats::dpois(y, lambda = exp(param[1]))+1e-10

         } ,


         "power1" = {     ## [R+]


           f = VaRES::dpower1(y, a = exp(param[1]))+1e-10

         } ,


         "power2" = {     ## [R+]


           f = VaRES::dpower2(y, b = exp(param[1]))+1e-10

         } ,


         "powerdistribution" = {     ## [R+, R+]


           f = extraDistr::dpower(y, alpha = exp(param[1]), beta = exp(param[1]))+1e-10

         } ,


         "powerexponential" = {     ## [R, R+, R+]


           f = rmutil::dpowexp(y, m = param[1], s = exp(param[2]), f = exp(param[3]))+1e-10

         } ,


         "rayleigh" = {     ## [R+]


           f = VGAM::drayleigh(y, scale = exp(param[1]))+1e-10

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


           f = VaRES::drgamma(y, a = exp(param[1]), theta = param[2],
                              phi = exp(param[3]))+1e-10

         } ,



         "rice" = {     ## [R+, R+]


           f = VGAM::drice(y, sigma = exp(param[1]), vee = exp(param[2]))+1e-10

         } ,


         "scaledchisquared" = {     ## [R+, R+]


           f = extraDistr::dinvchisq(y, nu = exp(param[1]), tau = exp(param[2]))+1e-10

         } ,


         "schabe" = {     ## [R+, R+]


           f = VaRES::dschabe(y, gamma = exp(param[1]), theta = exp(param[2]))+1e-10

         } ,



         "simplex" = {     ## [01, R+]


           f = rmutil::dsimplex(y, m = 1/(1+exp(-param[1])), s = exp(param[2]))+1e-10

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


           f = rmutil::dskewlaplace(y, m = param[1], s = exp(param[2]), f = exp(param[3]))+1e-10

         } ,


         "skewedt" = {     ## [R+, R+]


           f = skewt::dskt(y, df = exp(param[1]), gamma = exp(param[2]))+1e-10

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+ (<25)]

           f = sn::dst(y, xi = param[1], omega = exp(param[2]), alpha = param[3],
                       nu = 25/(1+exp(-param[4])) )+1e-10

         } ,


         "skewednormal" = {     ## [R, R+, R]

           f = sn::dsn(y, xi = param[1], omega = exp(param[2]), alpha = param[3])+1e-10

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           f = sgt::dsgt(y, mu = param[1], sigma = exp(param[2]),
                         lambda = (exp(2*param[3])-1)/(exp(2*param[3])+1),
                         p = 1+exp(-param[4]))+1e-10

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           f = gamlss.dist::dSEP(y, mu = param[1], sigma = exp(param[2]),
                                 nu = param[3], tau = exp(param[4]))+1e-10

         } ,



         "slash" = {     ## [R, R+]


           f = extraDistr::dslash(y, mu = param[1], sigma = exp(param[2]))+1e-10

         } ,


         "stacy" = {     ## [R+, R+, R+]


           f = VaRES::dstacygamma(y, gamma = exp(param[1]), c = exp(param[2]),
                                  theta = exp(param[3]))+1e-10

         } ,


         "t" = {     ## [R, R+, R+(<25)]


           f =  stats::dt((y-param[1])/exp(param[2]), df = 25/(1+exp(-param[3])) ) / exp(param[2])+1e-10

         } ,



         "tobit" = {     ## [R, R+]


           f = VGAM::dtobit(y, mean = param[1], sd = exp(param[2]))+1e-10

         } ,


         "topple" = {     ## [01]


           f = VGAM::dtopple(y, shape = 1/(1+exp(-param[1])))+1e-10

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]


           f = actuar::dtrbeta(y, shape1 = exp(param[1]), shape2 = exp(param[2]), shape3 = exp(param[3]),
                               scale = exp(param[4]))+1e-10

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]


           f = actuar::dtrgamma(y, shape1 = exp(param[1]), shape2 = exp(param[2]),
                                scale = exp(param[3]))+1e-10

         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


           f = extraDistr::dtnorm(y, mean = param[1], sd = exp(param[2]), a = param[3],
                                  b = param[4])+1e-10

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


           f = VGAM::dtruncpareto(y, lower = exp(param[2]), upper = exp(param[2]),
                                  shape = exp(param[3]))+1e-10

         } ,


         "twosidedpower" = {     ## [01, R+]


           f = rmutil::dtwosidedpower(y, m = 1/(1+exp(-param[1])), s =  exp(param[2]))+1e-10

         } ,


         "wald" = {     ## [R, R+]


           f = extraDistr::dwald(y, mu = param[1], lambda = exp(param[2]))+1e-10

         } ,


         "weibull" = {     ## [R+, R+]


           f = stats::dweibull(y, shape = exp(param[1]), scale = exp(param[2]))+1e-10

         } ,


         "xie" = {     ## [R+, R+, R+]


           f = VaRES::dxie(y, a = exp(param[1]), b = exp(param[2]), lambda = exp(param[3]))+1e-10

         } ,


         "yules" = {     ## [R+] mais > 0.5


           f = VGAM::dyules(y, shape = 0.5+exp(-param[1]))+1e-10

         } ,


  )

  f[is.nan(f)] <- 1e-10
  return(f)


}





