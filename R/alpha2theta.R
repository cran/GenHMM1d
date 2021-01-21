#'@title Transforms unconstrained parameters to constrainted  parameters
#'
#'@description This function computes the constrained  parameters theta of a univariate distribution from the unconstrainted parameter alpha.
#'
#'@param family distribution name; run the command distributions() for help
#'@param param  unconstrained parameters of the univariate distribution
#'
#'@return \item{theta}{constrained parameters}
#'
#'
#'@export
#'@keywords internal


alpha2theta<-function(family,param){


  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])
           theta[3] = 1/(1+exp(-param[3]))

         } ,


         "asymlaplace" = {    ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "asympower" = {    ## [01, R+, R+]

           theta = matrix(0,1,3)
           theta[1] = 1/(1+exp(-param[1]))
           theta[2:3] = exp(param[2:3])

         } ,


         "asymt" = {    ## [R+, R+, 01, R]

           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = 1/(1+exp(-param[3]))
           theta[4] = param[4]

         } ,


         "beard" = {    ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "benini" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "benford" = {     ## [1 ou 2]

           theta = matrix(0,1,1)
           if (param[1]>2){
             theta[1] = 2
           } else if (param[1]<1){
             theta[1] = 1
           } else {
             theta[1] = round(param[1])
           }

         } ,


         "bernoulli" = {     ## [01]

           theta = matrix(0,1,1)
           theta[1] = 1/(1+exp(-param[1]))

         } ,


         "beta" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "betabinomial" = {     ## [N+, R+, R+]

           theta = matrix(0,1,2)
           #theta[1] = round(param[1])
           theta[1:2] = exp(param[1:2])

         } ,


         "betageometric" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]

           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]


           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]
           theta[4] = exp(param[4])

         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]

           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]
           theta[4] = exp(param[4])

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]
           theta[4] = exp(param[4])

         } ,


         "betaprime" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "bhattacharjee" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "binomial" = {     ## [N+, 01]

           theta = matrix(0,1,1)
           theta[1] = 1 / (1+exp(-param[1]))

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]

           theta = matrix(0,1,3)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]

         } ,


         "boxcox" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "burr" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "burr2param" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "cauchy" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "chen" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "chi" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "chisquared" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "clg" = {     ## [R+, R+, R]


           theta = matrix(0,1,3)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]

         } ,


         "complementarybeta" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,



         "dagum" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "diffzeta" = {     ## [R+, >1]

           theta = matrix(0,1,2)
           theta[1] = exp(param[1])
           theta[2] = 1+exp(-param[2])

         } ,


         "discretegamma" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "discretelaplace" = {     ## [R, 01]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = 1/(1+exp(-param[2]))

         } ,


         "discretenormal" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "discreteweibull" = {     ## [01, R+]

           theta = matrix(0,1,2)
           theta[1] = 1/(1+exp(-param[1]))
           theta[2] = exp(param[2])

         } ,


         "doubleweibull" = {     ## [R+, R, R+]

           theta = matrix(0,1,3)
           theta[1] = exp(param[1])
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,


         "ev" = {

           ## [R, R+]
           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "exponential" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "exponentialextension" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,



         "exponentialgeometric" = {     ## [R+, 01]

           theta = matrix(0,1,2)
           theta[1] = exp(param[1])
           theta[2] = 1/(1+exp(-param[2]))

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


           theta = matrix(0,1,2)
           theta[1] = exp(param[1])
           theta[2] = 1/(1+exp(-param[2]))

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "F" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]

           theta = matrix(0,1,5)
           theta[1] = param[1]
           theta[2:5] = exp(param[2:5])

         } ,


         "fisk" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "foldednormal" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "frechet" = {     ## [R+, R, R+]

           theta = matrix(0,1,3)
           theta[1] = exp(param[1])
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,


         "gamma" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "gammapoisson" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "gaussian" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "gev" = {     ## [R, R+, R]

           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "genlogis" = {     ## [R+, R, R+]


           theta = matrix(0,1,3)
           theta[1] = exp(param[1])
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,


         "genlogis3" = {     ## [R+, R, R+]


           theta = matrix(0,1,3)
           theta[1] = exp(param[1])
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]
           theta[4] = exp(param[4])

         } ,


         "genpowerweibull" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "geometric" = {     ## [01]

           theta = matrix(0,1,1)
           theta[1] = 1/(1+exp(-param[1]))

         } ,




         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

           theta = matrix(0,1,5)
           theta[1] = param[1]
           theta[2] = max( c( 10e-12, exp(param[2]) ) )
           theta[3] = max( c( 10e-12, exp(param[3]) ) )
           theta[4] = theta[3] * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 )
           theta[5] = param[5]

         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]
           theta[4] = param[4]

         } ,



         "generalizedt" = {     ## [R, R+, R+, R+]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2:4] = exp(param[2:4])

         } ,



         "gompertz" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "gpd" = {     ## [R, R+, R]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]

         } ,


         "gumbel" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "gumbel2" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "halfcauchy" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "halflogistic" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "halfnormal" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "halft" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "hjorth" = {     ## [R+, R+, R+]

           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "hblaplace" = {     ## [01, R, R+]

           theta = matrix(0,1,3)
           theta[1] = 1/(1+exp(-param[1]))
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]

             theta = matrix(0,1,4)
             theta[1] = param[1]
             theta[2] = exp(param[2])
             theta[3] = exp(param[3])
             theta[4] = theta[3] * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 )

         } ,


         "huber" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "hzeta" = {     ## [R+]

           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "inversebeta" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "inversechisquared" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "inverseexponential" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inversegamma" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inverselomax" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inversepareto" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "inverseweibull" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "kumaraswamy" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3] = param[3]
           theta[4] = exp(param[4])

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2:4] = exp(param[2:4])

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,



         "laplace" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "levy" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "linearfailurerate" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "lindley" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "logcauchy" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,



         "loggamma" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "loggumbel" = {     ## [R, R+]

           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "loglaplace" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "loglog" = {     ## [R+, 01]


           theta = matrix(0,1,2)
           theta[1] = exp(param[1])
           theta[2] = 1+exp(-param[2])

         } ,


         "loglogistic" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "lognormal" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "lognormal3" = {     ## [R, R+, R]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]

         } ,


         "logistic" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "logisticexponential" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "logseries" = {     ## [01]

           theta = matrix(0,1,1)
           theta[1] = 1/(1+exp(-param[1]))

         } ,


         "lomax" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "makeham" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "maxwell" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "moexponential" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "moweibull" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "nakagami" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "ncchisquared" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "ncF" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "negativebinomial" = {     ## [N+, 01]

           theta = matrix(0,1,1)
           theta[1] = 1/(1+exp(-param[1]))

         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]

             theta = matrix(0,1,4)
             theta[1] = param[1]
             theta[2] = exp(param[2])
             theta[3] = exp(param[3])
             theta[4] = theta[3] * ( exp(2*param[4])-1 ) / ( exp(2*param[4])+1 )

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


           theta = matrix(0,1,4)
           theta[1:2] = exp(param[1:2])
           theta[3:4] = param[3:4]

         } ,



         "paralogistic" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "pareto" = {     ## [R+, R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "pareto1" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "pareto2" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "pareto3" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2:4] = exp(param[2:4])

         } ,


         "perks" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "pctalaplace" = {     ## [R+, R]

           theta = matrix(0,1,2)
           theta[1] = exp(param[1])
           theta[2] = param[2]

         } ,


         "poisson" = {     ## [R+]

           theta = matrix(0,1,1)
           theta[1] = exp(param[1])
         } ,



         "power1" = {     ## [R+]

           theta = matrix(0,1,1)
           theta[1] = exp(param[1])
         } ,


         "power2" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "powerdistribution" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "powerexponential" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "rayleigh" = {     ## [R+]


           theta = matrix(0,1,1)
           theta[1] = exp(param[1])

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


           theta = matrix(0,1,3)
           theta[1] = exp(param[1])
           theta[2] = param[2]
           theta[3] = exp(param[3])

         } ,



         "rice" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "scaledchisquared" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "schabe" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,



         "simplex" = {     ## [01, R+]


           theta = matrix(0,1,2)
           theta[1] = 1/(1+exp(-param[1]))
           theta[2] = exp(param[2])

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2:3] = exp(param[2:3])

         } ,


         "skewedt" = {     ## [R+, R+]

           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]
           theta[4] = 25/(1+exp(-param[4]))

         } ,


         "skewednormal" = {     ## [R, R+, R, R]

           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = (exp(2*param[3])-1)/(exp(2*param[3])+1)
           theta[4] = 1+exp(-param[4])

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = param[3]
           theta[4] = max( c( 10e-12, exp(param[4]) ) )

         } ,


         "slash" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "stacy" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "t" = {     ## [R, R+, R+ (<25)]


           theta = matrix(0,1,3)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3] = 25/(1+exp(-param[3]))

         } ,



         "tobit" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "topple" = {     ## [01]

           theta = matrix(0,1,1)
           theta[1] = 1/(1+exp(-param[1]))

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]

           theta = matrix(0,1,4)
           theta[1:4] = exp(param[1:4])

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]

           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


           theta = matrix(0,1,4)
           theta[1] = param[1]
           theta[2] = exp(param[2])
           theta[3:4] = param[3:4]

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "twosidedpower" = {     ## [01, R+]


           theta = matrix(0,1,2)
           theta[1] = 1/(1+exp(-param[1]))
           theta[2] = exp(param[2])

         } ,


         "wald" = {     ## [R, R+]


           theta = matrix(0,1,2)
           theta[1] = param[1]
           theta[2] = exp(param[2])

         } ,


         "weibull" = {     ## [R+, R+]


           theta = matrix(0,1,2)
           theta[1:2] = exp(param[1:2])

         } ,


         "xie" = {     ## [R+, R+, R+]


           theta = matrix(0,1,3)
           theta[1:3] = exp(param[1:3])

         } ,


         "yules" = {     ## [R+] >0.5

           theta = matrix(0,1,1)
           theta[1] = 0.5+exp(-param[1])

         } ,



  )


  return(theta)


}





