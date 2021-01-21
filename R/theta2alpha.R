#'@title Transform constrained parameters  to unconstrained parameters
#'
#'@description This function computes the unconstrained parameters alpha of a univariate distribution corresponding to constrainted parameters theta.
#'
#'@param family distribution name; run the command distributions() for help
#'@param param  constrained parameters of the univariate distribution
#'
#'@return \item{alpha}{constrained parameters}
#'
#'
#'
#'@export
#'@keywords internal


theta2alpha<-function(family,param){


  switch(family,

         "asymexppower" = {    ## [R+, R+, 01]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])
           alpha[3] = log(param[3]/(1-param[3]))

         } ,


         "asymlaplace" = {    ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "asympower" = {    ## [01, R+, R+]

           alpha = matrix(0,1,3)
           alpha[1] = log(param[1]/(1-param[1]))
           alpha[2:3] = log(param[2:3])

         } ,


         "asymt" = {    ## [R+, R+, 01]

           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = log(param[3]/(1-param[3]))
           alpha[4] = param[4]

         } ,


         "beard" = {    ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "benini" = {     ## [R, R+]

           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "benford" = {     ## [1 ou 2]

           alpha = matrix(0,1,1)
           if (param[1]>log(2)){
             alpha[1] = log(2)
           } else if (param[1]<0){
             alpha[1] = 0
           } else {
             alpha[1] = log(param[1])
           }

         } ,


         "bernoulli" = {     ## [01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,


         "beta" = {     ## [R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,



         "betabinomial" = {     ## [N+, R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "betageometric" = {     ## [R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "betanegativebinomial" = {     ## [N+, R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "betaburr" = {     ## [R+, R+, R+, R+]

           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "betaburr7" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "betaexponential" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "betafrechet" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "betagompertz" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "betagumbel" = {     ## [R+, R+, R, R+]


           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,


         "betagumbel2" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "betalognormal" = {     ## [R+, R+, R, R+]

           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,


         "betalomax" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,



         "betanormal" = {     ## [R+, R+, R, R+]


           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,


         "betaprime" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "betaweibull" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "bhattacharjee" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "binomial" = {     ## [N+, 01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,


         "birnbaumsaunders" = {     ## [R+, R+, R]

           alpha = matrix(0,1,3)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]

         } ,


         "boxcox" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "burr" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "burr2param" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "cauchy" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "chen" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "chi" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "chisquared" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "clg" = {     ## [R+, R+, R]


           alpha = matrix(0,1,3)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]

         } ,


         "complementarybeta" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,



         "dagum" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "diffzeta" = {     ## [R+, >1]

           alpha = matrix(0,1,2)
           alpha[1] = log(param[1])
           alpha[2] = log(exp(-param[2]) + 1)

         } ,


         "discretegamma" = {     ## [R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "discretelaplace" = {     ## [R, 01]

           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])/(1-param[2])

         } ,


         "discretenormal" = {     ## [R, R+]

           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "discreteweibull" = {     ## [01, R+]

           alpha = matrix(0,1,2)
           alpha[1] = log(param[1]/(1-param[1]))
           alpha[2] = log(param[2])

         } ,


         "doubleweibull" = {     ## [R+, R, R+]

           alpha = matrix(0,1,3)
           alpha[1] = log(param[1])
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,


         "ev" = {

           ## [R, R+]
           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "exponential" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "exponentialextension" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,



         "exponentialgeometric" = {     ## [R+, 01]

           alpha = matrix(0,1,2)
           alpha[1] = log(param[1])
           alpha[2] = log(param[2]/(1-param[2]))

         } ,


         "exponentiallogarithmic" = {     ## [R+, 01]


           alpha = matrix(0,1,2)
           alpha[1] = log(param[1])
           alpha[2] = log(param[2]/(1-param[2]))

         } ,


         "exponentialpoisson" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "exponentialpower" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "exponentiatedexponential" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "exponentiatedlogistic" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "exponentiatedweibull" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,



         "F" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]

           alpha = matrix(0,1,5)
           alpha[1] = param[1]
           alpha[2:5] = log(param[2:5])

         } ,


         "fisk" = {     ## [R+, R+]

            alpha = matrix(0,1,2)
            alpha[1:2] = log(param[1:2])
         } ,


         "foldednormal" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "frechet" = {     ## [R+, R, R+]

           alpha = matrix(0,1,3)
           alpha[1] = log(param[1])
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,


         "gamma" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "gammapoisson" = {     ## [R+, R+]

           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "gaussian" = {     ## [R, R+]

           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "gev" = {     ## [R, R+, R]

           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]

         } ,


         "geninvbeta" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "genlogis" = {     ## [R+, R, R+]


           alpha = matrix(0,1,3)
           alpha[1] = log(param[1])
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,


         "genlogis3" = {     ## [R+, R, R+]


           alpha = matrix(0,1,3)
           alpha[1] = log(param[1])
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,


         "genlogis4" = {     ## [R+, R+, R, R+]


           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,


         "geometric" = {     ## [01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,



         "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]  [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

           alpha = matrix(0,1,5)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = log(param[3])
           alpha[4] = log( (param[4]+param[3]) / (param[3]-param[4]) )
           alpha[5] = param[5]

         } ,


         "generalizedlambda" = {     ## [R, R+, R, R]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]
           alpha[4] = param[4]

         } ,


         "generalizedt" = {     ## [R, R+, R+, R+]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2:4] = log(param[2:4])

         } ,


         "gompertz" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "gpd" = {     ## [R, R+, R]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]

         } ,


         "gumbel" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "gumbel2" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "halfcauchy" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "halflogistic" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "halfnormal" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "halft" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "hjorth" = {     ## [R+, R+, R+]

           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "hblaplace" = {     ## [01, R, R+]

           alpha = matrix(0,1,3)
           alpha[1] = log(param[1]/(1-param[1]))
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,


         "hyperbolic" = {     ## [R, R+, R+, R]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = log(param[3])
           alpha[4] = log( (param[4]+param[3]) / (param[3]-param[4]) )

         } ,


         "hzeta" = {     ## [R+]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "huber" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "inversebeta" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inverseburr" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "inversechisquared" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "inverseexponential" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "inverseexpexponential" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inversegamma" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inverselomax" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inverseparalogistic" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inversepareto" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "inversetransformedgamma" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "inverseweibull" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "kumaraswamy" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "kumaraswamyexponential" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "kumaraswamygamma" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "kumaraswamygumbel" = {     ## [R+, R+, R, R+]


           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,


         "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "kumaraswamyloglogistic" = {     ## [R+, R+, R, R+]

           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "kumaraswamynormal" = {     ## [R, R+, R+, R+]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2:4] = log(param[2:4])

         } ,


         "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,



         "laplace" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "levy" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "linearfailurerate" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "lindley" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "libbynovickbeta" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "logcauchy" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,



         "loggamma" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "loggumbel" = {     ## [R, R+]

           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "loglaplace" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "loglog" = {     ## [R+, >1]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(exp(-param[2]) + 1)

         } ,


         "loglogistic" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "lognormal" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "lognormal3" = {     ## [R, R+, R]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]

         } ,


         "logistic" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "logisticexponential" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "logisticrayleigh" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "logseries" = {     ## [01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,


         "lomax" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "makeham" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "maxwell" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "mcgilllaplace" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "moexponential" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "moweibull" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "nakagami" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "ncchisquared" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "ncF" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "negativebinomial" = {     ## [N+, 01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,


         "normalinversegaussian" = {     ## [R, R+, R+, R]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = log(param[3])
           alpha[4] = log( (param[4]+param[3]) / (param[3]-param[4]) )

         } ,


         "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]


           alpha = matrix(0,1,4)
           alpha[1:2] = log(param[1:2])
           alpha[3:4] = param[3:4]

         } ,



         "paralogistic" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "pareto" = {     ## [R+, R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "paretopositivestable" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "pareto1" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "pareto2" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "pareto3" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "pareto4" = {     ## [R, R+, R+, R+]


           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2:4] = log(param[2:4])

         } ,


         "perks" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "pctalaplace" = {     ## [R+, R]


           alpha = matrix(0,1,2)
           alpha[1] = log(param[1])
           alpha[2] = param[2]

         } ,


         "poisson" = {     ## [R+]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])
         } ,



         "power1" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])
         } ,


         "power2" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "powerdistribution" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "powerexponential" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "rayleigh" = {     ## [R+]


           alpha = matrix(0,1,1)
           alpha[1] = log(param[1])

         } ,


         "reflectedgamma" = {     ## [R+, R, R+]


           alpha = matrix(0,1,3)
           alpha[1] = log(param[1])
           alpha[2] = param[2]
           alpha[3] = log(param[3])

         } ,



         "rice" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "scaledchisquared" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "schabe" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,



         "simplex" = {     ## [01, R+]


           alpha = matrix(0,1,2)
           alpha[1] = log(param[1]/ (1-param[1]))
           alpha[2] = log(param[2])

         } ,


         "skewedlaplace" = {     ## [R, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2:3] = log(param[2:3])

         } ,


         "skewedt" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "skewedtfourparam" = {     ## [R, R+, R, R+(<25)]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]
           alpha[4] = log(param[4]/(25-param[4]))

         } ,


         "skewednormal" = {     ## [R, R+, R, R]

           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]

         } ,


         "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = log( (param[3]+param[2]) / (param[2]-param[3]) )
           alpha[4] = log(exp(-param[4]) + 1)

         } ,


         "skewedexponentialpower" = {     ## [R, R+, R, R+]

           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = param[3]
           alpha[4] = log(param[4])

         } ,



         "slash" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "stacy" = {     ## [R+, R+, R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "t" = {     ## [R, R+, R+ (<25)]


           alpha = matrix(0,1,3)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3] = log(param[3]/(25-param[3]))

         } ,



         "tobit" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "topple" = {     ## [01]

           alpha = matrix(0,1,1)
           alpha[1] = log(param[1]/(1-param[1]))

         } ,


         "transformedbeta" = {     ## [R+, R+, R+, R+]

           alpha = matrix(0,1,4)
           alpha[1:4] = log(param[1:4])

         } ,


         "transformedgamma" = {     ## [R+, R+, R+]

           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,



         "truncatednormal" = {     ## [R, R+, R(min), R(max)]


           alpha = matrix(0,1,4)
           alpha[1] = param[1]
           alpha[2] = log(param[2])
           alpha[3:4] = param[3:4]

         } ,


         "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]


           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "twosidedpower" = {     ## [01, R+]


           alpha = matrix(0,1,2)
           alpha[1] = log(param[1]/(1-param[1]))
           alpha[2] = log(param[2])

         } ,


         "wald" = {     ## [R, R+]


           alpha = matrix(0,1,2)
           alpha[1] = param[1]
           alpha[2] = log(param[2])

         } ,


         "weibull" = {     ## [R+, R+]


           alpha = matrix(0,1,2)
           alpha[1:2] = log(param[1:2])

         } ,


         "xie" = {     ## [R+, R+, R+]

           alpha = matrix(0,1,3)
           alpha[1:3] = log(param[1:3])

         } ,


         "yules" = {     ## [R+] >0.5

           alpha = matrix(0,1,1)
           alpha[1] = log(exp(-param[1]) + 0.5)

         } ,


  )


  return(alpha)


}





