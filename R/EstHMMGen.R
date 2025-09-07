#'#############################################################################
#'@title Estimation of univariate hidden Markov model
#'
#'@description This function estimates the parameters from a univariate hidden Markov model
#'
#'@param y       observations; (n x 1)
#'@param ZI      1 if zero-inflated, 0 otherwise (default)
#'@param reg     number of regimes (including zero-inflated; must be > ZI)
#'@param family  distribution name; run the function distributions() for help
#'@param start   starting parameters for the estimation; (1 x p)
#'@param max_iter  maximum number of iterations of the EM algorithm; suggestion 10000
#'@param eps     precision (stopping criteria); suggestion 0.001.
#'@param size    additional parameter for some discrete distributions; run the command distributions() for help
#'@param theta0  initial parameters for each regimes; (r x p), default is NULL
#'@param graph   TRUE a graph, FALSE otherwise (default); only for continuous distributions
#'
#'@return \item{theta}{ estimated parameters; (r x p)}
#'@return \item{Q}{estimated transition matrix for the regimes; (r x r)}
#'@return \item{eta}{conditional probabilities of being in regime k at time t given observations up to time t; (n x r) }
#'@return \item{lambda}{conditional probabilities of being in regime k at time t given all observations; (n x r)}
#'@return \item{U}{pseudo-observations that should be uniformly distributed under the null hypothesis}
#'@return \item{cvm}{cramer-von-Mises statistic for goodness-of-fit}
#'@return \item{W}{matrix of Rosenblatt transforms; (n x r)}
#'@return \item{LL}{log-likelihood}
#'@return \item{nu}{stationary distribution}
#'@return \item{AIC}{Akaike information criterion}
#'@return \item{BIC}{Bayesian information criterion}
#'@return \item{CAIC}{consistent Akaike information criterion}
#'@return \item{AICcorrected}{Akaike information criterion corrected}
#'@return \item{HQC}{Hannan-Quinn information criterion}
#'@return \item{stats}{empirical means and standard deviation of each regimes using lambda}
#'@return \item{pred_l}{estimated regime using lambda}
#'@return \item{pred_e}{estimated regime using eta}
#'@return \item{runs_l}{estimated number of runs using lambda}
#'@return \item{runs_e}{estimated number of runs using eta}
#'
#'@export
#'
#'@examples
#'family = "gaussian"
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2) ;
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2) ;
#'y = SimHMMGen(theta, Q=Q, family=family,  n=100)$SimData
#'est = EstHMMGen(y, ZI=0,reg=2, family=family)
#'
#'
#'
#'
#'
#'
EstHMMGen<-function(y, ZI=0, reg, family, start=0, max_iter=10000, eps=1e-4, size=0, theta0=NULL, graph=FALSE){
  ninit=100   #minimum number of iterations
  n = length(y)

  un = 1+ZI;

  reg1 = reg-ZI;
  if(ZI==1)
  {
    y0=y[y!=0]
    n00 = length(y0)
  } else {
    n00=n
     y0=y}
  n0 = floor((n00/reg1))
  ind0   = 1:n0
  if(is.null(theta0)){
    p=0
  }  else {
    p = dim(theta0)[2]
  }



  if ( is.null(theta0) || p==1){

    switch(family,

           "asymexppower" = {    ## [R+, R+, 01]
             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(1, 5, 1), p3=c(0.1, 0.9, 0.1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)
             } else {
               theta00 = start
             }

             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "asymlaplace" = {    ## [R, R+, R+]
             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.25), p2=c(1, 10, 1), p3=c(1, 1, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)
             } else {
               theta00=start
             }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "asympower" = {    ## [01, R+, R+]
             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.1, 0.9, 0.1), p2=c(1, 5, 1), p3=c(0.5, 5, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)
             }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "asymt" = {    ## [R+, R+, 01, R]
             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(3, 5, 1), p3=c(0.1, 0.9, 0.1), p4=c(-2,2,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "beard" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "benini" = {     ## [R, R+]
             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (mean(y)>0){
               p1=c(min(y)-0.05, min(y)+0.05, 0.05)
             } else{
               p1=c(max(y)-0.05, max(y)+0.05, 0.05)
             }
             if (start == 0){
               params = list(p1=p1, p2=c(0.25, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "bernoulli" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1/0.9), log(0.9/0.1)))$minimum
             }
           } ,



           "beta" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 3, 1), p2=c(1, 3, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betabinomial" = {     ## [N+, R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start==0){
               params = list(p1=c(size, size, 1), p2=c(1, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta0 = GridSearchS0(family, y, params,size)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa, size)  ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betageometric" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa)  ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betanegativebinomial" = {     ## [N+, R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(size, size, 1), p2=c(1, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta0 = GridSearchS0(family, y, params,size)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa, size)  ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betaburr" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betaburr7" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betaexponential" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betafrechet" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betagompertz" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 4, 0.5), p2=c(1, 4, 1), p3=c(1,4,1), p4=c(1,4,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betagumbel" = {     ## [R+, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(-2,2,1), p4=c(0.5,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betagumbel2" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 0.5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betalognormal" = {     ## [R+, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(-1,1,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "betalomax" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,




           "betanormal" = {     ## [R+, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0,0,1), p4=c(1,3,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betaprime" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,3,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "betaweibull" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 0.5), p2=c(1, 5, 1), p3=c(1,5,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "bhattacharjee" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 1), p2=c(0.5,10, 1), p3=c(0.5, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "birnbaumsaunders" = {     ## [R+, R+, R]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 10, 1), p2=c(0.5,5, 1), p3=c(0, 0, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "binomial" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa, size) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1/0.9), log(0.9/0.1)))$minimum
             }
           } ,


           "boxcox" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5,5, 1), p3=c(0.5, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "burr" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 10, 2), p2=c(0.5,5, 2), p3=c(1, 3, 2), p4=c(1,1,2),
                             p5=c(1,1,2), p6=c(1,1,2))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "burr2param" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5,5, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "cauchy" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 5, 1), p2=c(0.5, 10, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "chen" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "chi" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 5, 1), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                           p5=c(1,1,1), p6=c(1,1,1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "chisquared" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 5, 1), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                           p5=c(1,1,1), p6=c(1,1,1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "clg" = {     ## [R+, R+, R]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5,5, 1), p3=c(-1, 1, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "complementarybeta" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5,5, 1), p3=c(0, 0, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "dagum" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 1, 1), p2=c(0.5,5, 1), p3=c(0.5, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "diffzeta" = {     ## [R+, >1]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1.1,5, 1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "discretegamma" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1,5, 1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "discretelaplace" = {     ## [R, 01]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1,1, 1), p2=c(0.1,0.9,0.1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "discretenormal" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1,1, 1), p2=c(1,5,1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "discreteweibull" = {     ## [01, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.1,0.9, 1), p2=c(1,5,1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "doubleweibull" = {     ## [R+, R, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(-1,1, 1), p3=c(0.5, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "exponential" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 10, 1), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                           p5=c(1,1,1), p6=c(1,1,1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "exponentialextension" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "exponentialgeometric" = {     ## [R+, 01]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.1, 0.9, 0.2), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "exponentiallogarithmic" = {     ## [R+, 01]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(0.1, 0.9, 0.2), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "exponentialpoisson" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "exponentialpower" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "ev" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 1), p2=c(0.5, 10, 0.5), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "F" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "fellerpareto" = {     ## [R(mini), R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(min(y)-0.5, min(y)-0.5, 0.5),
                           p2=c(0.5, 5, 1), p3=c(0.5, 5, 1), p4=c(0.5, 5, 1), p5=c(0.5, 5, 1),
                           p6=c(1, 1, 1))
             theta000 = GridSearchS0(family, y0, params)
             theta00[1,1:4] = theta000[2:5]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "fisk" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 0.5), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "foldednormal" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(0.5, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "frechet" = {     ## [R+, R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 5, 1), p2=c(min(y)-0.5, min(y)-0.5, 5), p3=c(0.5, 5, 2),
                           p4=c(1, 1, 1), p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta000 = GridSearchS0(family, y0, params)
             theta00[1,] = c(theta000[1], theta000[3])
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00[1,]), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "hzeta" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(1, 10))$minimum
             }
           } ,


           "gamma" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "gammapoisson" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1,5, 1), p3=c(0.5, 0.5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "gaussian" = {     ## [R, R+]
             p = 2 ;
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for (j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               theta00[j,] = c(mean(x), stats::sd(x))
             }
             alpha0[,1]= theta00[,1]
             alpha0[,2]= log(theta00[,2])
           } ,


           "gev" = {     ## [R, R+, R]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               theta00[j,] = EnvStats::egevd(as.numeric(y))$parameters} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
               alpha0[j,] = theta2alpha(family, theta00[j,])
             }
           } ,


           "geninvbeta" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(0.5, 5, 0.5), p3=c(1, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "genlogis" = {     ## [R+, R, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(-1, 1, 1), p3=c(1, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "genlogis3" = {     ## [R+, R, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(-1, 1, 1), p3=c(1, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "genlogis4" = {     ## [R+, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1, 5, 1), p3=c(-1, 1, 1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "genpowerweibull" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(1, 5, 1), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "geometric" = {     ## [01]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1/0.9), log(0.9/0.1)))$minimum
             }
           } ,



           "generalizedhyperbolic" = {     ## [R, R+, R+, R, R]   [mu, delta, alpha, beta, lambda] (avec alpha^2 > beta^2)

             p = 5
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-2, 2, 1), p2=c(0.5, 5, 0.5), p3=c(1, 5, 1), p4=c(0.1, 0.9, 0.1),
                             p5=c(1, 5, 1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "generalizedlambda" = {     ## [R, R+, R, R]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-2, 2, 1), p2=c(0.5, 5, 0.5), p3=c(0, 3, 1), p4=c(0, 3, 1),
                             p5=c(1, 1, 1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "generalizedt" = {     ## [R, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-2, 2, 1), p2=c(0.5, 5, 0.5), p3=c(1, 5, 1), p4=c(1, 5, 1),
                             p5=c(1, 1, 1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "gompertz" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "gpd" = {     ## [R+(min(y)), R+, R]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(min(y), min(y), 1), p2=c(0.5, 5, 0.5), p3=c(0, 0, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "gumbel" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(0.5, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "gumbel2" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "halfcauchy" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "halflogistic" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,



           "halfnormal" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "halft" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "hjorth" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "hblaplace" = {     ## [01, R, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.1, 0.9, 0.1), p2=c(-2, 2, 1), p3=c(0.5,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "huber" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 1), p2=c(0.5, 10, 0.5), p3=c(1.345, 1.345, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "hyperbolic" = {     ## [R, R+, R+, R]   [mu, delta, alpha, beta] (avec alpha^2 > beta^2)

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 1), p2=c(0.5, 5, 0.5), p3=c(1, 2, 1), p4=c(0.1, 0.9, 0.1),
                             p5=c(1, 1, 1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inversebeta" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 0.5), p3=c(1, 1, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inverseburr" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 0.5), p3=c(1, 1, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inversechisquared" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "inverseexponential" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "inverseexpexponential" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inversegamma" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inverselomax" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inverseparalogistic" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inversepareto" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inversetransformedgamma" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 10, 0.5), p3=c(1, 1, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "inverseweibull" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamy" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamyexponential" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamygamma" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamygumbel" = {     ## [R+, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(-1,1,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamyhalfnormal" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamyloglogistic" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamynormal" = {     ## [R, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "kumaraswamyweibull" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,5,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "laplace" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "levy" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "linearfailurerate" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "lindley" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "libbynovickbeta" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1, 5, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "logcauchy" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "loggamma" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 1), p2=c(0.5, 10, 0.5), p3=c(1, 1, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "loggumbel" = {     ## [R, R+,]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "loglaplace" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 1), p2=c(0.5, 10, 0.5), p3=c(1, 5, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "loglog" = {     ## [R+, >1]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(1.5, 10, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "loglogistic" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "lognormal" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "lognormal3" = {     ## [R, R+, R]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-2, 2, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "logistic" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "logisticexponential" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "logisticrayleigh" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "logseries" = {     ## [01]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1/0.9), log(0.9/0.1)))$minimum
             }
           } ,


           "lomax" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "makeham" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 10, 0.5), p3=c(0.5, 5, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "maxwell" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 10, 0.5), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "mcgilllaplace" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "moexponential" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "moweibull" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "nakagami" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "ncchisquared" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "ncF" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 10, 0.5), p3=c(0.5, 5, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "negativebinomial" = {     ## [01]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa, size) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1/0.9), log(0.9/0.1)))$minimum
             }
           } ,


           "normalinversegaussian" = {     ## [R, R+, R+, R]   [mu, delta, alpha, beta] (avec alpha^2 > beta^2)

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 1), p2=c(1, 5, 0.5), p3=c(1, 5, 1), p4=c(0.1, 0.9, 0.1),
                             p5=c(1, 1, 1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "nsbeta" = {     ## [R+, R+, R(min), R(maxi)]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 5, 1), p2=c(0.5, 10, 0.5), p3=c(1,1,1),
                           p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1))
             theta000 = GridSearchS0(family, y0, params)
             theta0 = theta000[1:2]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "paralogistic" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "pareto" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "pareto1" = {     ## [R+, R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "pareto2" = {     ## [R, R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(min(y)-0.5, min(y)-0.5, 1), p2=c(0.5, 10, 0.5),
                           p3=c(0.5, 5, 0.5), p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1))
             thetaaa0 = GridSearchS0(family, y0, params)
             theta00[1,] = thetaaa0[2:3]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00[1,]), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "pareto3" = {     ## [R, R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(min(y)-0.5, min(y)-0.5, 1), p2=c(0.5, 10, 0.5),
                           p3=c(0.5, 5, 0.5), p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1))
             thetaaa0 = GridSearchS0(family, y0, params)
             theta00[1,] = thetaaa0[2:3]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00[1,]), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "pareto4" = {     ## [R, R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(min(y)-0.5, min(y)-0.5, 1), p2=c(0.5, 10, 0.5),
                           p3=c(0.5, 5, 0.5), p4=c(0.5,5,1), p5=c(1,1,1), p6=c(1,1,1))
             thetaaa0 = GridSearchS0(family, y0, params)
             theta00[1,] = thetaaa0[2:4]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00[1,]), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "paretopositivestable" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 0.5), p3=c(0.5, 5, 0.5),
                           p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "pctalaplace" = {     ## [R+, R]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(-1, 1, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "perks" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "poisson" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(1), log(100)))$minimum
             }
           } ,


           "power1" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "power2" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "powerdistribution" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 10, 0.5), p2=c(0.5, 10, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "powerexponential" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 2, 0.5), p2=c(0.5, 2, 0.5), p3=c(0.5,1,0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "rayleigh" = {     ## [R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(1, 10, 1), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "reflectedgamma" = {     ## [R+, R, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(-1, 1, 1), p3=c(1,5,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "rice" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "scaledchisquared" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "schabe" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "simplex" = {     ## [01, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.1, 0.9, 0.2), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "skewedlaplace" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(0.5, 5, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "skewedt" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(0.5, 0.5, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "skewedtfourparam" = {     ## [R, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(-1, 1, 1), p4=c(1, 3, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "skewednormal" = {     ## [R, R+, R, R]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(-1, 1, 1), p4=c(-2, 5, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "skewedgeneralizedt" = {     ## [R, R+, -1+1, R+(>1), R+(>1)]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(2, 5, 1), p3=c(-0.9, 0.9, 0.1), p4=c(1, 10, 1),
                             p5=c(100, 100, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "skewedexponentialpower" = {     ## [R, R+, R, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(-1, 1, 0.5), p2=c(0.5, 5, 1), p3=c(-1, 1, 1), p4=c(1, 5, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,




           "slash" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           # "stable" = {     ## [02, -1+1, R+, R]
           #
           #   p = 4
           #   theta00 = matrix(0,reg,p)
           #   alpha0 = theta00
           #   params = list(p1=c(0.5, 2, 0.5), p2=c(0, 0, 1), p3=c(1, 5, 1), p4=c(0, 0, 1),
           #                 p5=c(1, 1, 1), p6=c(1, 1, 1))
           #   theta00 = GridSearchS0(family, y0, params)
           #   for( j in un:reg){
           #     ind = (j-un)*n0+ind0
           #     x=y0[ind]
           #     funParam0 <- function(thetaa){
           #       log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
           #       return(log_likelihood)
           #     }
           #     alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
           #                                method = "Nelder-Mead")$par
           #   }
           # } ,


           "stacy" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1, 5, 1), p4=c(0, 0, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "t" = {     ## [R, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(0.5, 5, 1), p3=c(0.5, 5, 0.5), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,




           "tobit" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0, 0, 0.5), p2=c(1, 10, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "topple" = {     ## [01]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0.1, 0.9, 0.2), p2=c(1, 1, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                           p5=c(1, 1, 1), p6=c(1, 1, 1))
             theta00 = GridSearchS0(family, y0, params)
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
           } ,


           "transformedbeta" = {     ## [R+, R+, R+, R+]

             p = 4
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(0.5, 5, 1), p4=c(1, 1, 1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "transformedgamma" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 1), p3=c(1, 5, 1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,



           "truncatednormal" = {     ## [R, R+, R(min), R(max)]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(0, 0, 1), p2=c(0.5, 10, 0.5), p3=c(min(y)-0.05, min(y)-0.05, 1),
                           p4=c(max(y)+0.05, max(y)+0.05, 0.5), p5=c(1,1,1), p6=c(1,1,1))
             thetaaa0 = GridSearchS0(family, y0, params)
             theta0 = thetaaa0[1:2]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "truncatedpareto" = {     ## [R+(mini), R+(maxi), R+]

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             params = list(p1=c(min(y)-0.05, min(y)-0.05, 1), p2=c(max(y)+0.05, max(y)+0.05, 0.5),
                           p3=c(0.5, 10, 0.5), p4=c(0.5,0.5,1), p5=c(1,1,1), p6=c(1,1,1))
             thetaaa0 = GridSearchS0(family, y0, params)
             theta0 = thetaaa0[3]
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(0.1), log(100)))$minimum
             }
             if (alpha0[1]==alpha0[2]){
               alpha0[2] = log(exp(alpha0[1]+5))
             }
           } ,


           "twosidedpower" = {     ## [01, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.1, 0.9, 0.2), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "wald" = {     ## [R, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "weibull" = {     ## [R+, R+]

             p = 2
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 0.5), p2=c(0.5, 5, 1), p3=c(1,1,1), p4=c(1,1,1),
                             p5=c(1,1,1), p6=c(1,1,1))
               theta00 = GridSearchS0(family, y0, params)} else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "xie" = {     ## [R+, R+, R+]

             p = 3
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             if (start == 0){
               params = list(p1=c(0.5, 5, 1), p2=c(0.5, 5, 0.5), p3=c(1, 5, 1), p4=c(1, 1, 1),
                             p5=c(1, 1, 1), p6=c(1, 1, 1))
               theta00 = GridSearchS0(family, y0, params) } else{
                 theta00 = start
               }
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optim(par= theta2alpha(family,theta00), funParam0,
                                          method = "Nelder-Mead")$par
             }
           } ,


           "yules" = {     ## [R+] > 0.5

             p = 1
             theta00 = matrix(0,reg,p)
             alpha0 = theta00
             for( j in un:reg){
               ind = (j-un)*n0+ind0
               x=y0[ind]
               funParam0 <- function(thetaa){
                 log_likelihood = (-sum(log(PDF_unc(family, x, thetaa) ) ) )
                 return(log_likelihood)
               }
               alpha0[j,] =  stats::optimize(funParam0, c(log(1), log(100)))$minimum
             }
           } ,

    )

  } else {
      alpha0 = matrix(0, reg, p)
      for (j in un:reg){
        alpha0[j,] = theta2alpha(family,theta00[j,])
      }
  }




  Q0 = matrix(1,reg,reg)/reg


  # warm-up
  for (k in 1:ninit){
    emstep= EMStep(y,ZI,alpha0,Q0,family, size)
    nu=emstep$nu ;
    alpha_new = emstep$alpha_new;
    Qnew = emstep$Qnew;
    eta = emstep$eta
    eta_bar =emstep$eta_bar;
    lambda= emstep$lambda;
    Lambda=emstep$Lambda ;
    LL=emstep$LL
    Q0 = Qnew
    alpha0 = alpha_new
  }




  # iterations
  for (k in 1:max_iter){
    emstep= EMStep(y,ZI,alpha0,Q0,family, size)
    nu = emstep$nu ;
    alpha_new = emstep$alpha_new;
    Qnew = emstep$Qnew;
    eta = emstep$eta
    eta_bar = emstep$eta_bar;
    lambda= emstep$lambda;
    Lambda = emstep$Lambda ;
    LL=emstep$LL

    sum1 = sum(abs(alpha0))
    sum2 = sum(abs(alpha_new-alpha0))
    if ( sum2 < (sum1 * reg * eps) ){
      break
    }
    Q0 = Qnew
    alpha0 = alpha_new
  }

  # output
  alpha = alpha0
  Q = Q0


  theta = matrix(0,reg,p)
  for (j in un:reg){
    theta[j,] = alpha2theta(family, alpha[j,])
  }


  numParam = (dim(theta)[1]-ZI)*dim(theta)[2] + reg*(reg-1)

  AIC = (2 * numParam - 2*LL  ) / n

  BIC = (log(n) * numParam - 2*LL  ) / n

  CAIC = ( (log(n)+1) * numParam - 2*LL  ) / n

  AICcorrected = AIC + (2*(numParam)*(numParam+1))/(n-numParam-1)

  HQC = (2*numParam*log(log(n)) - 2*LL) / n





  tryCatch({
             u_Ros = stats::runif(n)
             cdf_gof = matrix(0,n,reg)
             for(j in un:reg){
                 cdf_gof[,j] = CDF_est(family, y, theta[j,], size = size, u_Ros = u_Ros)
             }
             if(ZI==1) cdf_gof[,1] = as.numeric(y>0)+u_Ros*(y==0)

             eta00 = rep(1,reg)/reg
             w00 = rbind(eta00,eta) %*% Q
             W   = w00[1:(dim(w00)[1]-1),]
             if (reg == 1){
                W = data.frame(as.list(W))
             }
             if (reg>1){
               U = rowSums(W * cdf_gof)
             } else {
                 U = cdf_gof
                 W = eta
             }
             cvm = Snd1(U)
           },
           error = function(e) {
             message("Difficulties computing the CDF, results for the CDF are not valid")
             U = 0
             W = 0
             cvm = 1000
           }
          )


  pred_e = rep(0,n);
  pred_l = rep(0,n);
  pred_e = rep(0,n);
  pred_l = rep(0,n);
  for (i in 1:n){
    bb_e = sort(eta[i,], index.return=TRUE)$ix
    pred_e[i] = bb_e[reg]
    bb_l = sort(lambda[i,], index.return=TRUE)$ix
    pred_l[i] = bb_l[reg]
  }
  dp_e = diff(pred_e);
  dp_l = diff(pred_l);
  runs_e = 1+sum(dp_e != 0);
  runs_l = 1+sum(dp_l != 0);

  stats = matrix(0, nrow = reg, ncol = 2)
  for (i in 1:reg){
    stats[i,] = c(mean(y[pred_l==i]), stats::sd(y[pred_l==i]))
  }


  if (graph){
    graphs = graphEstim(y, ZI, reg, theta, family, pred_l, pred_e)
  }



  out = list(theta=theta, Q=Q, eta=eta, nu=nu, U=U, cvm=cvm, W=W,lambda=lambda, LL=LL,
             AIC=AIC, BIC=BIC, AICcorrected=AICcorrected, CAIC=CAIC, HQC=HQC, stats=stats,
             pred_l=pred_l, pred_e=pred_e, runs_e=runs_e, runs_l=runs_l)




  return(out)
}




EMStep <- function(y, ZI, alpha, Q, family, size){

  n = length(y) #length of series
  r = dim(alpha)[1] #number of regimes
  p = dim(alpha)[2]
  eta_bar = matrix(0,n,r)
  eta = matrix(0,n,r)
  lambda   = matrix(0,n,r)
  c   = matrix(0,n,r)
  Lambda   = array(0,c(r,r,n))
  Z = rep(0,n)


  un=1+ZI

  if(ZI==1)
  {
    ind0=(y==0)
    ind=(y!=0)
    c[ind0,1] = 1


    for (j in un:r){

      c[ind,j] = PDF_unc(family, y[ind], alpha[j,], size)
    }
  }else{
    ind=c(1:n)
    for (j in 1:r){
      c[,j] = PDF_unc(family, y, alpha[j,], size)  }
  }

  # eta_bar
  eta_bar[n,]=1/r

  for(k in 1:(n-1)){
    i = n-k        # backward
    j = i+1        # j is the index of period (t+1)
    v =  ( eta_bar[j,] * c[j,] ) %*% t(Q) #  numerator of eta_bar(i)
    eta_bar[i,] = v/sum(v);
  }


  # eta
  eta0 = rep(1,r)/r

  v = ( eta0 %*% Q) * c[1,]  #numerator of eta at  t = 1;
  eta[1,] = v/sum(v)

  for (i in 2:n){
    v = ( eta[i-1,] %*% Q) * c[i,]    # numerator of eta when t > 1;
    Z[i] = sum(v)
    eta[i,] = v/Z[i]        #value of eta when  t > 1;
  }
  LL = sum(log(Z)[2:length(Z)])


  #lambda
  v      = eta * eta_bar    # numerator of lambda
  sv0    = rowSums(v)       # computation of the denomirator of lambda

  for(j in 1:r){
    lambda[,j] = v[,j] / sv0        # computation of lambda
  }

  if (r==1){
    Q = matrix(1,1,1)
  }

  #lambda
  gc =  eta_bar * c   #les deux derniers calcul (plus \E0 droite) du numerateur de Lambda

  M = Q * (as.matrix(eta0)%*% gc[1,] ) # numerator of Lambda
  MM = sum(M)   #denominator of Lambda
  Lambda[,,1] = M/MM


  for(i in 2:n){
    M = Q * ( as.matrix(eta[i-1,]) %*% gc[i,]  )    # numerateur de Lambda
    MM = sum(M)    #denominateur de Lambda
    Lambda[,,i] = M/MM    #Lambda pour t < n
  }
  nu = colMeans(lambda)
  Qnew = Q

  if (r >= 2){
    for(j in 1:r){
      sv = rowSums(Lambda[j,,], dims=1)
      ssv = sum(sv)
      Qnew[j,] = sv/ssv
    }
  }

  y0=y[ind]

  alpha_new = matrix(0, nrow=r, ncol=p)
  for(i in un:r){

    fun <- function(alphaa){
      log_likelihood = (-sum(lambda[ind,i]* log(PDF_unc(family, y0, alphaa, size) ) ) )
      return(log_likelihood)
    }
    if (p>1){
      alpha_new[i,] =  stats::optim(par= alpha[i,], fn = fun, method = "Nelder-Mead")$par
    } else {
      alpha_new[i,] =  stats::optimize(fun, c(log(0.1), log(100)))$minimum
    }

  }


  out = list(nu=nu, alpha_new=alpha_new, Qnew=Qnew, eta=eta, eta_bar=eta_bar, lambda=lambda,
             Lambda=Lambda, LL=LL)
  return(out)


}
