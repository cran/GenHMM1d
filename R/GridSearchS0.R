#'@title Gridsearch
#'
#'@description This function performs a gridsearch to find a good starting value for the EM algorithm. A good starting value for the EM algorithm is one for which all observations have strictly positive density (the higher the better)
#'
#'@param family distribution name; run the function distributions() for help
#'@param y  observations
#'@param params list of six vectors named (p1, p2, p3, p4, p5, p6). Each corresponding to a parameter of the distribution (additionnal parameters will be ignored).
#' For example : params = list(p1=c(0.5, 5, 0.5), p2=c(1, 5, 1), p3=c(0.1, 0.9, 0.1), p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1)) where p1 is the grid of value for the first parameter.
#'@param size additional parameter for some discrete distributions; run the command distributions() for help
#'@param lbpdf  minimal acceptable value of the density; (should be >= 0)
#'
#'@return \item{goodStart}{accepted parameter set}
#'
#'
#'@examples
#'family = "gaussian"
#'
#' Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2) ;
#' theta = matrix(c(-1.5, 1.7, 1, 1),2,2) ;
#' sim = SimHMMGen(theta, size=0, Q, ZI=0,"gaussian", 50)$SimData ;
#' params = list(p1=c(-2, 2, 0.5), p2=c(1, 5, 1), p3=c(1, 1, 1), p4=c(1,1,1), p5=c(1,1,1), p6=c(1,1,1))
#' accepted_params = GridSearchS0(family, sim, params)
#'
#'
#'@export


GridSearchS0<-function(family, y, params, size=0, lbpdf=0){


    param1 = seq(from = params$p1[1], to = params$p1[2], by=params$p1[3])

    param2 = seq(from = params$p2[1], to = params$p2[2], by=params$p2[3])

    param3 = seq(from = params$p3[1], to = params$p3[2], by=params$p3[3])

    param4 = seq(from = params$p4[1], to = params$p4[2], by=params$p4[3])

    param5 = seq(from = params$p5[1], to = params$p5[2], by=params$p5[3])

    param6 = seq(from = params$p6[1], to = params$p6[2], by=params$p6[3])



    paramtest1 = c(param1[length(param1)], param2[1], param3[1], param4[1], param5[1], param6[1])
    f = PDF(family, y, paramtest1,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest1
      return(goodStart)
    }
    paramtest2 = c(param1[1], param2[length(param2)], param3[1], param4[1], param5[1], param6[1])
    f = PDF(family, y, paramtest2,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest2
      return(goodStart)
    }
    paramtest3 = c(param1[1], param2[1], param3[length(param3)], param4[1], param5[1], param6[1])
    f = PDF(family, y, paramtest3,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest3
      return(goodStart)
    }
    paramtest4 = c(param1[1], param2[1], param3[1], param4[length(param4)], param5[1], param6[1])
    f = PDF(family, y, paramtest4,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest4
      return(goodStart)
    }
    paramtest5 = c(param1[1], param2[1], param3[1], param4[1], param5[length(param5)], param6[1])
    f = PDF(family, y, paramtest5,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest5
      return(goodStart)
    }
    paramtest6 = c(param1[1], param2[1], param3[1], param4[1], param5[1], param6[length(param6)])
    f = PDF(family, y, paramtest6,size)
    if (!is.na(min(f)) && min(f) > lbpdf){
      goodStart = paramtest6
      return(goodStart)
    }



    for (i in 1:length(param1)){
        for (j in 1:length(param2)){
            for (k in 1:length(param3)){
                for (l in 1:length(param4)){
                    for (m in 1:length(param5)){
                        for (n in 1:length(param6)){
                            param = c(param1[i], param2[j], param3[k], param4[l], param5[m], param6[n])
                            f = PDF(family, y, param,size)
                            if (!is.na(min(f)) && min(f) > lbpdf){
                               goodStart = param
                               return(goodStart)
                            }
                        }
                    }
                }
            }
        }
    }







}




