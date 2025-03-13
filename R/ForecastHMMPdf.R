#'@title Forecasted density function of a univariate HMM at time n+k1, n+k2, ...
#'
#'@description This function computes the probability forecasted density function (with respect to Dirac(0)+Lesbesgue) of a univariate HMM for multiple horizons, given observations up to time n
#'
#'@param y          points at which the pdf function is computed
#'@param ZI         1 if zero-inflated, 0 otherwise (default)
#'@param family     distribution name; run the function distributions() for help
#'@param theta      parameters; (r  x p)
#'@param Q          probability transition  matrix for the regimes; (r  x r)
#'@param eta        vector of the estimated probability of each regime at time n; (1  x r)
#'@param k          prediction times (may be a vector of integers)
#'@param size       additional parameter for some discrete distributions; run the command distributions() for help
#'@param graph      TRUE to produce plots (FALSE is default)
#'
#'@return \item{pdf}{values of the pdf function}
#'
#'
#'@examples
#'family = "gaussian"
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.06, 0.94)
#'x=seq(from=-6, to=6, by=0.1)
#' k=c(1,5,10,20)
#'pdf = ForecastHMMPdf(x, 1, family, theta, Q, eta, k=k, graph=TRUE)
#'
#'
#'@export


ForecastHMMPdf<-function(y, ZI=0, family, theta, Q, eta, size=0, k=1, graph=FALSE){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }

  un = 1+ZI;

  pdf = matrix(0, nrow=length(y), ncol=length(k))

  for (d in 1:length(k)){
    Q_prime = matrixcalc::matrix.power(Q, as.integer(k[d]))
    for (l in un:r){
      for (j in 1:r){
        pdf[,d] = pdf[,d] + eta[j] * Q_prime[j, l] * PDF(family, y, theta[l,],size)
      }
    }
    if(ZI==1)
    {
      for (j in 1:r){
        pdf[,d] = pdf[,d] + eta[j] * Q_prime[j, 1] * (y==0)
      }
    }
  }


  if (graph){
    Horizon=NULL
    index= NULL
    value=NULL
    df = data.frame(pdf)
    colnames(df) = as.character(k)
    df$index = y
    df = reshape2::melt(df ,  id.vars = 'index', variable.name = 'Horizon')
    g=ggplot2::ggplot(df, ggplot2::aes(index,value)) + ggplot2::geom_line(ggplot2::aes(colour = Horizon)) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(fill = "white"),
              panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
              legend.position = c(.02, .95),
              legend.justification = c("left", "top"),
              legend.margin = ggplot2::margin(6, 6, 6, 6),
              legend.box.background = ggplot2::element_rect(color="black", size=1),
              legend.box.margin = ggplot2::margin(1, 1, 1, 1),
              plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Forecasted probability density functions", x="values", y="densities" )

       print(g)
  }


  return(pdf)


}


