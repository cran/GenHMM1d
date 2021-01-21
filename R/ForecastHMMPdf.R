#'@title Forecasted density function of a univariate HMM at time n+k1, n+k2, ...
#'
#'@description This function computes the probability forecasted density function of a univariate HMM for multiple horizons, given observations up to time n
#'
#'@param y  points at which the pdf function is computed
#'@param family    distribution name; run the function distributions() for help
#'@param theta  parameters; (r  x p)
#'@param Q    probability transition  matrix; (r  x r)
#'@param eta    vector of the estimated probability of each regime at time n; (1  x r)
#'@param k  prediction times (may be a vector of integers).
#'@param graph (0 or else) produce plots
#'
#'@return \item{pdf}{values of the pdf function}
#'
#'
#'@examples
#'family = "gaussian"
#'
#'lb = -6
#'ub = 6
#'
#'theta = matrix(c(-1.5, 1.7, 1, 1),2,2)
#'Q = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
#'eta = c(0.96091218, 0.03908782)
#'
#'
#'forecastedhmmpdf = ForecastHMMPdf(y=seq(from=lb, to=ub, by=0.1), family=family,
#' theta=theta, Q=Q, eta=eta, k=c(1,5,10,20), graph=1)
#'
#'
#'@export


ForecastHMMPdf<-function(y, family, theta, Q, eta, k=1, graph=0){

  if(is.null(dim(Q))){
    QQ0 = matrix(Q)
    r = dim(QQ0)[1]
  } else {
    r = dim(Q)[2]
  }


  pdf = matrix(0, nrow=length(y), ncol=length(k))

  for (d in 1:length(k)){
    Q_prime = matrixcalc::matrix.power(Q, as.integer(k[d]))
    for (l in 1:r){
      for (j in 1:r){
        pdf[,d] = pdf[,d] + eta[j] * Q_prime[j, l] * PDF(family, y, theta[l,])
      }
    }
  }


  if (graph != 0){
    df = data.frame(pdf)
    colnames(df) = as.character(k)
    df$index = y
    df = reshape2::melt(df ,  id.vars = 'index', variable.name = 'horizon')
    print(
      ggplot2::ggplot(df, ggplot2::aes(df$index,df$value)) + ggplot2::geom_line(ggplot2::aes(colour = df$horizon)) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(fill = "white"),
              panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
              legend.position = c(.02, .95),
              legend.justification = c("left", "top"),
              legend.margin = ggplot2::margin(6, 6, 6, 6),
              legend.box.background = ggplot2::element_rect(color="black", size=1),
              legend.box.margin = ggplot2::margin(1, 1, 1, 1),
              plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Forecasted probability density function", x="observations", y="densities" )
    )
  }


  return(pdf)


}


