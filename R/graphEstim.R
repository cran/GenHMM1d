#'@title Graphs
#'
#'@description This function shows the graphs resulting from the estimation of a HMM model
#'
#'@param y   observations
#'@param ZI  1 if zero-inflated, 0 otherwise (default)
#'@param reg    number of regimes
#'@param theta estimated parameters; (r x p)
#'@param family  distribution name; run the function distributions() for help
#'@param pred_l estimated regime using lambda
#'@param pred_e estimated regime using eta
#'
#'@return No returned value; produces figures of interest for the HMM model
#'
#'
#'@export
#'
graphEstim<-function(y, ZI=0, reg, theta, family, pred_l, pred_e){
  ninit=100   #minimum number of iterations
  n = length(y)

  un = 1+ZI;
    dd1 = data.frame("observations" = y)
    dd1$pred_l = as.factor(pred_l)
    dd1$pred_e = as.factor(pred_e)
    f = matrix(0, nrow = length(y), ncol = reg)
    for (j in un:reg){
      f[,j] = PDF(family, y, theta[j,])
    }

    if (reg==1){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$index = c(1:length(y))
    } else if (reg==2){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$index = c(1:length(y))
    } else if (reg==3){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$index = c(1:length(y))
    } else if (reg==4){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$index = c(1:length(y))
    } else if (reg==5){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$index = c(1:length(y))
    } else if (reg==6){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$index = c(1:length(y))
    } else if (reg==7){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$index = c(1:length(y))
    } else if (reg==8){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$index = c(1:length(y))
    } else if (reg==9){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$index = c(1:length(y))
    } else if (reg==10){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$index = c(1:length(y))
    } else if (reg==11){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$index = c(1:length(y))
    } else if (reg==12){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$regime12 = f[,12]
      dd1$index = c(1:length(y))
    } else if (reg==13){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$regime12 = f[,12]
      dd1$regime13 = f[,13]
      dd1$index = c(1:length(y))
    } else if (reg==14){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$regime12 = f[,12]
      dd1$regime13 = f[,13]
      dd1$regime14 = f[,14]
      dd1$index = c(1:length(y))
    } else if (reg==15){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$regime12 = f[,12]
      dd1$regime13 = f[,13]
      dd1$regime14 = f[,14]
      dd1$regime15 = f[,15]
      dd1$index = c(1:length(y))
    } else if (reg==16){
      d = stats::density(y)
      dens = data.frame("kernel"= d$y)
      dens$x = d$x
      dd1$regime1 = f[,1]
      dd1$regime2 = f[,2]
      dd1$regime3 = f[,3]
      dd1$regime4 = f[,4]
      dd1$regime5 = f[,5]
      dd1$regime6 = f[,6]
      dd1$regime7 = f[,7]
      dd1$regime8 = f[,8]
      dd1$regime9 = f[,9]
      dd1$regime10 = f[,10]
      dd1$regime11 = f[,11]
      dd1$regime12 = f[,12]
      dd1$regime13 = f[,13]
      dd1$regime14 = f[,14]
      dd1$regime15 = f[,15]
      dd1$regime16 = f[,16]
      dd1$index = c(1:length(y))
    }



    if (reg==1){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = "darkblue")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==2){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red')) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==3){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==4){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==5){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==6){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==7){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==8){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==9){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==10){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==11){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==12){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime12, color="regime 12"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue",
          "regime 12" = "lightblue3")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==13){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime12, color="regime 12"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime13, color="regime 13"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue",
          "regime 12" = "lightblue3",
          "regime 13" = "hotpink2")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==14){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime12, color="regime 12"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime13, color="regime 13"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime14, color="regime 14"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue",
          "regime 12" = "lightblue3",
          "regime 13" = "hotpink2",
          "regime 14" = "lightpink2")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==15){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime12, color="regime 12"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime13, color="regime 13"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime14, color="regime 14"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime15, color="regime 15"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue",
          "regime 12" = "lightblue3",
          "regime 13" = "hotpink2",
          "regime 14" = "lightpink2",
          "regime 15" = "lightsalmon2")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    } else if (reg==16){
      b = ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = dens$x, y = dens$kernel, color="kernel"), size=1.1, linetype="solid", data=dens) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime1, color="regime 1"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime2, color="regime 2"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime3, color="regime 3"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime4, color="regime 4"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime5, color="regime 5"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime6, color="regime 6"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime7, color="regime 7"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime8, color="regime 8"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime9, color="regime 9"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime10, color="regime 10"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime11, color="regime 11"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime12, color="regime 12"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime13, color="regime 13"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime14, color="regime 14"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime15, color="regime 15"), data=dd1) +
        ggplot2::geom_line(ggplot2::aes(x = dd1$observations, y = dd1$regime16, color="regime 16"), data=dd1) +
        ggplot2::scale_color_manual(values = c(
          "kernel" = "black",
          "regime 1" = 'darkblue',
          "regime 2" = 'red',
          "regime 3" = "orange",
          "regime 4" = "darkolivegreen4",
          "regime 5" = "darkorchid4",
          "regime 6" = "darkmagenta",
          "regime 7" = "aquamarine2",
          "regime 8" = "blue3",
          "regime 9" = "coral3",
          "regime 10" = "burlywood3",
          "regime 11" = "deepskyblue",
          "regime 12" = "lightblue3",
          "regime 13" = "hotpink2",
          "regime 14" = "lightpink2",
          "regime 15" = "lightsalmon2",
          "regime 16" = "rosybrown2")) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                       legend.position = c(.02, .95),
                       legend.justification = c("left", "top"),
                       legend.margin = ggplot2::margin(6, 6, 6, 6),
                       legend.box.background = ggplot2::element_rect(color="black", size=1),
                       legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                       legend.title = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title="Estimated densities" )
    }

    print(b)





    c = ggplot2::ggplot(dd1, ggplot2::aes(x = dd1$index, y=dd1$observations, color=pred_e)) +
      ggplot2::geom_point(size=2) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                     legend.position = c(.02, .95),
                     legend.justification = c("left", "top"),
                     legend.margin = ggplot2::margin(6, 6, 6, 6),
                     legend.box.background = ggplot2::element_rect(color="black", size=1),
                     legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                     legend.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ eta ), x="index" )
    print(c)



    d = ggplot2::ggplot(dd1, ggplot2::aes(x = dd1$index, y=dd1$observations, color=pred_l)) +
      ggplot2::geom_point(size=2) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                     legend.position = c(.02, .95),
                     legend.justification = c("left", "top"),
                     legend.margin = ggplot2::margin(6, 6, 6, 6),
                     legend.box.background = ggplot2::element_rect(color="black", size=1),
                     legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                     legend.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ lambda ), x="index" )
    print(d)



    print(ggplot2::ggplot(dd1, ggplot2::aes(pred_e, dd1$observations, color=pred_e)) +
            ggplot2::geom_boxplot() +
            ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2)) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = "white"),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                           legend.position = c(.02, .95),
                           legend.justification = c("left", "top"),
                           legend.margin = ggplot2::margin(6, 6, 6, 6),
                           legend.box.background = ggplot2::element_rect(color="black", size=1),
                           legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                           legend.title = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ eta ), x="regimes" )
    )



    print(ggplot2::ggplot(dd1, ggplot2::aes(pred_l, dd1$observations, color=pred_l)) +
            ggplot2::geom_boxplot() +
            ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2)) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = "white"),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                           legend.position = c(.02, .95),
                           legend.justification = c("left", "top"),
                           legend.margin = ggplot2::margin(6, 6, 6, 6),
                           legend.box.background = ggplot2::element_rect(color="black", size=1),
                           legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                           legend.title = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ lambda ), x="regimes")
    )



    print(ggplot2::ggplot(dd1, ggplot2::aes(x = dd1$index, y = pred_e, color = pred_e)) +
            ggplot2::geom_point() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = "white"),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                           legend.position = c(.02, .95),
                           legend.justification = c("left", "top"),
                           legend.margin = ggplot2::margin(6, 6, 6, 6),
                           legend.box.background = ggplot2::element_rect(color="black", size=1),
                           legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                           legend.title = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ eta ), x="index", y="regimes")
    )

    print(ggplot2::ggplot(dd1, ggplot2::aes(x = dd1$index, y = pred_l, color = pred_l)) +
            ggplot2::geom_point() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = "white"),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                           legend.position = c(.02, .95),
                           legend.justification = c("left", "top"),
                           legend.margin = ggplot2::margin(6, 6, 6, 6),
                           legend.box.background = ggplot2::element_rect(color="black", size=1),
                           legend.box.margin= ggplot2::margin(1, 1, 1, 1),
                           legend.title = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::labs(title=bquote("Most likely regime of each observation using" ~ lambda ), x="index", y="regimes")
    )

  #return(1)
}




