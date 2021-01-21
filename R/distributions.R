#'@title The names and descriptions of the univariate distributions
#'
#'
#'@description This function allows the users to find the details on the available distributions.
#'
#'
#'
#'@return No returned value, allows the users to know the differents distributions and parameters
#'
#'@export
#'
#'


distributions <- function(){

  distr = utils::read.csv("https://www.dropbox.com/s/dxmo875sg2jq502/Distribution.csv?dl=1")
  utils::View(distr)

}
