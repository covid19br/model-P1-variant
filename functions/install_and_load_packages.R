#' install.and.load.packages
#' 
#' A function that evaluates if a package is already installed
#' and install it otherwise. Next, it loads the package to the R environment.
#'
#' @param package.name character. Name of the package that will be installed and required.

install.and.load.packages <- function(package.name) {
  
  if(class(package.name) != "character") warning("An object of class 'character' must be provided")
  
  if (!require(package.name, character.only = TRUE)) install.packages(package.name, character.only = TRUE)
  library(package.name, character.only = TRUE)
}