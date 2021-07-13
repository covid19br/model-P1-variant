#' Fits an exponential model for n cases as a function of n days elapsed since the start of the series.
#' @details Fits a mixed generalized model for counts (Poisson) of the values of a time series as a 
#' function of the number and days elapsed since the beginning of the series. As the Poisson model 
#' has a log link function, it is equivalent to adjusting an exponential growth function to counts 
#' with Poisson errors.
#' @param zoo.obj zoo class object with univariate time series
#' @param family name of the error distribution to use in the linear model
#'     (see function help at stats::family)
#' @param only.coef if TRUE the function returns a vector with the regression 
#' coefficients and their confidence intervals. If FALSE returns the adjusted 
#' class glm model object

if (!require(zoo)) install.packages(zoo); library(zoo)
fitP.exp <- function(zoo.obj, only.coef = TRUE){
    ## Does not work with roll aply function.
    ## if(class(zoo.obj)!="zoo"|!is.null(ncol(zoo.obj)))
    ##    stop("'zoo.obj' must be an object of class zoo with a single variable")
    ndias <- as.vector( rev(max(as.Date(time(zoo.obj))) - as.Date(time(zoo.obj)) ))
    fit <- try(glm(as.integer(zoo.obj) ~ ndias, family = poisson))
    if(only.coef){
            ci <- try(confint(fit))
            if((any(class(fit)=="try-error")||any(class(ci)=="try-error")))
                results <- rep(NA,6)
            else
                results <- c(coef(fit),ci[1,], ci[2,])
            names(results) <- c("intercept", "coef", "int.low", "int.upp", "coef.low", "coef.upp")
    }
    if(!only.coef){
        if(any(class(fit)=="try-error"))
            results <- NA
        else
            results  <-  fit
    }
    return(results)
}
