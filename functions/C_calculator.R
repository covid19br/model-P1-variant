if (!require(pracma)) install.packages("pracma"); library(pracma)

C_calculator <- function(Sol, EXPOSURE.PERIOD.DAYS.V2, SEVERITY.FRAC.V2){
    TIMES <- Sol[,1]
    POP.E1 <- Sol[,c(5:7)]
    POP.E2 <- Sol[,c(26:28)]
    
    new.Sol <- data.frame(Sol)
    POP.C1 <- data.frame(matrix(NA, ncol = 3, nrow = length(TIMES)))
    POP.C2 <- data.frame(matrix(NA, ncol = 3, nrow = length(TIMES)))
    
    for(i in 1:3){
        POP.C1[, i] <- cumtrapz(TIMES, POP.E1[, i]*SEVERITY.FRAC.V2[i]/EXPOSURE.PERIOD.DAYS.V2)
        POP.C2[, i] <- cumtrapz(TIMES, POP.E2[, i]*SEVERITY.FRAC.V2[i]/EXPOSURE.PERIOD.DAYS.V2)
    }
    new.Sol[, c(17:19)] <- POP.C1
    new.Sol[, c(38:40)] <- POP.C2
    return(new.Sol)
}