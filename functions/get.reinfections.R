

get.reinfections <- function(params, ...) {
    
    if(class(params) == "data.frame") params <- unlist(params)
    
    full.sol <- model_solution(params, full_solution=T, ...)
    
    T <- nrow(full.sol)
    C1f <- sum(full.sol[T, c(paste0("POP.C1", 1:3))] / SEVERITY.FRAC)
    C2f <- sum(full.sol[T, c(paste0("POP.C2", 1:3))] / SEVERITY.FRAC)
    T1i <- sum(full.sol[1, c(paste0("POP.A1", 1:3))],
               full.sol[1, c(paste0("POP.I1", 1:3))],
               full.sol[1, c(paste0("POP.H1", 1:3))],
               full.sol[1, c(paste0("POP.R1", 1:3))])
    T1f <- sum(full.sol[T, c(paste0("POP.A1", 1:3))],
               full.sol[T, c(paste0("POP.I1", 1:3))],
               full.sol[T, c(paste0("POP.H1", 1:3))],
               full.sol[T, c(paste0("POP.R1", 1:3))])
    R <- C1f - (T1f - T1i)
    return(c(reinf = R, R1 = T1i, C1 = C1f, C2 = C2f))
}
