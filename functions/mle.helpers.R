## convert parameters using logit link functions for those bound between zero
## and one.  That is, logit link function for INIT.V2.FRAC, PROB.REINFEC.V2,
## and prevalence
conv.logit <- function(params) {
    params[c(1,3,4)] <- inv.logit(params[c(1,3,4)])
    params[4] <- 100*params[4] ## the model is fitted with prevalence in 0-100 scale
    params
}

# 5 pars, original scale
helper.mle.fun <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r) {
    l <- as.list(match.call())
    params <- unlist(l[2:length(l)])
    return(neg_loglike_sol(model_solution(params)))
}
# 5 pars, logit scale
helper.mle.fun.2 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r) {
    l <- as.list(match.call())
    params <- conv.logit(unlist(l[2:length(l)]))
    return(neg_loglike_sol(model_solution(params)))
}
# 5 pars, IHR.V2.PROP x 2, original scale
helper.mle.fun.3 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r) {
    l <- as.list(match.call())
    params <- c(unlist(l[2:length(l)]), IHR.V2.PROP = 2)
    return(neg_loglike_sol(model_solution(params)))
}
# 5 pars, IHR.V2.PROP x 2, logit scale
helper.mle.fun.4 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r) {
    l <- as.list(match.call())
    params <- c(conv.logit(unlist(l[2:length(l)])), IHR.V2.PROP = 2)
    return(neg_loglike_sol(model_solution(params)))
}
# 6 pars, original scale
helper.mle.fun.5 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r, IHR.V2.PROP) {
    l <- as.list(match.call())
    params <- unlist(l[2:length(l)])
    return(neg_loglike_sol(model_solution(params)))
}
# 6 pars, logit scale
helper.mle.fun.6 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                           prevalence, r, IHR.V2.PROP) {
    l <- as.list(match.call())
    params <- conv.logit(unlist(l[2:length(l)]))
    return(neg_loglike_sol(model_solution(params)))
}

# 4 pars, logit scale, fixed prevalence
helper.mle.fun.7 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2, r)
{
    l <- as.list(match.call())
    params <- unlist(l[2:length(l)])
    params[c(1,3)] <- inv.logit(params[c(1,3)])
    params <- c(params, prevalence = 21.7)
    return(neg_loglike_sol(model_solution(params)))
}

