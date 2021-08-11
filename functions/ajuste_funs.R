## Load and install packages if required
source('functions/install_and_load_packages.R')
list.of.packages <- c("doParallel", "dplyr", "lubridate", "minpack.lm", 
                      "zoo", "rARPACK", "pse", "bbmle", "gridExtra", "optimx", "zeallot", "deSolve")
for(i in list.of.packages) install.and.load.packages(i)

## Load sources

source('functions/beta_init_condits2.R')
source('functions/logit.R')
source('functions/end.of.epiweek.R')
source('functions/C_calculator.R')
source('functions/SEIR_run.R')

#' Implement of numeric solver ODE from devtools into the 2-Strains model
#'
#' @param params A list of parameters to be used to solve the model 'params'. Default taken from params_epi.R
#' @param ... Any other parameters necessary to a different fitting '...'.
#'
#' @return A data.frame of Cumulative case of each strain.
#' 

model_solution <- function(params, full_solution = FALSE, return_compartments = FALSE, ...) {

    if ("prevalence" %in% names(params) | "r" %in% names(params)) {
        if ("prevalence" %in% names(params))
            PREVALENCE <- unlist(params["prevalence"]) / 100 * c(1,1,1)
        if ("r" %in% names(params))
            r <- params["r"]
        # set initial conditions
        init.conds <- init_condits(r, new.hosp, PREVALENCE = PREVALENCE, POP.DISTR,
                     CONTACT.M, EXPOSURE.PERIOD.DAYS, SICKNESS.PERIOD.DAYS,
                     SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC, SEVERE.CONT.REDUC.FRAC,
                     REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC, SEVERITY.FRAC, DEATH.FRAC,
                     V2.FRAC = 0)
        POP0 <- init.conds$POP0
        # Average contact rate and success between unvaccinated susceptibles and
        # infectious individuals \beta
        # calculated from hospitalized cases in the last 3 weeks
        BETA.RATE1 <- init.conds$BETA.RATE
        names(POP0) <- paste0(rep(c("POP.S", "POP.E1", "POP.A1", "POP.I1", "POP.H1",
                                     "POP.C1", "POP.R1", "POP.D1",
                                     "POP.E2", "POP.A2", "POP.I2", "POP.H2",
                                     "POP.C2", "POP.R2", "POP.D2"), each=3), rep(1:3, 15))
    }
    if ("IHR.V2.PROP" %in% names(params)) {
        ihr2 <- ihr[,2]/100
        ihr2 <- 1 / (1 + 1/(params["IHR.V2.PROP"] * ihr2/(1-ihr2)))
        weighted_ihr2 <- ihr2 * age.distr
        SEVERITY.FRAC.V2 = c(sum(weighted_ihr2[1:4]) / sum(age.distr[1:4]),
                             sum(weighted_ihr2[5:12]) / sum(age.distr[5:12]),
                             sum(weighted_ihr2[13:length(age.distr)]) / sum(age.distr[13:length(age.distr)]))
    }

    # introduce new variant in initial conditions
    new.POP0 <- POP0
    new.POP0[4:15] <- POP0[4:15] * (1-params["INIT.V2.FRAC"])
    new.POP0[25:36] <- POP0[4:15] * params["INIT.V2.FRAC"]

    diffEqs <- func.factory(REL_BETA_RATE2 = params["REL.BETA.RATE2"],
                            PROB_REINFEC_V2 = params["PROB.REINFEC.V2"],
                            BETA_RATE1 = BETA.RATE1,
                            SEVERITY_FRAC_V2 = SEVERITY.FRAC.V2,
                            ...)
    
    #SOLUTION <- ode(y = new.POP0, times = TIME.VECTOR, func = diffEqs, parms = c(maxsteps = 1000000))
    SOLUTION <- ode(y = new.POP0, times = TIME.VECTOR, func = diffEqs, parms = c())
    SOLUTION <- C_calculator(SOLUTION, EXPOSURE.PERIOD.DAYS.V2, SEVERITY.FRAC.V2)
    
    if (full_solution)
        return(SOLUTION)
    sol <- data.frame(time = SOLUTION[-1,"time"],
                      C1 = diff(SOLUTION[,"POP.C11"] + SOLUTION[,"POP.C12"] + SOLUTION[,"POP.C13"]),
                      C2 = diff(SOLUTION[,"POP.C21"] + SOLUTION[,"POP.C22"] + SOLUTION[,"POP.C23"]),
                      
                      S = diff(SOLUTION[,"POP.S1"] + SOLUTION[,"POP.S2"] + SOLUTION[,"POP.S3"]),
                      
                      E1 = diff(SOLUTION[,"POP.E11"] + SOLUTION[,"POP.E12"] + SOLUTION[,"POP.E13"]),
                      E2 = diff(SOLUTION[,"POP.E21"] + SOLUTION[,"POP.E22"] + SOLUTION[,"POP.E23"]),
                      
                      A1 = diff(SOLUTION[,"POP.A11"] + SOLUTION[,"POP.A12"] + SOLUTION[,"POP.A13"]),
                      A2 = diff(SOLUTION[,"POP.A21"] + SOLUTION[,"POP.A22"] + SOLUTION[,"POP.A23"]),
                      
                      I1 = diff(SOLUTION[,"POP.I11"] + SOLUTION[,"POP.I12"] + SOLUTION[,"POP.I13"]),
                      I2 = diff(SOLUTION[,"POP.I21"] + SOLUTION[,"POP.I22"] + SOLUTION[,"POP.I23"]),
                      
                      H1 = diff(SOLUTION[,"POP.H11"] + SOLUTION[,"POP.H12"] + SOLUTION[,"POP.H13"]),
                      H2 = diff(SOLUTION[,"POP.H21"] + SOLUTION[,"POP.H22"] + SOLUTION[,"POP.H23"]),
                      
                      R1 = diff(SOLUTION[,"POP.R11"] + SOLUTION[,"POP.R12"] + SOLUTION[,"POP.R13"]),
                      R2 = diff(SOLUTION[,"POP.R21"] + SOLUTION[,"POP.R22"] + SOLUTION[,"POP.R23"]),
                      
                      D1 = diff(SOLUTION[,"POP.D11"] + SOLUTION[,"POP.D12"] + SOLUTION[,"POP.D13"]),
                      D2 = diff(SOLUTION[,"POP.D21"] + SOLUTION[,"POP.D22"] + SOLUTION[,"POP.D23"]))
    sol$time <- d0 + (sol$time -1)

    # aggregating data by the frequency time windows
    freq <- rep(NA, nrow(d))
    for (i in 1:nrow(d)) {
        soli <- sol[sol$time >= d[i,"t_start"] & sol$time <= d[i,"t_end"],]
        freq[i] <- sum(soli$C2) / sum(soli$C1, soli$C2)
    }

    # aggregating by epidemiological week
    if(!return_compartments) {
    sol <- sol %>%
        mutate(week = end.of.epiweek(time, end = 0)) %>%
        group_by(week) %>%
        summarise(C1 = sum(C1), C2 = sum(C2)) %>%
        as.data.frame()
    sol.zoo <- zoo(sol[,c("C1","C2")], sol$week) } else {
        
    sol <- sol %>%
            mutate(week = end.of.epiweek(time, end = 0)) %>%
            group_by(week) %>%
            summarise(C1 = sum(C1), C2 = sum(C2),
                      S = sum(S),
                      E1 = sum(E1), E2 = sum(E2),
                      A1 = sum(A1), A2 = sum(A2),
                      I1 = sum(I1), I2 = sum(I2),
                      H1 = sum(H1), H2 = sum(H2),
                      R1 = sum(R1), R2 = sum(R2),
                      D1 = sum(D1), D2 = sum(D2)) %>%
            as.data.frame()
        sol.zoo <- zoo(sol[,c("C1","C2","S","E1","E2","A1","A2","I1","I2","H1","H2","R1","R2","D1","D2")], sol$week)     
    }
    
    # first point is not a full week
    sol.zoo <- sol.zoo[-1]
    
    index(covid.zoo) <- as.Date(index(covid.zoo))
    d.cases <- merge.zoo(sol.zoo, covid.zoo, all = FALSE)
    return(list(cases = d.cases, freq = cbind(d, pred.freq=freq)))
}

#' A function to calculate the residual of the solution
#'
#' @param sol A data.frame of the solution given by model_solution, 'sol'.
#' @param weight.freq Weight to be considered in frequency series 'weight.freq'. Default 1.
#' @param use.logit Logical, define if uses Logit transform 'use.logit'. Default TRUE.
#'
#' @return A vector of Cases and frequencies. 
#' 
residual_sol <- function(sol, weight.freq = 1, use.logit = TRUE) {
    cases <- sol$cases
    freq <- sol$freq$pred.freq
    cases <- (cases - mean(cases$covid.zoo)) / sd(cases$covid.zoo)
    cases <- cases$C1 + cases$C2 - cases$covid.zoo
    if (use.logit) {
        freq.logit <- logit(freq)
        if (any(! is.finite(c(freq.logit, d$freq.logit))))
            warning("Infinite logit value!")
        d.freq <- (freq.logit - d$freq.logit) / sd(d$freq.logit)
    } else {
        d.freq <- (freq - d$freq) / sd(d$freq - mean(d$freq))
    }
    return(c(as.numeric(cases), weight.freq * as.numeric(d.freq)))
}

#' Negative Loglikelihood calculation from the solution model.
#' Negative Log-Likelihood calculation to the a Poisson distribution of cummulative cases solution 
#' and a Binomial distribuition to the new variant frequency.
#'
#' @param sol Data.frame from model_solution 'sol'.
#'
#' @return A value to the Negative Log-Likelihood.
#' 
neg_loglike_sol <- function(sol) {
    cases <- sol$cases
    freq <- sol$freq$pred.freq
    return(-(sum(dpois(cases$covid.zoo, cases$C1 + cases$C2, log = TRUE)) +
             sum(dbinom(x = d$P.1, size = d$total, prob = freq, log = TRUE))))
}

#' Performs the model fitting analysis according to a range of parameters
#'
#' @param lower A vector with the lower limits of parameters that will be fitted
#' @param upper A vector with the upper limits of parameters that will be fitted
#' @param conv.guess A function that converts initial guesses to the proper format (optional)
#' @param LHS Logical, Generates a Latin Hypercube Sampling for initial parameter guesses
#' @param n The number of parameter combinations that will be used for the fitting analysis
#' @param `...` Any other parameters necessary to a different fitting '...'.
#'
#' @return A data.frame with the fitting results for every parameter combination evaluated
#'

all.params.exploration <- function(lower, upper, conv.guess, LHS = TRUE, n = 100, ...) {
    if (LHS) {
        # use LHS to sample the parameter space
        q <- rep('qunif', length(lower))
        q.args <- list()
        for (p in names(lower))
            q.args[[p]] <- list(min=unname(lower[p]), max=unname(upper[p]))
        uncoupledLHS <- LHS(model=NULL, factors=names(lower), N=n,
                            q=q, q.arg=q.args, method='random')
        starters <- uncoupledLHS$data
    } else {
        N <- ceiling(n**(1/length(lower)))
        arg.ranges <- list()
        for (i in 1:length(lower))
            arg.ranges[[i]] <- seq(lower[i], upper[i], length.out = N)
        starters <- do.call(expand.grid, arg.ranges)
        colnames(starters) <- names(lower)
    }

    if (! missing(conv.guess)) {
        starters <- conv.guess(starters)
    }

    ExtraArgs <- list(...)
    if ("IHR.V2.PROP" %in% names(ExtraArgs)) {
        starters$IHR.V2.PROP  <- ExtraArgs[["IHR.V2.PROP"]]
        ExtraArgs["IHR.V2.PROP"] <- NULL
    }
    if ("prevalence" %in% names(ExtraArgs)) {
        starters$prevalence <- ExtraArgs[["prevalence"]]
        ExtraArgs["prevalence"] <- NULL
    }

    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(starters), .combine = rbind) %dopar% {
        sol <- do.call(model_solution, c(list(params=unlist(starters[i,])),
                                         ExtraArgs))
        #residual <- residual_sol(sol)
        neg.loglike <- neg_loglike_sol(sol)

        res <- c(unlist(starters[i,]), #ss = sum(residual**2),
                 neg.loglike = neg.loglike)
    }
    results <- as.data.frame(results)
}

#' Helper function that returns the residual for the nls regression model
#' 
#' @return a vector of residuals

helper.nls.fun <- function(...) {
    return(residual(model_solution(...)))
}

#' Performs in parallel the nonlinear least square (nls) analysis of the model
#' 
#' @param guesses a data.frame with the starting values for the nls method
#' @param lower A vector with the lower limits of fitted parameters for the nls method
#' @param upper A vector with the upper limits of fitted parameters for the nls method
#' @param control a  list of control settings
#' @param ... additional optional arguments
#' 
#' @return a data.frame with the coefficients, residuals, and negative log-likelihood from the analysis

fit.nls.guesses <- function(guesses, lower, upper, control, ...) {
    if (missing(control)) {
        control <- nls.lm.control(ptol=1.490116e-12, ftol=1.490116e-12, maxiter=200,
                          maxfev=4*100*(length(lower) + 1))
    }
    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(guesses), .combine = rbind) %dopar% {
        parms <- unlist(guesses[i,])
        fitval <- nls.lm(par = parms, fn = helper.nls.fun, upper = upper,
                         lower = lower, control = control, ...)
        residuo <- sum(residuals(fitval)*residuals(fitval))
        neg.likelihood <- neg_loglike_sol(model_solution(coef(fitavl)))
        res <- c(coef(fitval), residuo = residuo, neg.loglike = neg.likelihood)
    }
    return(as.data.frame(results))
}

#' Parallel method of the model`s Maximum Likelihood Estimation
#' 
#' @param guesses a data.frame with the parameter values for optimizer
#' @param helper helper function to calculate negative log-likelihood
#' @param ... additional optional arguments
#' 
#' @return a data.frame with the coefficients and log-likelihood for each set of parameters

fit.mle.guesses <- function(guesses, helper, ...) {
    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(guesses), .combine = rbind) %dopar% {
        guess <- guesses[i,]
        fit <- try({mle2(helper, start = as.list(guess), ...)})
        if(class(fit) == "try-error") {
            res <- c(guess, loglike = NA)
        } else {
            res <- c(coef(fit), loglike = logLik(fit))
        }
    }
    return(as.data.frame(results))
}

## Plots of predicted x observed
#' Helper functions to plot solution
#' 
#' @param sol object with solution of the model
#' 
#' @return plot (ggplot format) with the solution of the model


plot.3 <- function(sol){
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2) %>%
        ggplot(aes(x = Index)) +
        geom_point(aes(y=covid.zoo)) +
        geom_line(aes(y = C1), col = "#00B4B8", lty =2) +
        geom_line(aes(y=C2), col = "#F65D57", lty =2) +
        geom_line(aes(y=total)) +
        theme_bw() +
        ylab("Weekly hospitalized cases") +
        xlab("Week")
    p2 <- sol$freq %>%
        fortify() %>%
        filter(source == "fiocruz") %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("P.1 Proportion") +
        xlab("Month")
    p3 <- sol$freq %>%
        fortify() %>%
        filter(source == "cadde") %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("P.1 Proportion") +
        xlab("Dia")
    grid.arrange(p1, p2, p3, ncol =3, nrow =1)
}

plot.2 <- function(sol){
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2) %>%
        ggplot(aes(x = Index)) +
        geom_point(aes(y=covid.zoo)) +
        geom_line(aes(y = C1), col = "#00B4B8", lty =2) +
        geom_line(aes(y=C2), col = "#F65D57", lty =2) +
        geom_line(aes(y=total)) +
        theme_bw() +
        ylab("Weekly hospitalized cases") +
        xlab("Week")
    p2 <- sol$freq %>%
        fortify() %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("P.1 Proportion") +
        xlab("Week")
    grid.arrange(p1, p2, ncol = 2, nrow = 1)
}

# plot with confidence interval
# to produce plots for the paper, run before:
# Sys.setlocale("LC_ALL","en_US") # labels in english
# and save it using:
# dev.copy(png, "plot.CI.fit3.f2.png", width=1000, height=1800, res=250, pointsize=10)
# OR in EPS (open BEFORE plotting):
# cairo_ps("plot.CI.fit3.f2.eps", width=6, height=8, pointsize=12)
# dev.off()

plot.3.CI <- function(sol){
    tmin <- min(sol$freq$t_start[1], time(sol$cases)[1])
    tmax <- max(sol$freq$t_end[length(sol$freq$t_end)], max(time(sol$cases)))
    sol$freq$t_plot <- sol$freq$t_start + difftime(sol$freq$t_end, sol$freq$t_start) / 2
    color1 <- "blue" # #00B4B8
    color2 <- "red" # #F65D57
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2,
               total.inf = C1.q2.5 + C2.q2.5,
               total.sup = C1.q97.5 + C2.q97.5) %>%
        ggplot(aes(x = Index)) +
        geom_line(aes(y = C1), col = color1, lty =2) +
        geom_line(aes(y=C2), col = color2, lty =2) +
        geom_line(aes(y=total)) +
        geom_ribbon(aes(ymin = C1.q2.5, ymax = C1.q97.5), fill = "blue", alpha = 0.3) +
        geom_ribbon(aes(ymin = C2.q2.5, ymax = C2.q97.5), fill = "red", alpha = 0.3) +
        geom_ribbon(aes(ymin = total.inf, ymax = total.sup), fill = "black", alpha = 0.3) +
        geom_point(aes(y=covid.zoo)) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("New hospitalized cases") +
        xlab("")
    p2 <- sol$freq %>%
        fortify() %>%
        filter(source == "fiocruz") %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("")
    p3 <- sol$freq %>%
        fortify() %>%
        filter(source == "cadde") %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("date")
    grid.arrange(p1, p2, p3, ncol =1, nrow =3)
}

plot.2.CI <- function(sol){
    tmin <- min(sol$freq$t_start[1], time(sol$cases)[1])
    tmax <- max(sol$freq$t_end[length(sol$freq$t_end)], max(time(sol$cases)))
    sol$freq$t_plot <- sol$freq$t_start + difftime(sol$freq$t_end, sol$freq$t_start) / 2
    color1 <- "blue" # #00B4B8
    color2 <- "red" # #F65D57
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2,
               total.inf = C1.q2.5 + C2.q2.5,
               total.sup = C1.q97.5 + C2.q97.5) %>%
        ggplot(aes(x = Index)) +
        geom_line(aes(y = C1), col = color1, lty =2) +
        geom_line(aes(y=C2), col = color2, lty =2) +
        geom_line(aes(y=total)) +
        geom_ribbon(aes(ymin = C1.q2.5, ymax = C1.q97.5), fill = "blue", alpha = 0.3) +
        geom_ribbon(aes(ymin = C2.q2.5, ymax = C2.q97.5), fill = "red", alpha = 0.3) +
        geom_ribbon(aes(ymin = total.inf, ymax = total.sup), fill = "black", alpha = 0.3) +
        geom_point(aes(y=covid.zoo)) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("New hospitalized cases") +
        xlab("")
    p2 <- sol$freq %>%
        fortify() %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("date")
    grid.arrange(p1, p2, ncol = 1, nrow = 2)
}


## Ad hoc helper functions to convert back and forth parameters to logit scale ##
## Convert guesses that are probabilities to logit scale
prep.guess <- function(guess){
    if ("prevalence" %in% names(guess))
        guess["prevalence"] <- logit(guess["prevalence"]/100)
    guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- logit(guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(guess)
}
## convert guess back to prob. scale
prep.guess.inv <- function(guess) {
    if ("prevalence" %in% names(guess))
        guess["prevalence"] <- 100 * inv.logit(guess["prevalence"])
    guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- inv.logit(guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(guess)
}

## Convert coefficients that are at the logit scale back to
## proportions
prep.coef <- function(fit){
    fit.cf <- coef(fit)
    fit.parms <- fit.cf
    if ("prevalence" %in% names(fit.parms))
        fit.parms["prevalence"] <- inv.logit(fit.parms["prevalence"])*100
    if("IHR.V2.PROP" %in% names(fit.parms))
        fit.parms["IHR.V2.PROP"] <- exp(fit.parms["IHR.V2.PROP"])
    fit.parms[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- inv.logit(fit.parms[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(fit.parms)
}
## Convert confidence intervals that are at the logit scale
## back to proportions
prep.conf <- function(fit, ...){
    fit.conf <- confint(fit, ...)
    fit.conf[c("INIT.V2.FRAC", "PROB.REINFEC.V2"),] <- inv.logit(fit.conf[c("INIT.V2.FRAC", "PROB.REINFEC.V2"),])
    if ("prevalence" %in% names(fit.conf))
        fit.conf["prevalence",] <- inv.logit(fit.conf["prevalence",])*100
    return(fit.conf)
}
