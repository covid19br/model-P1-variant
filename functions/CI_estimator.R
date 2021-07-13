source("./functions/logit.R")

if (!require(doParallel)) install.packages(doParallel); library(doParallel)
if (!require(mvtnorm)) install.packages(mvtnorm); library(mvtnorm)

solution_CI <- function(fit, n = 1000, conv.logit=c(), extra_parms) {
    ###### SAMPLE GENERATION ##########
    coefs <- coef(fit)
    sigma <- vcov(fit)
    samples <- rmvnorm(n, mean = coefs, sigma = sigma, method = "svd")
    if (length(conv.logit) > 0) {
        samples[, conv.logit] <- inv.logit(samples[, conv.logit]) ###inverse logit
        coefs[conv.logit] <- inv.logit(coefs[conv.logit])
        if ("prevalence" %in% names(coefs)) {
            samples[, "prevalence"] <- 100 * samples[, "prevalence"]
            coefs["prevalence"] <- 100 * coefs["prevalence"]
        }
    }
    if (! missing(extra_parms)) {
        samples <- cbind(samples, as.data.frame(as.list(extra_parms)))
        coefs <- c(coefs, extra_parms)
    }

    # run simulations
    registerDoParallel(cores=detectCores()-1)
    sols <- foreach (i = 1:dim(samples)[1], .combine = rbind) %dopar% {
        sol <- model_solution(unlist(samples[i,]))
        c(as.numeric(sol$cases$C1), as.numeric(sol$cases$C2), sol$freq$pred.freq)
    }
    aggregated <- data.frame(mean = apply(sols, 2, mean),
                             sd = apply(sols, 2, sd),
                             q2.5 = apply(sols, 2, quantile, 0.025),
                             q50 = apply(sols, 2, quantile, 0.5),
                             q97.5 = apply(sols, 2, quantile, 0.975))
    # this is ugly, but how else?
    sol1 <- model_solution(coefs)
    t1 <- length(time(sol1$cases))
    # cases
    sol1$cases$C1.mean <- aggregated$mean[1:t1]
    sol1$cases$C2.mean <- aggregated$mean[(t1+1):(2*t1)]
    sol1$cases$C1.sd <- aggregated$sd[1:t1]
    sol1$cases$C2.sd <- aggregated$sd[(t1+1):(2*t1)]
    sol1$cases$C1.q2.5 <- aggregated$q2.5[1:t1]
    sol1$cases$C2.q2.5 <- aggregated$q2.5[(t1+1):(2*t1)]
    sol1$cases$C1.q50 <- aggregated$q50[1:t1]
    sol1$cases$C2.q50 <- aggregated$q50[(t1+1):(2*t1)]
    sol1$cases$C1.q97.5 <- aggregated$q97.5[1:t1]
    sol1$cases$C2.q97.5 <- aggregated$q97.5[(t1+1):(2*t1)]
    # freq
    sol1$freq$pred.freq.mean <- aggregated$mean[(2*t1+1):length(aggregated$mean)]
    sol1$freq$q2.5 <- aggregated$q2.5[(2*t1+1):length(aggregated$mean)]
    sol1$freq$q50 <- aggregated$q50[(2*t1+1):length(aggregated$mean)]
    sol1$freq$q97.5 <- aggregated$q97.5[(2*t1+1):length(aggregated$mean)]
    # freq error bars
    error_freq <- binconf(sol1$freq$P.1,sol1$freq$total,method = "wilson")
    sol1$freq$obs_freq_lower <- error_freq[,"Lower"]
    sol1$freq$obs_freq_upper <- error_freq[,"Upper"]
    return(sol1)
}

